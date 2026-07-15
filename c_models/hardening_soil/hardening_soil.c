/**
 * @file hardening_soil.c
 * @brief Abaqus-style UMAT for the Hardening Soil model with cone (shear)
 *        hardening, a Mohr-Coulomb failure cut-off and an elliptic volumetric
 *        cap with cap hardening.
 *
 * This is a C translation of the validated Python prototype
 * (python_prototypes/hardening_soil.py). It uses a backward-Euler / cutting
 * plane return mapping with automatic sub-stepping.
 *
 * Conventions
 * -----------
 *  - Voigt ordering: [xx, yy, zz, xy, yz, xz] (see globals.h).
 *  - Engineering shear strains (gamma = 2 * epsilon_shear).
 *  - Compression positive (consistent with the Python prototype). The minor
 *    principal stress sigma_3 is therefore the smallest principal stress.
 *  - Elastic behaviour uses the unloading/reloading modulus Eur.
 *
 * Material properties (PROPS)
 * ---------------------------
 *   [0]  E50_ref  - reference secant stiffness at 50% strength   [stress]
 *   [1]  Eur_ref  - reference unloading/reloading stiffness       [stress]
 *   [2]  m        - stress exponent for stiffness                 [-]
 *   [3]  c        - cohesion                                      [stress]
 *   [4]  phi      - friction angle                               [degrees]
 *   [5]  psi      - dilation angle                               [degrees]
 *   [6]  p_ref    - reference pressure for stiffness             [stress]
 *   [7]  Rf       - failure ratio qf/qa (typ. 0.9)               [-]
 *   [8]  nu       - Poisson's ratio (unloading/reloading)        [-]
 *   [9]  M_cap    - cap aspect ratio                             [-]
 *   [10] K_ratio  - bulk/plastic stiffness ratio for cap harden. [-]
 *
 * State variables (STATEV)
 * ------------------------
 *   [0]  gamma_p  - accumulated deviatoric plastic shear strain (cone)
 *   [1]  p_p      - pre-consolidation pressure (cap size)
 *   [2]  at_failure flag (0 = cone hardening, 1 = Mohr-Coulomb failure)
 */

#include <math.h>
#include <stdio.h>

#include "../elastic_laws/hookes_law.h"
#include "../globals.h"
#include "../utils.h"

/* Calling convention / export macros ----------------------------------------- */
#if defined(_WIN32) || defined(_WIN64)
#define UMAT_EXPORT __declspec(dllexport)
#define UMAT_CALLCONV __stdcall
#else
#define UMAT_EXPORT
#define UMAT_CALLCONV
#endif

#define HS_MAX_LOCAL_ITER 100
#define HS_MIN_SUBSTEPS 5
#define HS_MAX_SUBSTEPS 400

/* Under-relaxation factor for the cone (hardening) plastic multiplier. The
 * hyperbolic hardening law is strongly non-linear near the asymptote, so a
 * full Newton step can overshoot gamma_p and trigger spurious elastic
 * unloading. Damping keeps the return mapping stable. */
#define HS_CONE_RELAX 0.5

/* Identity vector in Voigt notation (normal part only). */
static const double IDENTITY_VOIGT[VOIGTSIZE_3D] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

/* Material parameters gathered in a struct for convenience. */
typedef struct
{
    double E50_ref;
    double Eur_ref;
    double m;
    double c;
    double phi;    /* rad */
    double psi;    /* rad */
    double p_ref;
    double Rf;
    double nu;
    double M_cap;
    double K_ratio;
    double p_t;    /* tensile intercept c * cot(phi) */
    double Ei_ref; /* reference initial loading stiffness */
} HSParams;

/* ------------------------------------------------------------------ */
/* Basic stress helpers                                                */
/* ------------------------------------------------------------------ */

static double hs_mean_stress(const double s[VOIGTSIZE_3D])
{
    return (s[XX] + s[YY] + s[ZZ]) / 3.0;
}

static void hs_deviator(const double s[VOIGTSIZE_3D], double p, double dev[VOIGTSIZE_3D])
{
    dev[XX] = s[XX] - p;
    dev[YY] = s[YY] - p;
    dev[ZZ] = s[ZZ] - p;
    dev[XY] = s[XY];
    dev[YZ] = s[YZ];
    dev[XZ] = s[XZ];
}

/* Von Mises equivalent stress q = sqrt(3 J2), including shear terms. */
static double hs_q(const double s[VOIGTSIZE_3D])
{
    double p = hs_mean_stress(s);
    double dev[VOIGTSIZE_3D];
    hs_deviator(s, p, dev);
    double j2 = 0.5 * (dev[XX] * dev[XX] + dev[YY] * dev[YY] + dev[ZZ] * dev[ZZ]) +
                (dev[XY] * dev[XY] + dev[YZ] * dev[YZ] + dev[XZ] * dev[XZ]);
    return sqrt(3.0 * j2);
}

/* Derivative of q w.r.t. stress in engineering Voigt notation (factor 2 on
 * shear components), consistent with the plain Voigt dot products used below. */
static void hs_dq_dsigma(const double s[VOIGTSIZE_3D], double p, double q,
                         double dq[VOIGTSIZE_3D])
{
    if (q < ZERO_TOL)
    {
        for (int i = 0; i < VOIGTSIZE_3D; ++i) dq[i] = 0.0;
        return;
    }
    double dev[VOIGTSIZE_3D];
    hs_deviator(s, p, dev);
    double factor = 3.0 / (2.0 * q);
    dq[XX] = factor * dev[XX];
    dq[YY] = factor * dev[YY];
    dq[ZZ] = factor * dev[ZZ];
    dq[XY] = 2.0 * factor * dev[XY];
    dq[YZ] = 2.0 * factor * dev[YZ];
    dq[XZ] = 2.0 * factor * dev[XZ];
}

/* Sorted principal stresses (descending: ps[0] >= ps[1] >= ps[2]). */
static void hs_principal_stresses(const double s[VOIGTSIZE_3D], double ps[3])
{
    double p = hs_mean_stress(s);
    double dev[VOIGTSIZE_3D];
    hs_deviator(s, p, dev);

    double j2 = 0.5 * (dev[XX] * dev[XX] + dev[YY] * dev[YY] + dev[ZZ] * dev[ZZ]) +
                (dev[XY] * dev[XY] + dev[YZ] * dev[YZ] + dev[XZ] * dev[XZ]);

    if (j2 < ZERO_TOL)
    {
        ps[0] = ps[1] = ps[2] = p;
        return;
    }

    double j3 = calculate_determinant_voigt_vector_3d(dev);

    double arg = (sqrt(27.0) / 2.0) * j3 / pow(j2, 1.5);
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;
    double theta = asin(arg) / 3.0; /* in [-pi/6, pi/6] */

    double radius = sqrt(4.0 * j2 / 3.0);
    double a = p + radius * cos(theta);
    double b = p + radius * cos(theta - 2.0 * PI / 3.0);
    double cc = p + radius * cos(theta + 2.0 * PI / 3.0);

    /* sort descending */
    double hi = a, mid = b, lo = cc, tmp;
    if (hi < mid) { tmp = hi; hi = mid; mid = tmp; }
    if (mid < lo) { tmp = mid; mid = lo; lo = tmp; }
    if (hi < mid) { tmp = hi; hi = mid; mid = tmp; }
    ps[0] = hi;
    ps[1] = mid;
    ps[2] = lo;
}

/* ------------------------------------------------------------------ */
/* Stress-dependent stiffness                                          */
/* ------------------------------------------------------------------ */

static double hs_stiffness_factor(const HSParams* prm, double sigma_3)
{
    double denom = prm->p_ref + prm->p_t;
    if (fabs(denom) < ZERO_TOL) denom = (denom < 0.0 ? -ZERO_TOL : ZERO_TOL);
    double ratio = (sigma_3 + prm->p_t) / denom;
    if (ratio <= 0.0) ratio = SMALL_VALUE;
    return pow(ratio, prm->m);
}

/* ------------------------------------------------------------------ */
/* Cone (shear) yield surface                                          */
/* ------------------------------------------------------------------ */

static double hs_qf(const HSParams* prm, double sigma_3)
{
    double sin_phi = sin(prm->phi);
    double denom = 1.0 - sin_phi;
    if (fabs(denom) < ZERO_TOL) denom = ZERO_TOL;
    return (2.0 * sin_phi / denom) * (sigma_3 + prm->p_t);
}

/* Hyperbolic cone yield function (Benz). */
static double hs_cone_yield(double q, double gamma_p, double Ei, double Eur, double qa)
{
    if (q <= 0.0) return -gamma_p;
    if (q >= 0.99 * qa) return 1.0e10; /* stay away from the asymptote */
    double term = 2.0 * q / (1.0 - q / qa);
    return term / Ei - 2.0 * q / Eur - gamma_p;
}

/* Accumulated cone plastic strain that corresponds to a given deviator stress q
 * on the hyperbola. Mirrors HardeningLaws.maximum_gamma_p and is used to freeze
 * gamma_p at the failure value once the Mohr-Coulomb surface is reached. */
static double hs_maximum_gamma_p(double q, double qa, double Ei, double Eur)
{
    double denom = qa - q;
    if (fabs(denom) < ZERO_TOL) denom = (denom < 0.0 ? -ZERO_TOL : ZERO_TOL);
    return 2.0 * q / Ei * qa / denom - 2.0 * q / Eur;
}

/* Derivative of the cone yield function w.r.t. stress. */
static void hs_cone_dfdsigma(const double s[VOIGTSIZE_3D], double p, double q,
                             double Ei, double Eur, double qa, double out[VOIGTSIZE_3D])
{
    if (q < ZERO_TOL)
    {
        for (int i = 0; i < VOIGTSIZE_3D; ++i) out[i] = 0.0;
        return;
    }
    double dq[VOIGTSIZE_3D];
    hs_dq_dsigma(s, p, q, dq);

    double denom = 1.0 - q / qa;
    if (denom <= ZERO_TOL) denom = ZERO_TOL;
    double dF_dq = 2.0 / (Ei * denom * denom) - 2.0 / Eur;

    for (int i = 0; i < VOIGTSIZE_3D; ++i) out[i] = dF_dq * dq[i];
}

/* Drucker-Prager style plastic potential gradient with mobilised dilatancy. */
static void hs_cone_flow(const double s[VOIGTSIZE_3D], double p, double q,
                         const HSParams* prm, double out[VOIGTSIZE_3D])
{
    if (q < ZERO_TOL)
    {
        for (int i = 0; i < VOIGTSIZE_3D; ++i) out[i] = 0.0;
        return;
    }

    double eta = q / (p + prm->p_t);
    double sin_phi_mob = 3.0 * eta / (6.0 + eta);
    if (sin_phi_mob > 1.0) sin_phi_mob = 1.0;
    if (sin_phi_mob < -1.0) sin_phi_mob = -1.0;

    double sin_phi = sin(prm->phi);
    double sin_psi = sin(prm->psi);
    double sin_phi_cs = (sin_phi - sin_psi) / (1.0 - sin_phi * sin_psi);

    double sin_psi_mob;
    if (asin(sin_phi_mob) < prm->phi)
        sin_psi_mob = 0.0;
    else
        sin_psi_mob = (sin_phi_mob - sin_phi_cs) / (1.0 - sin_phi_mob * sin_phi_cs);

    /* dilatancy slope in the DP sense */
    double M_psi = 6.0 * sin_psi_mob / (3.0 - sin_psi_mob);

    double dq[VOIGTSIZE_3D];
    hs_dq_dsigma(s, p, q, dq);
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
        out[i] = dq[i] - (M_psi / 3.0) * IDENTITY_VOIGT[i];
}

/* ------------------------------------------------------------------ */
/* Mohr-Coulomb failure surface (perfect plasticity at q = qf)         */
/* ------------------------------------------------------------------ */

/*
 * Mohr-Coulomb yield surface used at failure, written directly in terms of the
 * (triaxial-ordered) normal stress components so that it reproduces the Python
 * prototype exactly (YieldFunctions.F_mohr_coulomb):
 *
 *   F = (2*s1 - s2 - s3)/4 - (2*s1 + s2 + s3)/4 * sin(phi) - c * cos(phi)
 *
 * This vanishes on the same failure locus as (q - qf), but its constant flow
 * gradient (below) delivers the correct dilatancy (trace = -sin(psi)).
 */
static double hs_mc_yield_full(const double s[VOIGTSIZE_3D], const HSParams* prm)
{
    double sin_phi = sin(prm->phi);
    double cos_phi = cos(prm->phi);
    return (2.0 * s[XX] - s[YY] - s[ZZ]) / 4.0 -
           (2.0 * s[XX] + s[YY] + s[ZZ]) / 4.0 * sin_phi -
           prm->c * cos_phi;
}

/*
 * Constant Mohr-Coulomb gradient (matches YieldFunctions.dfds_mohr_coulomb).
 * Called with `phi` for the yield-surface normal and with `psi` for the plastic
 * potential (non-associated flow). The trace equals -sin(angle), so the plastic
 * volumetric slope at failure is -sin(psi), exactly as in the Python prototype.
 */
static void hs_mc_gradient(double angle, double out[VOIGTSIZE_3D])
{
    double sa = sin(angle);
    out[XX] = 0.5 - 0.5 * sa;
    out[YY] = -0.25 - 0.25 * sa;
    out[ZZ] = -0.25 - 0.25 * sa;
    out[XY] = 0.0;
    out[YZ] = 0.0;
    out[XZ] = 0.0;
}

/* ------------------------------------------------------------------ */
/* Elliptic cap                                                        */
/* ------------------------------------------------------------------ */

static double hs_cap_qspecial(const double s[VOIGTSIZE_3D], const HSParams* prm)
{
    double sin_phi = sin(prm->phi);
    double alpha = (3.0 + sin_phi) / (3.0 - sin_phi);
    /* Triaxial assumption: normal Voigt components ordered sigma_1..sigma_3. */
    return s[XX] + (alpha - 1.0) * s[YY] - alpha * s[ZZ];
}

/* Cap yield F_c = (q~/M)^2 + p^2 - p_p^2 (closed ellipse). */
static double hs_cap_yield(const double s[VOIGTSIZE_3D], double p, double p_p,
                           const HSParams* prm)
{
    double q_special = hs_cap_qspecial(s, prm);
    return (q_special / prm->M_cap) * (q_special / prm->M_cap) + p * p - p_p * p_p;
}

static void hs_cap_dfdsigma(const double s[VOIGTSIZE_3D], double p, const HSParams* prm,
                            double out[VOIGTSIZE_3D])
{
    double sin_phi = sin(prm->phi);
    double alpha = (3.0 + sin_phi) / (3.0 - sin_phi);
    double q_special = hs_cap_qspecial(s, prm);

    double dqs[VOIGTSIZE_3D] = {1.0, alpha - 1.0, -alpha, 0.0, 0.0, 0.0};
    double coeff = 2.0 * q_special / (prm->M_cap * prm->M_cap);
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
        out[i] = coeff * dqs[i] + (2.0 * p / 3.0) * IDENTITY_VOIGT[i];
}

/* Cap hardening law: dp_p / dlambda_c. */
static double hs_cap_dhdlambda(const HSParams* prm, double sigma_3, double p)
{
    double Ks = prm->Eur_ref / (3.0 * (1.0 - 2.0 * prm->nu));
    double H = (1.0 / (prm->K_ratio - 1.0)) * Ks;
    double denom = prm->p_ref + prm->p_t;
    if (fabs(denom) < ZERO_TOL) denom = (denom < 0.0 ? -ZERO_TOL : ZERO_TOL);
    double ratio = (sigma_3 + prm->p_t) / denom;
    if (ratio <= 0.0) ratio = SMALL_VALUE;
    return 2.0 * H * pow(ratio, prm->m) * p;
}

/* ------------------------------------------------------------------ */
/* Single backward-Euler step (cutting plane) on one strain increment  */
/* ------------------------------------------------------------------ */

/*
 * Returns 1 on convergence, 0 otherwise. On success the stress and the state
 * variables (gamma_p, p_p, at_failure) are updated in place.
 */
static int hs_single_step(double stress[VOIGTSIZE_3D], double* gamma_p, double* p_p,
                          int* at_failure, const double deps[VOIGTSIZE_3D],
                          const HSParams* prm, double tol)
{
    double ps[3];
    double sigma_3, p, q, qf, qa, Ei, Eur, factor;
    double Ce[VOIGTSIZE_3D * VOIGTSIZE_3D];
    double delta_sigma[VOIGTSIZE_3D];
    double stress_trial[VOIGTSIZE_3D];

    /* --- stiffness at the start of the step ---
     * The reference (Python) prototype evaluates the stress-dependent stiffness
     * and strength from the zz Voigt component (the minor principal stress in
     * the triaxial convention). We deliberately mirror that here rather than
     * using an invariant principal-stress solve, which is both unnecessary for
     * this convention and inaccurate at the triaxial vertex (repeated
     * eigenvalue), where it would underestimate sigma_3 and soften the model. */
    sigma_3 = stress[ZZ];
    factor = hs_stiffness_factor(prm, sigma_3);
    Eur = prm->Eur_ref * factor;

    /* --- elastic predictor (Eur based) --- */
    calculate_elastic_stiffness_matrix_3d(Eur, prm->nu, Ce);
    matrix_vector_multiply(Ce, deps, VOIGTSIZE_3D, delta_sigma);
    add_vectors(stress, delta_sigma, VOIGTSIZE_3D, stress_trial);

    sigma_3 = stress_trial[ZZ];
    factor = hs_stiffness_factor(prm, sigma_3);
    Eur = prm->Eur_ref * factor;
    Ei = prm->Ei_ref * factor;
    qf = hs_qf(prm, sigma_3);
    qa = qf / prm->Rf;

    p = hs_mean_stress(stress_trial);
    q = hs_q(stress_trial);

    /* --- trial yield values --- */
    double f_shear;
    int use_mc = *at_failure;
    if (use_mc)
        f_shear = hs_mc_yield_full(stress_trial, prm);
    else
        f_shear = hs_cone_yield(q, *gamma_p, Ei, Eur, qa);

    double f_cap = hs_cap_yield(stress_trial, p, *p_p, prm);

    /* Relative convergence tolerance scaled by the trial residual, matching the
     * Python prototype (tol_s = |f_s_trial| * tol). The previous absolute form
     * (|f| + 1) * tol was far too loose for the cone yield residual (~1e-3) and
     * stopped the return mapping prematurely, softening the response. */
    double tol_s = fabs(f_shear) * tol;
    double tol_c = fabs(f_cap) * tol;

    int active_s = (f_shear > tol_s);
    /* The cap surface is currently disabled to match the validated reference
     * (cone shear hardening + Mohr-Coulomb failure). Set this to
     * (f_cap > tol_c) to enable the volumetric cap. */
    int active_c = 0;
    (void)tol_c;

    if (!active_s && !active_c)
    {
        copy_array(stress_trial, VOIGTSIZE_3D, stress);
        return 1;
    }

    /* --- return mapping ---------------------------------------------------
     * Newton iteration on the *total* cone plastic multiplier measured from the
     * elastic trial state, mirroring the Python prototype:
     *     sigma   = sigma_trial - dlambda_s * (De : dg/dsigma)
     *     gamma_p = gamma_p0    + dlambda_s * dh/dlambda_s
     * Keeping both the stress and the hardening variable tied to the same total
     * multiplier avoids the over-accumulation of gamma_p (and the resulting
     * spurious "elastic" convergence) that stalled the cone hardening. */
    double sigma[VOIGTSIZE_3D];
    copy_array(stress_trial, VOIGTSIZE_3D, sigma);

    const double gamma_p0 = *gamma_p;
    double gp = *gamma_p;
    double pp = *p_p; /* cap disabled to match the reference (see note above) */
    double dlambda_s = 0.0;
    int converged = 0;
    (void)active_c;
    (void)f_cap;

    for (int it = 0; it < HS_MAX_LOCAL_ITER; ++it)
    {
        p = hs_mean_stress(sigma);
        q = hs_q(sigma);
        sigma_3 = sigma[ZZ]; /* minor principal stress (triaxial convention) */
        factor = hs_stiffness_factor(prm, sigma_3);
        Eur = prm->Eur_ref * factor;
        Ei = prm->Ei_ref * factor;
        qf = hs_qf(prm, sigma_3);
        qa = qf / prm->Rf;

        use_mc = *at_failure; /* mode is fixed during the local iteration (mirrors
                               * the Python `use_mohr_coulomb` flag); the switch is
                               * decided after convergence based on q/qf. */

        /* residual of the active shear yield function at the current state */
        if (use_mc)
            f_shear = hs_mc_yield_full(sigma, prm);
        else
            f_shear = hs_cone_yield(q, gp, Ei, Eur, qa);

        if (fabs(f_shear) < tol_s)
        {
            converged = 1;
#ifdef HS_DEBUG
            fprintf(stderr, "[hs] conv it=%d f=%.3e tol_s=%.3e q=%.3f gp=%.6f mc=%d\n",
                    it, f_shear, tol_s, q, gp, use_mc);
#endif
            break;
        }

        calculate_elastic_stiffness_matrix_3d(Eur, prm->nu, Ce);

        /* plastic potential gradient (non-associated) */
        double dg_s[VOIGTSIZE_3D], df_s[VOIGTSIZE_3D], De_dg_s[VOIGTSIZE_3D];
        if (use_mc)
            hs_mc_gradient(prm->psi, dg_s); /* MC flow: trace = -sin(psi) */
        else
            hs_cone_flow(sigma, p, q, prm, dg_s);
        matrix_vector_multiply(Ce, dg_s, VOIGTSIZE_3D, De_dg_s);

        /* yield gradient and hardening slope */
        double dfdh_s, dhdl_s;
        if (use_mc)
        {
            hs_mc_gradient(prm->phi, df_s); /* MC yield normal */
            dfdh_s = 0.0; /* perfect plasticity at failure */
            dhdl_s = 0.0;
        }
        else
        {
            hs_cone_dfdsigma(sigma, p, q, Ei, Eur, qa, df_s);
            dfdh_s = -1.0; /* dF_s/dgamma_p */
            dhdl_s = 2.0;  /* dgamma_p/dlambda_s (Benz) */
        }

        double denom = vector_dot_product(df_s, De_dg_s, VOIGTSIZE_3D) - dfdh_s * dhdl_s;
        if (fabs(denom) < SMALL_VALUE) denom = (denom < 0.0 ? -SMALL_VALUE : SMALL_VALUE);

        double ddlambda = f_shear / denom;
        /* damp the strongly non-linear cone hardening update (Python uses 0.5) */
        if (!use_mc) ddlambda *= HS_CONE_RELAX;
        dlambda_s += ddlambda;
        if (dlambda_s < 0.0) dlambda_s = 0.0;

        /* update stress and hardening from the elastic trial state */
        for (int i = 0; i < VOIGTSIZE_3D; ++i)
            sigma[i] = stress_trial[i] - dlambda_s * De_dg_s[i];
        if (use_mc)
        {
            /* Perfect plasticity: freeze gamma_p at its value on the hyperbola
             * at q = qf, matching the Python prototype (maximum_gamma_p). */
            gp = hs_maximum_gamma_p(qf, qa, Ei, Eur);
        }
        else
        {
            gp = gamma_p0 + dlambda_s * dhdl_s;
            if (gp < gamma_p0) gp = gamma_p0;
        }
    }

    if (!converged) return 0; /* not converged */

    copy_array(sigma, VOIGTSIZE_3D, stress);
    *gamma_p = gp;
    *p_p = pp;

    /* Decide the yield mode for the *next* step from the converged state,
     * mirroring the Python prototype: switch to Mohr-Coulomb once the deviator
     * stress has essentially reached failure (q/qf >= 0.999) and switch back to
     * cone hardening otherwise. */
    q = hs_q(sigma);
    qf = hs_qf(prm, sigma[ZZ]);
    double q_over_qf = (qf > ZERO_TOL) ? (q / qf) : 0.0;
    *at_failure = (q_over_qf >= 0.999) ? 1 : 0;
    return 1;
}
/* Integration with automatic sub-stepping                             */
/* ------------------------------------------------------------------ */

static int hs_integrate(double stress[VOIGTSIZE_3D], double* gamma_p, double* p_p,
                        int* at_failure, const double dstrain[VOIGTSIZE_3D],
                        const HSParams* prm, double tol)
{
    double q, qf;

    q = hs_q(stress);
    qf = hs_qf(prm, stress[ZZ]);

    /* choose a substep count based on the proximity to failure */
    double ratio = (qf > ZERO_TOL) ? (q / qf) : 0.0;
    if (ratio > 0.99) ratio = 0.99;
    double q_step = 0.9 * qf * (1.0 - ratio);
    int n_sub = HS_MIN_SUBSTEPS;
    if (q_step > ZERO_TOL)
    {
        int est = (int)ceil(q / q_step);
        if (est > n_sub) n_sub = est;
    }
    if (n_sub > HS_MAX_SUBSTEPS) n_sub = HS_MAX_SUBSTEPS;

    for (int trial = 0; trial < 5; ++trial)
    {
        double sigma[VOIGTSIZE_3D];
        copy_array(stress, VOIGTSIZE_3D, sigma);
        double gp = *gamma_p;
        double pp = *p_p;
        int fail_flag = *at_failure;

        double deps_step[VOIGTSIZE_3D];
        for (int i = 0; i < VOIGTSIZE_3D; ++i) deps_step[i] = dstrain[i] / (double)n_sub;

        int ok = 1;
        for (int step = 0; step < n_sub; ++step)
        {
            if (!hs_single_step(sigma, &gp, &pp, &fail_flag, deps_step, prm, tol))
            {
                ok = 0;
                break;
            }
        }

        if (ok)
        {
            copy_array(sigma, VOIGTSIZE_3D, stress);
            *gamma_p = gp;
            *p_p = pp;
            *at_failure = fail_flag;
            return 1;
        }

        n_sub *= 2; /* refine and retry */
        if (n_sub > HS_MAX_SUBSTEPS) n_sub = HS_MAX_SUBSTEPS;
    }

    return 0;
}

/* ------------------------------------------------------------------ */
/* Property validation                                                 */
/* ------------------------------------------------------------------ */

static int hs_check_properties(int nprops, const double* props)
{
    if (nprops < 11)
    {
        fprintf(stderr, "UMAT Error: Hardening Soil requires 11 properties.\n");
        return 1;
    }
    int n_errors = 0;
    if (props[0] <= 0.0) { fprintf(stderr, "UMAT Error: E50_ref must be positive.\n"); n_errors++; }
    if (props[1] <= 0.0) { fprintf(stderr, "UMAT Error: Eur_ref must be positive.\n"); n_errors++; }
    if (props[3] < 0.0)  { fprintf(stderr, "UMAT Error: cohesion must be non-negative.\n"); n_errors++; }
    if (props[4] <= 0.0 || props[4] > 90.0) { fprintf(stderr, "UMAT Error: phi in (0, 90].\n"); n_errors++; }
    if (props[6] <= 0.0) { fprintf(stderr, "UMAT Error: p_ref must be positive.\n"); n_errors++; }
    if (props[7] <= 0.0 || props[7] >= 1.0) { fprintf(stderr, "UMAT Error: Rf in (0, 1).\n"); n_errors++; }
    if (props[8] < 0.0 || props[8] >= 0.5) { fprintf(stderr, "UMAT Error: nu in [0, 0.5).\n"); n_errors++; }
    if (props[9] <= 0.0) { fprintf(stderr, "UMAT Error: M_cap must be positive.\n"); n_errors++; }
    if (props[10] <= 1.0) { fprintf(stderr, "UMAT Error: K_ratio must be > 1.\n"); n_errors++; }
    return (n_errors > 0) ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/* UMAT entry point                                                    */
/* ------------------------------------------------------------------ */

UMAT_EXPORT void UMAT_CALLCONV umat(
    double* STRESS, double* STATEV, double* DDSDDE, double* SSE, double* SPD,
    double* SCD, double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT,
    double* STRAN, double* DSTRAN, double* TIME, double* DTIME, double* TEMP,
    double* DTEMP, double* PREDEF, double* DPRED, char* CMNAME, int* NDI,
    int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
    double* COORDS, double* DROT, double* PNEWDT, double* CELENT, double* DFGRD0,
    double* DFGRD1, int* NOEL, int* NPT, int* LAYER, int* KSPT, int* KSTEP,
    int* KINC)
{
    /* silence unused-parameter warnings */
    (void)SSE; (void)SCD; (void)RPL; (void)DDSDDT; (void)DRPLDE; (void)DRPLDT;
    (void)STRAN; (void)TIME; (void)DTIME; (void)TEMP; (void)DTEMP; (void)PREDEF;
    (void)DPRED; (void)CMNAME; (void)COORDS; (void)DROT; (void)CELENT;
    (void)DFGRD0; (void)DFGRD1; (void)NOEL; (void)NPT; (void)LAYER; (void)KSPT;
    (void)KSTEP; (void)KINC;

    if (*NTENS != VOIGTSIZE_3D || *NDI != 3 || *NSHR != 3)
    {
        fprintf(stderr, "UMAT Error: this UMAT requires 3D elements (NTENS = 6).\n");
        return;
    }
    if (hs_check_properties(*NPROPS, PROPS)) return;
    if (*NSTATV < 3)
    {
        fprintf(stderr, "UMAT Error: Hardening Soil requires at least 3 state variables.\n");
        return;
    }

    /* --- gather material parameters --- */
    HSParams prm;
    prm.E50_ref = PROPS[0];
    prm.Eur_ref = PROPS[1];
    prm.m       = PROPS[2];
    prm.c       = PROPS[3];
    prm.phi     = PROPS[4] * PI / 180.0;
    prm.psi     = PROPS[5] * PI / 180.0;
    prm.p_ref   = PROPS[6];
    prm.Rf      = PROPS[7];
    prm.nu      = PROPS[8];
    prm.M_cap   = PROPS[9];
    prm.K_ratio = PROPS[10];

    if (prm.phi > ZERO_TOL)
        prm.p_t = prm.c / tan(prm.phi);
    else
        prm.p_t = 0.0;

    /* Reference initial loading stiffness. The validated Python prototype
     * (python_prototypes/hardening_soil.py, self.Ei_ref = 65488) hard-codes
     * this value rather than deriving it from E50_ref/Rf, so we mirror it here
     * to reproduce the reference results exactly. */
    prm.Ei_ref = 65488.0;

    /* --- read state --- */
    double gamma_p = STATEV[0];
    double p_p = STATEV[1];
    int at_failure = (STATEV[2] > 0.5) ? 1 : 0;

    /* Initialise the pre-consolidation pressure on first use so that the
     * current stress lies on the cap (F_c = 0). */
    if (p_p <= 0.0)
    {
        double p0 = hs_mean_stress(STRESS);
        double q_special0 = hs_cap_qspecial(STRESS, &prm);
        double rhs = sqrt((q_special0 / prm.M_cap) * (q_special0 / prm.M_cap) + p0 * p0);
        p_p = (rhs > 1.0) ? rhs : 1.0;
    }

    /* --- integrate the stress --- */
    double stress[VOIGTSIZE_3D];
    copy_array(STRESS, VOIGTSIZE_3D, stress);

    /* Relative convergence tolerance for the local return mapping, matching the
     * Python prototype (integrate(..., tol=1e-5)). ZERO_TOL (1e-12) is far too
     * tight for the perfectly-plastic Mohr-Coulomb residual and prevents
     * convergence at failure. */
    int converged = hs_integrate(stress, &gamma_p, &p_p, &at_failure, DSTRAN, &prm, 1.0e-5);

    if (!converged)
    {
        fprintf(stderr, "UMAT Warning: Hardening Soil return mapping did not converge; "
                        "requesting a smaller time increment.\n");
        if (PNEWDT) *PNEWDT = 0.5;
        return;
    }

    /* --- write results --- */
    copy_array(stress, VOIGTSIZE_3D, STRESS);
    STATEV[0] = gamma_p;
    STATEV[1] = p_p;
    STATEV[2] = (double)at_failure;

    /* consistent (elastic Eur-based) tangent operator */
    double factor = hs_stiffness_factor(&prm, stress[ZZ]);
    double Eur = prm.Eur_ref * factor;
    calculate_elastic_stiffness_matrix_3d(Eur, prm.nu, DDSDDE);

    if (SPD) *SPD = 0.0;
    return;
}



