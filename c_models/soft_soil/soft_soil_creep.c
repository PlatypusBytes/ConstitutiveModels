/**
 * @file soft_soil_creep.c
 * @brief Abaqus-style UMAT for the soft soil creep model of Stolle, Vermeer & Bonnier (1999),
 *        "Time integration of a constitutive law for soft clays", Commun. Numer. Meth. Engng 15,
 *        603-609.
 *
 * Equation numbers refer to that paper.
 *
 * Model
 * -----
 *  - Equivalent pressure on the Modified Cam-Clay ellipse (Eq. 1): p_eq = p + q^2 / (M^2 p).
 *  - Volumetric creep (Eq. 2): d eps_c / dt = C / tau (p_eq / p_p)^(B / C), with the
 *    pre-consolidation pressure p_p = p_p0 exp(eps_c / B) (Eq. 3) and p_p0 = OCR0 p_eq0 (Eq. 4).
 *    Integrated over a time increment at constant p_eq (Eq. 5):
 *        deps_c = C ln(1 + dt / tau (p_eq / p_p)^(B / C))
 *  - Elasticity with a pressure dependent bulk modulus K = p / A and a constant Poisson's ratio
 *    (Eq. 6).
 *  - Creep flow normal to the ellipse, scaled by its volumetric part (Eq. 7):
 *        deps^c = deps_c / (d p_eq/d p) d p_eq / d sigma
 *  - Mohr-Coulomb failure with a zero dilatancy angle (Eq. 15), so no negative excess pore
 *    pressures develop under undrained conditions.
 *  The parameters A, B and C correspond to kappa*, lambda* - kappa* and mu* of the PLAXIS soft
 *  soil creep model.
 *
 * Integration (the "modified" procedure of the paper)
 * ---------------------------------------------------
 *  - Elastic trial in p-q space; the stress at the end of the increment follows from (Eqs. 11-12)
 *        p = p_e - K deps_c,   q = q_e - 3 G deps_c 2 p q / (M^2 p^2 - q^2)
 *    with deps_c = deps_c(p, q) of Eq. 5 at the pre-consolidation pressure of the start of the
 *    increment. This is solved as a scalar equation in deps_c with Newton's method (Eq. 13),
 *    starting from deps_c at the stress at the start of the increment. deps_c(p, q) decreases with
 *    deps_c, so the root is unique and bracketed by [0, deps_c(p_e, q_e)], which safeguards the
 *    iteration. q is found from Eq. 12 at given p and deps_c with 0 <= q <= min(q_e, M p).
 *  - The stress follows from sigma = p I + q / q_e s_e (Eq. 14) and p_p is updated.
 *  - When the stress violates the Mohr-Coulomb criterion it is corrected with a zero dilatancy
 *    angle (Eq. 15), after the creep update as in the paper.
 *  - DDSDDE is the tangent of Eq. 9,
 *        D^ec = D - D a b^T D / (d p_eq/d p + b^T D a),   a = d p_eq / d sigma, b = d deps_c / d sigma
 *    evaluated at the stress after the creep update, followed by the elasto-plastic Mohr-Coulomb
 *    tangent when the stress is on the failure surface.
 *
 * Additions (not in the paper)
 * ----------------------------
 *  - The elastic volumetric law d p = (p / A) d eps_v^e is integrated exactly,
 *    p = p0 exp((deps_v - deps_c) / A); Eq. 11 is its linearisation. The mean stress thereby stays
 *    positive for any strain increment. The shear modulus is evaluated at the start of the
 *    increment, G = 3 K0 (1 - 2 nu) / (2 (1 + nu)) with K0 = p0 / A, so the elastic trial
 *    deviatoric stress s_e = s0 + 2 G de is linear in the strain increment.
 *  - The mean stress at the start of an increment is limited to p_min from below (a unit stress
 *    by default, as in PLAXIS), so the stiffness stays positive at zero stress.
 *  - The Mohr-Coulomb correction includes the triaxial compression / extension edges of the
 *    Mohr-Coulomb pyramid (return_mappings/mohr_coulomb_return_mapping.h).
 *  - The multiaxial generalisation uses p = tr(sigma) / 3 and q = sqrt(3 J2).
 *
 * Building blocks
 * ---------------
 *  - elastic_laws/logarithmic_elasticity.h            K = p / A (Eqs. 6 and 11)
 *  - elastic_laws/hookes_law.h                         elastic stiffness from K and G
 *  - yield_surfaces/modified_cam_clay_surface.h        equivalent pressure (Eq. 1)
 *  - flow_rules/modified_cam_clay_flow.h               deviatoric creep strain (Eqs. 7 and 12)
 *  - creep_laws/isotache_creep.h                       creep strain increment (Eqs. 2 and 5)
 *  - hardening_rules/exponential_volumetric_hardening.h pre-consolidation pressure (Eq. 3)
 *  - return_mappings/mohr_coulomb_return_mapping.h     failure (Eq. 15)
 *
 * Conventions
 * -----------
 *  - Voigt ordering: [xx, yy, zz, xy, yz, xz] (see globals.h), engineering shear strains.
 *  - The UMAT interface (STRESS, STRAN, DSTRAN) is tension positive, like Abaqus and the other
 *    models in this library. Internally the model is formulated compression positive, as the
 *    paper. Stresses are effective stresses.
 *  - The time increment is DTIME, in the unit of tau. For DTIME = 0 the response is elastic
 *    (with Mohr-Coulomb failure).
 *
 * Material properties (PROPS)
 * ---------------------------
 *   [0] A      - modified swelling index kappa* (K = p / A)              [-]
 *   [1] B      - lambda* - kappa*, hardening index of p_p (Eq. 3)        [-]
 *   [2] C      - modified creep index mu* (Eq. 2)                        [-]
 *   [3] nu     - Poisson's ratio                                         [-]
 *   [4] tau    - reference time (Eq. 2)                                  [time]
 *   [5] M      - slope of the critical state line, M = 6 sin(phi_cs) / (3 - sin(phi_cs)) (Eq. 1) [-]
 *   [6] phi    - Mohr-Coulomb friction angle, 0 disables the failure surface [degrees]
 *   [7] OCR0   - initial isotropic over-consolidation ratio p_p0 / p_eq0 (Eq. 4) [-]
 *   [8] p_min  - (optional) lower limit of the mean stress for the stiffness, default 1 [stress]
 *
 * State variables (STATEV)
 * ------------------------
 *   [0] eps_c      - accumulated volumetric creep strain (Eq. 3)
 *   [1] p_p        - isotropic pre-consolidation pressure. When 0 on the first call, it is
 *                    initialised as OCR0 p_eq of the initial stress (Eq. 4).
 *   [2] p_eq       - (optional output) equivalent pressure (Eq. 1)
 *   [3] OCR        - (optional output) over-consolidation ratio p_p / p_eq
 *   [4] at_failure - (optional output) 1 when the stress is on the Mohr-Coulomb surface
 */

#include <math.h>
#include <stdio.h>

#include "../creep_laws/isotache_creep.h"
#include "../elastic_laws/hookes_law.h"
#include "../elastic_laws/logarithmic_elasticity.h"
#include "../flow_rules/modified_cam_clay_flow.h"
#include "../globals.h"
#include "../hardening_rules/exponential_volumetric_hardening.h"
#include "../return_mappings/mohr_coulomb_return_mapping.h"
#include "../strain_utils.h"
#include "../stress_utils.h"
#include "../utils.h"
#include "../yield_surfaces/modified_cam_clay_surface.h"

/* Calling convention / export macros ----------------------------------------- */
#if defined(_WIN32) || defined(_WIN64)
#define UMAT_EXPORT __declspec(dllexport)
#define UMAT_CALLCONV __stdcall
#else
#define UMAT_EXPORT
#define UMAT_CALLCONV
#endif

#define SSC_MAX_LOCAL_ITER 200

/* Tolerance on the creep strain increment, relative to the creep index C. */
#define SSC_CREEP_TOL 1.0e-12

/* Default lower limit of the mean stress (a unit stress, as in PLAXIS). */
#define SSC_DEFAULT_MIN_PRESSURE 1.0

/* Fraction of the elastic stiffness kept in the tangent at the apex of the Mohr-Coulomb pyramid,
 * where the exact tangent vanishes. With the zero dilatancy angle of the paper the mean stress is
 * not changed by the Mohr-Coulomb correction, so the apex (p = 0) is not reached. */
#define SSC_APEX_STIFFNESS_FRACTION 1.0e-2

/* Material parameters gathered in a struct for convenience. */
typedef struct
{
    double A;
    double B;
    double C;
    double nu;
    double tau;
    double M;
    double phi; /* rad */
    double OCR0;
    double p_min;

    /* derived quantities */
    double sin_phi;
    double cos_phi;
    int use_failure;
} SSCParams;

/* Stress state of the creep update (Eqs. 11-12) for a given creep strain increment deps_c. */
typedef struct
{
    double p;
    double q;
    double creep;      /* deps_c(p, q) of Eq. 5 */
    double dcreep_dx;  /* total derivative of deps_c(p, q) w.r.t. the iterate */
} SSCCreepIterate;

/* Quantities that are fixed during the creep update of one increment. */
typedef struct
{
    const SSCParams* prm;
    double p0;     /* mean stress at the start of the increment */
    double q_e;    /* elastic trial deviatoric stress */
    double deps_v; /* volumetric strain increment */
    double G;      /* shear modulus */
    double p_p;    /* pre-consolidation pressure at the start of the increment */
    double dt;
} SSCStep;

/* ------------------------------------------------------------------ */
/* Creep update in p-q space                                           */
/* ------------------------------------------------------------------ */

/*
 * Stress for a creep strain increment x (Eqs. 11-12) and the creep strain increment of Eq. 5 at
 * that stress, with the derivative needed for the Newton iteration of Eq. 13.
 */
static void ssc_evaluate_creep(const SSCStep* st, double x, SSCCreepIterate* it)
{
    const SSCParams* prm = st->prm;

    /* Eq. 11, with the logarithmic elastic law integrated exactly */
    it->p = logarithmic_mean_stress(st->p0, st->deps_v - x, prm->A, prm->p_min);
    double dp_dx = -it->p / prm->A;

    /* Eq. 12 */
    double dq_dp, dq_dx_partial;
    it->q = modified_cam_clay_deviatoric_return(st->q_e, it->p, st->G, x, prm->M, &dq_dp,
                                                &dq_dx_partial);
    double dq_dx = dq_dp * dp_dx + dq_dx_partial;

    /* Eqs. 1 and 5 */
    double dpeq_dp, dpeq_dq, dcreep_dpeq;
    double p_eq = modified_cam_clay_equivalent_pressure(it->p, it->q, prm->M, &dpeq_dp, &dpeq_dq);
    it->creep = isotache_creep_increment(p_eq, st->p_p, prm->B, prm->C, prm->tau, st->dt,
                                         &dcreep_dpeq);
    it->dcreep_dx = dcreep_dpeq * (dpeq_dp * dp_dx + dpeq_dq * dq_dx);
}

/*
 * Solves x = deps_c(p(x), q(x)) for the creep strain increment x (Eq. 13). The residual
 * F(x) = deps_c(p(x), q(x)) - x decreases strictly, with F(0) >= 0 and F(F(0)) <= 0, so the root
 * is unique and bracketed; Newton steps that leave the bracket are replaced by bisection.
 */
static int ssc_solve_creep(const SSCStep* st, double x_guess, double* x_out, SSCCreepIterate* it)
{
    double tol = SSC_CREEP_TOL * st->prm->C;

    ssc_evaluate_creep(st, 0.0, it);
    double lo = 0.0, hi = it->creep;
    if (hi <= tol)
    {
        *x_out = 0.0;
        return 1;
    }

    double x = (x_guess > lo && x_guess < hi) ? x_guess : hi;
    for (int iter = 0; iter < SSC_MAX_LOCAL_ITER; ++iter)
    {
        ssc_evaluate_creep(st, x, it);
        double F = it->creep - x;
        if (fabs(F) <= tol)
        {
            *x_out = x;
            return 1;
        }
        if (F > 0.0)
            lo = x;
        else
            hi = x;

        double x_new = x - F / (it->dcreep_dx - 1.0);
        if (!(x_new > lo && x_new < hi)) x_new = 0.5 * (lo + hi);
        if (hi - lo <= tol)
        {
            ssc_evaluate_creep(st, x_new, it);
            *x_out = x_new;
            return 1;
        }
        x = x_new;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Tangent                                                             */
/* ------------------------------------------------------------------ */

/*
 * Tangent of Eq. 9 at the stress after the creep update, with the pre-consolidation pressure of
 * the start of the increment (as in Eq. 5):
 *   D^ec = D - c (D a)(D a)^T / (d p_eq/d p + c a^T D a),   a = d p_eq / d sigma,
 * with c = d deps_c / d p_eq, so that b = c a. When the stress is on the failure surface the
 * Mohr-Coulomb return is applied to it.
 */
static void ssc_tangent(const SSCStep* st, const double sigma_creep[VOIGTSIZE_3D],
                        const double De[VOIGTSIZE_3D * VOIGTSIZE_3D], const MohrCoulombReturn* mc,
                        double Q[3][3], double ddsdde[VOIGTSIZE_3D * VOIGTSIZE_3D])
{
    const SSCParams* prm = st->prm;
    copy_array(De, VOIGTSIZE_3D * VOIGTSIZE_3D, ddsdde);

    double a[VOIGTSIZE_3D], dpeq_dp, c;
    double p_eq = modified_cam_clay_equivalent_pressure_3d(sigma_creep, prm->M, a, &dpeq_dp);
    isotache_creep_increment(p_eq, st->p_p, prm->B, prm->C, prm->tau, st->dt, &c);
    if (c > 0.0)
    {
        double Da[VOIGTSIZE_3D];
        matrix_vector_multiply(De, a, VOIGTSIZE_3D, Da);
        double denominator = dpeq_dp + c * vector_dot_product(a, Da, VOIGTSIZE_3D);
        if (denominator > SMALL_VALUE)
            for (int i = 0; i < VOIGTSIZE_3D; ++i)
                for (int j = 0; j < VOIGTSIZE_3D; ++j)
                    ddsdde[i * VOIGTSIZE_3D + j] -= c * Da[i] * Da[j] / denominator;
    }

    if (mc->type == MC_RETURN_APEX)
    {
        for (int i = 0; i < VOIGTSIZE_3D * VOIGTSIZE_3D; ++i)
            ddsdde[i] = SSC_APEX_STIFFNESS_FRACTION * De[i];
    }
    else if (mc->type != MC_RETURN_ELASTIC)
    {
        mohr_coulomb_project_tangent(mc, Q, De, prm->sin_phi, 0.0, ddsdde);
    }
}

/* ------------------------------------------------------------------ */
/* Stress update                                                       */
/* ------------------------------------------------------------------ */

/*
 * Integrates the strain increment over the time increment dt. On success the stress, eps_c and
 * p_p are updated in place and, if ddsdde is given, the tangent is stored in it.
 */
static int ssc_integrate(const SSCParams* prm, double stress[VOIGTSIZE_3D], double* eps_c,
                         double* p_p, int* at_failure, const double dstrain[VOIGTSIZE_3D],
                         double dt, double* ddsdde)
{
    SSCStep st;
    double s0[VOIGTSIZE_3D], s_e[VOIGTSIZE_3D], ds[VOIGTSIZE_3D];
    double D_dev[VOIGTSIZE_3D * VOIGTSIZE_3D], De[VOIGTSIZE_3D * VOIGTSIZE_3D];

    /* elasticity at the start of the increment */
    st.prm = prm;
    st.p0 = calculate_mean_stress(stress);
    st.G = calculate_shear_modulus_from_bulk_modulus(
        logarithmic_bulk_modulus(st.p0, prm->A, prm->p_min), prm->nu);
    st.deps_v = calculate_volumetric_strain(dstrain);
    st.p_p = *p_p;
    st.dt = dt;

    /* elastic trial deviatoric stress s_e = s0 + 2 G de */
    calculate_deviatoric_stress(stress, st.p0, s0);
    calculate_elastic_stiffness_matrix_3d_bulk_shear(0.0, st.G, D_dev);
    matrix_vector_multiply(D_dev, dstrain, VOIGTSIZE_3D, ds);
    add_vectors(s0, ds, VOIGTSIZE_3D, s_e);
    st.q_e = calculate_von_mises_stress(s_e);

    /* creep update, starting from the creep strain at the stress at the start of the increment */
    double p_start = fmax(st.p0, prm->p_min);
    double p_eq0 = modified_cam_clay_equivalent_pressure(p_start, calculate_von_mises_stress(s0),
                                                         prm->M, NULL, NULL);
    double x_guess = isotache_creep_increment(p_eq0, st.p_p, prm->B, prm->C, prm->tau, dt, NULL);

    SSCCreepIterate it;
    double x;
    if (!ssc_solve_creep(&st, x_guess, &x, &it)) return 0;

    /* Eq. 14 */
    double sigma_creep[VOIGTSIZE_3D];
    double q_ratio = (st.q_e > 0.0) ? it.q / st.q_e : 0.0;
    for (int i = 0; i < VOIGTSIZE_3D; ++i) sigma_creep[i] = q_ratio * s_e[i];
    for (int i = 0; i < 3; ++i) sigma_creep[i] += it.p;

    /* elastic stiffness at the end of the increment: tangent of the logarithmic law */
    double K = logarithmic_bulk_modulus(it.p, prm->A, prm->p_min);
    calculate_elastic_stiffness_matrix_3d_bulk_shear(K, st.G, De);

    /* Mohr-Coulomb failure with zero dilatancy, Eq. 15 */
    MohrCoulombReturn mc;
    mc.type = MC_RETURN_ELASTIC;
    mc.n_active = 0;
    double Q[3][3];
    double sigma[VOIGTSIZE_3D];
    copy_array(sigma_creep, VOIGTSIZE_3D, sigma);
    if (prm->use_failure)
    {
        double s_creep[3], s[3], D_principal[9];
        calculate_principal_system(sigma_creep, s_creep, Q);
        calculate_elastic_stiffness_matrix_principal_bulk_shear(K, st.G, D_principal);
        if (!mohr_coulomb_return_mapping(s_creep, D_principal, prm->sin_phi, prm->cos_phi, 0.0,
                                         0.0, s, &mc))
            return 0;
        if (mc.type != MC_RETURN_ELASTIC) calculate_stress_from_principal_system(s, Q, sigma);
    }

    /* update the stress and the state; hardening, Eq. 3 */
    copy_array(sigma, VOIGTSIZE_3D, stress);
    *eps_c += x;
    *p_p = exponential_volumetric_hardening(*p_p, x, prm->B);
    *at_failure = (mc.type != MC_RETURN_ELASTIC) ? 1 : 0;

    if (ddsdde) ssc_tangent(&st, sigma_creep, De, &mc, Q, ddsdde);
    return 1;
}

/* ------------------------------------------------------------------ */
/* Parameters and initial state                                        */
/* ------------------------------------------------------------------ */

static int ssc_check_properties(int nprops, const double* props)
{
    if (nprops < 8)
    {
        fprintf(stderr, "UMAT Error: Soft Soil Creep requires 8 properties.\n");
        return 1;
    }
    int n_errors = 0;
    if (props[0] <= 0.0) { fprintf(stderr, "UMAT Error: A (kappa*) must be positive.\n"); n_errors++; }
    if (props[1] <= 0.0) { fprintf(stderr, "UMAT Error: B (lambda* - kappa*) must be positive.\n"); n_errors++; }
    if (props[2] <= 0.0) { fprintf(stderr, "UMAT Error: C (mu*) must be positive.\n"); n_errors++; }
    if (props[3] < 0.0 || props[3] >= 0.5) { fprintf(stderr, "UMAT Error: nu in [0, 0.5).\n"); n_errors++; }
    if (props[4] <= 0.0) { fprintf(stderr, "UMAT Error: tau must be positive.\n"); n_errors++; }
    if (props[5] <= 0.0) { fprintf(stderr, "UMAT Error: M must be positive.\n"); n_errors++; }
    if (props[6] < 0.0 || props[6] >= 90.0) { fprintf(stderr, "UMAT Error: phi in [0, 90) (0 disables failure).\n"); n_errors++; }
    if (props[7] <= 0.0) { fprintf(stderr, "UMAT Error: OCR0 must be positive.\n"); n_errors++; }
    if (nprops > 8 && props[8] <= 0.0) { fprintf(stderr, "UMAT Error: p_min must be positive.\n"); n_errors++; }
    return (n_errors > 0) ? 1 : 0;
}

static void ssc_set_parameters(int nprops, const double* props, SSCParams* prm)
{
    prm->A = props[0];
    prm->B = props[1];
    prm->C = props[2];
    prm->nu = props[3];
    prm->tau = props[4];
    prm->M = props[5];
    prm->phi = props[6] * PI / 180.0;
    prm->OCR0 = props[7];
    prm->p_min = (nprops > 8) ? props[8] : SSC_DEFAULT_MIN_PRESSURE;

    prm->sin_phi = sin(prm->phi);
    prm->cos_phi = cos(prm->phi);
    prm->use_failure = (prm->phi > 0.0) ? 1 : 0;
}

/* p_eq (Eq. 1) of a stress, with the mean stress limited to p_min. */
static double ssc_equivalent_pressure(const SSCParams* prm, const double stress[VOIGTSIZE_3D])
{
    double p = fmax(calculate_mean_stress(stress), prm->p_min);
    return modified_cam_clay_equivalent_pressure(p, calculate_von_mises_stress(stress), prm->M,
                                                 NULL, NULL);
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
    (void)SSE; (void)SPD; (void)SCD; (void)RPL; (void)DDSDDT; (void)DRPLDE; (void)DRPLDT;
    (void)STRAN; (void)TIME; (void)TEMP; (void)DTEMP; (void)PREDEF; (void)DPRED; (void)CMNAME;
    (void)COORDS; (void)DROT; (void)CELENT; (void)DFGRD0; (void)DFGRD1; (void)NOEL; (void)NPT;
    (void)LAYER; (void)KSPT; (void)KSTEP; (void)KINC;

    if (*NTENS != VOIGTSIZE_3D || *NDI != 3 || *NSHR != 3)
    {
        fprintf(stderr, "UMAT Error: this UMAT requires 3D elements (NTENS = 6).\n");
        return;
    }
    if (ssc_check_properties(*NPROPS, PROPS)) return;
    if (*NSTATV < 2)
    {
        fprintf(stderr, "UMAT Error: Soft Soil Creep requires at least 2 state variables.\n");
        return;
    }

    /* --- gather material parameters --- */
    SSCParams prm;
    ssc_set_parameters(*NPROPS, PROPS, &prm);

    /* --- read state --- */
    double eps_c = STATEV[0];
    double p_p = STATEV[1];
    int at_failure = 0;
    double dt = (DTIME && *DTIME > 0.0) ? *DTIME : 0.0;

    /* --- convert to the sign convention of the model ---
     * The UMAT interface uses the Abaqus convention of the other models in this library (tension
     * positive); the model is formulated compression positive, as the paper. The tangent
     * d sigma / d eps is the same in both conventions. */
    double stress[VOIGTSIZE_3D], dstrain[VOIGTSIZE_3D];
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        stress[i] = -STRESS[i];
        dstrain[i] = -DSTRAN[i];
    }

    /* first call: initial pre-consolidation pressure from the initial stress, Eq. 4 */
    if (p_p <= 0.0)
    {
        eps_c = 0.0;
        p_p = prm.OCR0 * ssc_equivalent_pressure(&prm, stress);
    }

    /* --- integrate the stress --- */
    if (!ssc_integrate(&prm, stress, &eps_c, &p_p, &at_failure, dstrain, dt, DDSDDE))
    {
        fprintf(stderr, "UMAT Warning: Soft Soil Creep stress update did not converge; "
                        "requesting a smaller time increment.\n");
        /* keep the stress and state at the start of the increment, with the elastic tangent */
        double K = logarithmic_bulk_modulus(calculate_mean_stress(stress), prm.A, prm.p_min);
        calculate_elastic_stiffness_matrix_3d_bulk_shear(
            K, calculate_shear_modulus_from_bulk_modulus(K, prm.nu), DDSDDE);
        if (PNEWDT) *PNEWDT = 0.5;
        return;
    }

    /* --- write results --- */
    for (int i = 0; i < VOIGTSIZE_3D; ++i) STRESS[i] = -stress[i];
    STATEV[0] = eps_c;
    STATEV[1] = p_p;
    double p_eq = ssc_equivalent_pressure(&prm, stress);
    if (*NSTATV > 2) STATEV[2] = p_eq;
    if (*NSTATV > 3) STATEV[3] = p_p / p_eq;
    if (*NSTATV > 4) STATEV[4] = (double)at_failure;
}
