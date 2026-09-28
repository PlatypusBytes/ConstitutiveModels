/**
 * @file hardening_soil.c
 * @brief Abaqus-style UMAT for the Hardening Soil model of Schanz, Vermeer & Bonnier (1999),
 *        "The hardening soil model: Formulation and verification".
 *
 * Equation numbers refer to that paper.
 *
 * Model
 * -----
 *  - Stress dependent stiffness (Eqs. 3-5): E50, Eur ~ ((sigma_3 + a) / (p_ref + a))^m with
 *    a = c cot(phi). Elasticity is isotropic with Eur and nu_ur.
 *  - Shear hardening (Eqs. 7-9):
 *        f_ij = 2/Ei * q_ij / (1 - q_ij/qa) - 2 q_ij / Eur - gamma_p,   q_ij = sigma_i - sigma_j
 *    with qa = qf / Rf (Eq. 2) and Ei = 2 E50 / (2 - Rf), so that E50 is the secant stiffness at
 *    q = qf / 2 as defined in Sec. 2.1 (for Rf -> 1 this is the qa/E50 form of Eqs. 7-8).
 *  - Non-associated flow (Eqs. 11-15): g_ij = (sigma_i - sigma_j)/2 - (sigma_i + sigma_j)/2 sin(psi_m)
 *    with the mobilised dilatancy angle psi_m of Rowe's stress-dilatancy theory. The plastic
 *    shear strain is gamma_p = eps1_p - eps2_p - eps3_p (Eq. 9), i.e. d(gamma_p) = dLambda_ij.
 *  - Mohr-Coulomb failure q <= qf (Eq. 2, end of Sec. 3), with the same (mobilised) flow as the
 *    hardening surface, which equals the dilation angle psi at failure (Eq. 13).
 *  - Elliptic cap with associated flow and hardening of the pre-consolidation stress p_c
 *    (Eqs. 27-35): dp_c = H deps_v^pc with H = Ks / (Ks/Kc - 1) ((p_c + a) / (p_ref + a))^m, so
 *    that the stiffness increases with the stress level like Eoed (Eq. 37).
 *  - Optional dilatancy cut-off at the critical void ratio (Eqs. 38-39).
 *
 * Integration (Sec. 3)
 * --------------------
 *  - The stiffness and the mobilised dilatancy are evaluated at the stress at the start of a
 *    (sub-)step (Euler explicit); qa is evaluated at the end of the step (Eq. 23).
 *  - The return mapping is carried out in principal stress space on the elastic trial stress
 *    (Eqs. 20-21) and rotated back with the principal directions of the trial stress. At the
 *    triaxial compression / extension corners two shear surfaces are active (Koiter, Eq. 26),
 *    and q~ of the cap has a corner as well (Fig. 3).
 *  - Automatic sub-stepping on the size of the elastic stress increment.
 *  - DDSDDE is the elasto-plastic tangent of the active yield surfaces at the end of the
 *    increment (the paper uses the elastic stiffness in the global iterations, Eq. 16, which
 *    converges slowly in a full Newton scheme).
 *
 * Additions for boundary value problems (not in the paper)
 * --------------------------------------------------------
 *  - The stress level in the stiffness laws is limited to
 *    sigma_3 + a >= HS_MIN_STRESS_RATIO (p_ref + a), so the stiffness stays positive at zero
 *    confinement and in tension.
 *  - Trial stresses in tension beyond the apex of the Mohr-Coulomb pyramid are returned to the
 *    apex.
 *
 * Conventions
 * -----------
 *  - Voigt ordering: [xx, yy, zz, xy, yz, xz] (see globals.h), engineering shear strains.
 *  - The UMAT interface (STRESS, STRAN, DSTRAN) is tension positive, like Abaqus and the other
 *    models in this library. Internally the model is formulated compression positive, as the
 *    paper. Principal stresses are ordered sigma_1 >= sigma_2 >= sigma_3 (compression positive),
 *    so sigma_3 is the minor principal stress. For this convention Eq. 12 reads
 *    sin(phi_m) = (sigma_1 - sigma_3) / (sigma_1 + sigma_3 + 2a), and dilatancy is a negative
 *    plastic volumetric strain, eps_v^p = -sin(psi_m) gamma_p (Eq. 15).
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
 *   [9]  M_cap    - cap aspect ratio M (Eq. 27), 0 disables the cap [-]
 *   [10] K_ratio  - Ks/Kc, swelling over elasto-plastic bulk modulus in isotropic
 *                   compression (Eq. 32)                          [-]
 *   [11] e0       - (optional) initial void ratio (Eq. 39)        [-]
 *   [12] e_cv     - (optional) critical void ratio of the dilatancy cut-off (Eq. 38),
 *                   0 disables the cut-off                         [-]
 *
 * State variables (STATEV)
 * ------------------------
 *   [0]  gamma_p    - plastic shear strain (Eq. 9)
 *   [1]  p_c        - pre-consolidation stress (cap size). When 0 on the first call, gamma_p and
 *                     p_c are initialised such that the initial stress lies on the shear
 *                     hardening surface and on the cap (normally consolidated state).
 *   [2]  at_failure - 1 when the stress is on the Mohr-Coulomb failure surface, 0 otherwise
 */

#include <math.h>
#include <stdio.h>

#include "../elastic_laws/hookes_law.h"
#include "../globals.h"
#include "../stress_utils.h"
#include "../utils.h"

/* Calling convention / export macros ----------------------------------------- */
#if defined(_WIN32) || defined(_WIN64)
#define UMAT_EXPORT __declspec(dllexport)
#define UMAT_CALLCONV __stdcall
#else
#define UMAT_EXPORT
#define UMAT_CALLCONV
#endif

#define HS_MAX_LOCAL_ITER 50
#define HS_MAX_LINE_SEARCH 30
#define HS_MAX_ACTIVE_SET_ITER 10
#define HS_MAX_SUBSTEPS 1000
#define HS_MAX_SUBSTEP_TRIALS 5

/* Maximum elastic stress increment per sub-step, relative to sigma_3 + a. The stiffness and the
 * dilatancy are frozen during a sub-step (Euler explicit, Sec. 3), so this bounds the explicit
 * integration error. */
#define HS_SUBSTEP_STRESS_FRACTION 0.01

/* Relative tolerance on the yield functions, and the tolerance at which a stagnating local
 * iteration (round-off level) is still accepted. */
#define HS_YIELD_TOL 1.0e-12
#define HS_YIELD_TOL_ACCEPT 1.0e-8

/* Lower limit of (sigma_3 + a) / (p_ref + a) in the stress dependency of the stiffness (Eqs. 3-4,
 * 35). Not part of the paper: without it the stiffness vanishes at zero confinement (free surfaces,
 * c = 0) and in tension, which makes the global stiffness matrix singular. */
#define HS_MIN_STRESS_RATIO 1.0e-2

/* Fraction of the elastic stiffness kept in the tangent for deformation modes in which the stress
 * does not change: the corner mode at a triaxial corner (see hs_elastoplastic_tangent) and all
 * modes at the tension apex. The exact tangent has zero stiffness there, which can make the global
 * stiffness matrix singular. */
#ifndef HS_CORNER_STIFFNESS_FRACTION
#define HS_CORNER_STIFFNESS_FRACTION 1.0e-2
#endif

/* At most the shear surface, the cap and the corner equality are active at the same time. */
#define HS_MAX_ACTIVE 3

#define HS_DEBUG

/* Material parameters gathered in a struct for convenience. */
typedef struct
{
    double E50_ref;
    double Eur_ref;
    double m;
    double c;
    double phi; /* rad */
    double psi; /* rad */
    double p_ref;
    double Rf;
    double nu;
    double M_cap;
    double K_ratio;
    double e0;
    double e_cv;

    /* derived quantities */
    double sin_phi;
    double cos_phi;
    double sin_psi;
    double sin_phi_cv; /* critical state friction angle, Eq. 13 */
    double a;          /* c cot(phi), Eq. 24 */
    double k_f;        /* qf = k_f * (sigma_3 + a), Eqs. 2 and 23 */
    double alpha;      /* cap shape factor, Eq. 29 */
    double Ei_ref;     /* reference initial loading stiffness 2 E50_ref / (2 - Rf) */
    double H_ref;      /* reference cap hardening modulus, Eq. 32 */
    int use_cap;
    int use_cutoff;
} HSParams;

/* Quantities that are frozen during one (sub-)step (Euler explicit, Sec. 3). */
typedef struct
{
    const HSParams* prm;
    double Eur;
    double Ei;
    double D[9];      /* elasticity in principal stress space (row-major 3x3) */
    double sin_psi_m; /* flow of the shear surfaces, Eqs. 11 and 38; equals psi at failure (Eq. 13) */
    double H;         /* cap hardening modulus, Eqs. 32 and 35 */
    double gamma_p0;
    double p_c0;
} HSStep;

typedef enum
{
    HS_CONE,  /* shear hardening surface f_13, Eqs. 7-8 */
    HS_MC,    /* Mohr-Coulomb failure surface */
    HS_CAP,   /* cap, Eq. 27 */
    HS_CORNER /* principal stress equality at a triaxial corner, see hs_return_mapping */
} HSSurfaceType;

/* Triaxial corners of the yield surfaces in principal stress space. */
typedef enum
{
    HS_CORNER_NONE,
    HS_CORNER_COMPRESSION, /* sigma_2 = sigma_3 */
    HS_CORNER_EXTENSION    /* sigma_1 = sigma_2 */
} HSCorner;

typedef struct
{
    HSSurfaceType type;
    HSCorner corner; /* shear and cap surfaces: corner at which the flow is averaged */
    double w[3];     /* cap: q~ = w . sigma (Eq. 28, or its average at a corner) */
} HSSurface;

/* Stress and hardening state of the return mapping for given plastic multipliers. */
typedef struct
{
    double s[3]; /* principal stresses */
    double gamma_p;
    double p_c;
    double A_inv[9];            /* inverse of the (cap) system matrix, see hs_evaluate_iterate */
    double n[HS_MAX_ACTIVE][3]; /* flow directions dg/dsigma */
    double f[HS_MAX_ACTIVE];    /* yield function values */
} HSIterate;

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

/* Minor principal stress sigma_3 (smallest, compression positive). */
static double hs_minor_principal_stress(const double s[VOIGTSIZE_3D])
{
    double ps[3];
    double Q[3][3];
    calculate_principal_system(s, ps, Q);
    return ps[2];
}

/* ------------------------------------------------------------------ */
/* Small dense linear algebra (principal stress space)                 */
/* ------------------------------------------------------------------ */

static void hs_matvec3(const double A[9], const double x[3], double y[3])
{
    for (int r = 0; r < 3; ++r) y[r] = A[3 * r] * x[0] + A[3 * r + 1] * x[1] + A[3 * r + 2] * x[2];
}

static int hs_invert3(const double A[9], double A_inv[9])
{
    double c00 = A[4] * A[8] - A[5] * A[7];
    double c01 = A[5] * A[6] - A[3] * A[8];
    double c02 = A[3] * A[7] - A[4] * A[6];
    double det = A[0] * c00 + A[1] * c01 + A[2] * c02;
    if (fabs(det) < SMALL_VALUE) return 0;

    double inv = 1.0 / det;
    A_inv[0] = c00 * inv;
    A_inv[1] = (A[2] * A[7] - A[1] * A[8]) * inv;
    A_inv[2] = (A[1] * A[5] - A[2] * A[4]) * inv;
    A_inv[3] = c01 * inv;
    A_inv[4] = (A[0] * A[8] - A[2] * A[6]) * inv;
    A_inv[5] = (A[2] * A[3] - A[0] * A[5]) * inv;
    A_inv[6] = c02 * inv;
    A_inv[7] = (A[1] * A[6] - A[0] * A[7]) * inv;
    A_inv[8] = (A[0] * A[4] - A[1] * A[3]) * inv;
    return 1;
}

/* Solves J x = b for n <= HS_MAX_ACTIVE (Gaussian elimination with partial pivoting). J and b
 * are overwritten. */
static int hs_solve_linear(int n, double J[HS_MAX_ACTIVE * HS_MAX_ACTIVE], double b[HS_MAX_ACTIVE],
                           double x[HS_MAX_ACTIVE])
{
    for (int col = 0; col < n; ++col)
    {
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (fabs(J[r * HS_MAX_ACTIVE + col]) > fabs(J[piv * HS_MAX_ACTIVE + col])) piv = r;
        if (fabs(J[piv * HS_MAX_ACTIVE + col]) < SMALL_VALUE) return 0;

        if (piv != col)
        {
            for (int k = 0; k < n; ++k)
            {
                double tmp = J[col * HS_MAX_ACTIVE + k];
                J[col * HS_MAX_ACTIVE + k] = J[piv * HS_MAX_ACTIVE + k];
                J[piv * HS_MAX_ACTIVE + k] = tmp;
            }
            double tmp = b[col];
            b[col] = b[piv];
            b[piv] = tmp;
        }

        for (int r = col + 1; r < n; ++r)
        {
            double factor = J[r * HS_MAX_ACTIVE + col] / J[col * HS_MAX_ACTIVE + col];
            for (int k = col; k < n; ++k) J[r * HS_MAX_ACTIVE + k] -= factor * J[col * HS_MAX_ACTIVE + k];
            b[r] -= factor * b[col];
        }
    }

    for (int r = n - 1; r >= 0; --r)
    {
        double sum = b[r];
        for (int k = r + 1; k < n; ++k) sum -= J[r * HS_MAX_ACTIVE + k] * x[k];
        x[r] = sum / J[r * HS_MAX_ACTIVE + r];
    }
    return 1;
}

/* Isotropic elasticity acting on the principal stresses / strains. */
static void hs_principal_elastic_matrix(double E, double nu, double D[9])
{
    double G = E / (2.0 * (1.0 + nu));
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    for (int r = 0; r < 3; ++r)
        for (int col = 0; col < 3; ++col) D[3 * r + col] = lambda + ((r == col) ? 2.0 * G : 0.0);
}

/* ------------------------------------------------------------------ */
/* Stress-dependent stiffness and dilatancy                            */
/* ------------------------------------------------------------------ */

/* Stress dependency of E50 and Eur, Eqs. 3-4: ((sigma_3 + a) / (p_ref + a))^m. */
static double hs_stiffness_factor(const HSParams* prm, double sigma_3)
{
    double ratio = (sigma_3 + prm->a) / (prm->p_ref + prm->a);
    if (ratio < HS_MIN_STRESS_RATIO) ratio = HS_MIN_STRESS_RATIO;
    return pow(ratio, prm->m);
}

/* Cap hardening modulus H = Ks Kc / (Ks - Kc) (Eq. 32), stress dependent through the
 * pre-consolidation stress as in Eq. 35. */
static double hs_cap_modulus(const HSParams* prm, double p_c)
{
    double ratio = (p_c + prm->a) / (prm->p_ref + prm->a);
    if (ratio < HS_MIN_STRESS_RATIO) ratio = HS_MIN_STRESS_RATIO;
    return prm->H_ref * pow(ratio, prm->m);
}

/* Plastic shear strain on the hyperbola for a deviator stress q < qa (Eq. 8 with f = 0). */
static double hs_hyperbolic_gamma_p(double q, double qa, double Ei, double Eur)
{
    return 2.0 / Ei * q / (1.0 - q / qa) - 2.0 * q / Eur;
}

/* Mobilised dilatancy angle from Rowe's stress-dilatancy theory (Eqs. 11-12), including the
 * dilatancy cut-off of Eq. 38. Contractant (negative) for phi_m < phi_cv. */
static double hs_sin_psi_mobilised(const HSParams* prm, const double s[3], double void_ratio)
{
    if (prm->use_cutoff && void_ratio >= prm->e_cv) return 0.0;

    double sin_phi_m = prm->sin_phi;
    double denom = s[0] + s[2] + 2.0 * prm->a;
    if (denom > ZERO_TOL) sin_phi_m = (s[0] - s[2]) / denom;
    if (sin_phi_m < 0.0) sin_phi_m = 0.0;
    if (sin_phi_m > prm->sin_phi) sin_phi_m = prm->sin_phi;

    return (sin_phi_m - prm->sin_phi_cv) / (1.0 - sin_phi_m * prm->sin_phi_cv);
}

/* ------------------------------------------------------------------ */
/* Yield surfaces and plastic potentials (principal stresses)          */
/* ------------------------------------------------------------------ */

/*
 * Shear hardening yield function f13 (Eq. 8):
 *
 *   f = 2/Ei * q / (1 - q/qa) - 2 q / Eur - gamma_p,   q = sigma_1 - sigma_3,
 *   qa = k_f (sigma_3 + a) / Rf                                       (Eqs. 2, 23)
 *
 * At the triaxial corners f12 (Eq. 7) and f23 coincide with f13. The function is evaluated
 * multiplied by (qa - q), which removes the pole at the asymptote and gives the quadratic form
 * mentioned in Sec. 3 (Eq. 25). For q < qa the sign is unchanged, and for q >= qa the scaled
 * function is strictly positive, so a stress beyond the asymptote is always detected as yielding.
 * E50 and Eur are taken at the start of the step, qa at the current stress.
 */
static double hs_cone_function(const HSStep* st, const double s[3], double gamma_p, double df_ds[3],
                               double* df_dgamma)
{
    const HSParams* prm = st->prm;
    double k_a = prm->k_f / prm->Rf;
    double q = s[0] - s[2];
    double qa = k_a * (s[2] + prm->a);
    double strain = 2.0 * q / st->Eur + gamma_p;

    if (df_ds)
    {
        double df_dq = 2.0 / st->Ei * qa - 2.0 / st->Eur * (qa - q) + strain;
        double df_dqa = 2.0 / st->Ei * q - strain;
        df_ds[0] = df_dq;
        df_ds[1] = 0.0;
        df_ds[2] = -df_dq + df_dqa * k_a;
    }
    if (df_dgamma) *df_dgamma = -(qa - q);

    return 2.0 / st->Ei * q * qa - strain * (qa - q);
}

/* Mohr-Coulomb failure surface, equivalent to q = sigma_1 - sigma_3 <= qf (Eq. 2):
 *   f = (sigma_1 - sigma_3)/2 - (sigma_1 + sigma_3)/2 sin(phi) - c cos(phi) */
static double hs_mc_function(const HSParams* prm, const double s[3], double df_ds[3])
{
    if (df_ds)
    {
        df_ds[0] = 0.5 - 0.5 * prm->sin_phi;
        df_ds[1] = 0.0;
        df_ds[2] = -0.5 - 0.5 * prm->sin_phi;
    }
    return 0.5 * (s[0] - s[2]) - 0.5 * (s[0] + s[2]) * prm->sin_phi - prm->c * prm->cos_phi;
}

/* Cap yield function (Eqs. 27-29), also the plastic potential (associated flow, Eq. 30):
 *   fc = q~^2 / M^2 + (p + a)^2 - (p_c + a)^2,   q~ = w . sigma */
static double hs_cap_function(const HSParams* prm, const double w[3], const double s[3], double p_c,
                              double df_ds[3], double* df_dpc)
{
    double M2 = prm->M_cap * prm->M_cap;
    double q_tilde = w[0] * s[0] + w[1] * s[1] + w[2] * s[2];
    double p = (s[0] + s[1] + s[2]) / 3.0;

    if (df_ds)
    {
        double dp_term = 2.0 * (p + prm->a) / 3.0;
        for (int r = 0; r < 3; ++r) df_ds[r] = 2.0 * q_tilde / M2 * w[r] + dp_term;
    }
    if (df_dpc) *df_dpc = -2.0 * (p_c + prm->a);

    return q_tilde * q_tilde / M2 + (p + prm->a) * (p + prm->a) - (p_c + prm->a) * (p_c + prm->a);
}

/*
 * q~ = w . sigma of the cap: w = [1, alpha - 1, -alpha] for ordered principal stresses (Eq. 28).
 * q~ is not smooth at the triaxial corners; there the average over both orderings of the equal
 * principal stresses is used, which gives the values of the paper at the corners:
 * q~ = sigma_1 - sigma_3 (compression) and q~ = alpha (sigma_1 - sigma_3) (extension).
 */
static void hs_cap_weights(const HSParams* prm, HSCorner corner, double w[3])
{
    double alpha = prm->alpha;
    w[0] = 1.0;
    w[1] = alpha - 1.0;
    w[2] = -alpha;
    if (corner == HS_CORNER_COMPRESSION)
    {
        w[1] = -0.5;
        w[2] = -0.5;
    }
    else if (corner == HS_CORNER_EXTENSION)
    {
        w[0] = 0.5 * alpha;
        w[1] = 0.5 * alpha;
    }
}

/* Principal stress equality at a triaxial corner, sigma_2 - sigma_3 = 0 or sigma_1 - sigma_2 = 0. */
static double hs_corner_function(HSCorner corner, const double s[3], double df_ds[3])
{
    int i = (corner == HS_CORNER_COMPRESSION) ? 1 : 0;
    if (df_ds)
    {
        df_ds[0] = df_ds[1] = df_ds[2] = 0.0;
        df_ds[i] = 1.0;
        df_ds[i + 1] = -1.0;
    }
    return s[i] - s[i + 1];
}

static double hs_surface_function(const HSStep* st, const HSSurface* sf, const double s[3],
                                  double gamma_p, double p_c, double df_ds[3], double* df_dgamma,
                                  double* df_dpc)
{
    if (df_dpc && sf->type != HS_CAP) *df_dpc = 0.0;
    if (df_dgamma && sf->type != HS_CONE) *df_dgamma = 0.0;

    switch (sf->type)
    {
    case HS_CONE:
        return hs_cone_function(st, s, gamma_p, df_ds, df_dgamma);
    case HS_MC:
        return hs_mc_function(st->prm, s, df_ds);
    case HS_CAP:
        return hs_cap_function(st->prm, sf->w, s, p_c, df_ds, df_dpc);
    default:
        return hs_corner_function(sf->corner, s, df_ds);
    }
}

/*
 * Plastic potential gradient of a shear surface, Eqs. 14-15:
 *   g_ij = (sigma_i - sigma_j)/2 - (sigma_i + sigma_j)/2 sin(psi)
 * Its plastic shear strain increment is d(eps_i^p) - d(eps_j^p) = dLambda_ij (Eq. 9).
 */
static void hs_pair_flow(double sin_psi, int i, int j, double n[3])
{
    n[0] = n[1] = n[2] = 0.0;
    n[i] = 0.5 - 0.5 * sin_psi;
    n[j] = -0.5 - 0.5 * sin_psi;
}

/*
 * Flow direction of the shear surfaces and of the corner equality. Away from the corners this is
 * g13; at a triaxial corner it is the average of g13 and g12 (compression) or g13 and g23
 * (extension), such that the multiplier is dLambda_13 + dLambda_12 (Eq. 15). The unequal split
 * between both surfaces is carried by the corner equality, whose flow is antisymmetric in the two
 * equal principal stresses and does not contribute to gamma_p.
 */
static void hs_shear_flow(const HSStep* st, const HSSurface* sf, double n[3])
{
    if (sf->type == HS_CORNER)
    {
        int i = (sf->corner == HS_CORNER_COMPRESSION) ? 1 : 0;
        n[0] = n[1] = n[2] = 0.0;
        n[i] = 0.5;
        n[i + 1] = -0.5;
        return;
    }

    double sin_psi = st->sin_psi_m;
    hs_pair_flow(sin_psi, 0, 2, n);
    if (sf->corner != HS_CORNER_NONE)
    {
        double n_corner[3];
        if (sf->corner == HS_CORNER_COMPRESSION)
            hs_pair_flow(sin_psi, 0, 1, n_corner);
        else
            hs_pair_flow(sin_psi, 1, 2, n_corner);
        for (int r = 0; r < 3; ++r) n[r] = 0.5 * (n[r] + n_corner[r]);
    }
}

/* ------------------------------------------------------------------ */
/* Return mapping in principal stress space                            */
/* ------------------------------------------------------------------ */

/*
 * Stress and hardening variables for given plastic multipliers of the active surfaces
 * (Eqs. 19, 21, 26 and 33):
 *   sigma   = sigma_tr - sum_k dLambda_k D n_k
 *   gamma_p = gamma_p0 + sum_{shear k} dLambda_k
 *   p_c     = p_c0 + H deps_v^pc = p_c0 + 2 H sum_{cap k} dLambda_k (p + a)
 * The flow directions of the shear (and corner) surfaces are constant during the step. The
 * associated cap flow n_k = N_k sigma + 2a/3 [1 1 1] is affine in sigma, so the stress follows from
 *   (I + sum_{cap k} dLambda_k D N_k) sigma
 *       = sigma_tr - sum_{other k} dLambda_k D n_k - sum_{cap k} dLambda_k 2a/3 D [1 1 1]
 * with N_k = 2/M^2 w_k w_k^T + 2/9 [1 1 1][1 1 1]^T.
 */
static int hs_evaluate_iterate(const HSStep* st, const HSSurface* surf, int n_act,
                               const double s_tr[3], const double dlambda[], HSIterate* it)
{
    const HSParams* prm = st->prm;
    double A[9], rhs[3], Dn[3], D_one[3];
    double dlambda_cap = 0.0;

    for (int r = 0; r < 3; ++r)
    {
        rhs[r] = s_tr[r];
        D_one[r] = st->D[3 * r] + st->D[3 * r + 1] + st->D[3 * r + 2];
        for (int col = 0; col < 3; ++col) A[3 * r + col] = (r == col) ? 1.0 : 0.0;
    }

    it->gamma_p = st->gamma_p0;
    for (int k = 0; k < n_act; ++k)
    {
        if (surf[k].type == HS_CAP)
        {
            double M2 = prm->M_cap * prm->M_cap;
            double D_w[3];
            hs_matvec3(st->D, surf[k].w, D_w);
            for (int r = 0; r < 3; ++r)
            {
                for (int col = 0; col < 3; ++col)
                    A[3 * r + col] +=
                        dlambda[k] * (2.0 / M2 * D_w[r] * surf[k].w[col] + 2.0 / 9.0 * D_one[r]);
                rhs[r] -= dlambda[k] * 2.0 * prm->a / 3.0 * D_one[r];
            }
            dlambda_cap += dlambda[k];
            continue;
        }
        hs_shear_flow(st, &surf[k], it->n[k]);
        hs_matvec3(st->D, it->n[k], Dn);
        for (int r = 0; r < 3; ++r) rhs[r] -= dlambda[k] * Dn[r];
        if (surf[k].type != HS_CORNER) it->gamma_p += dlambda[k];
    }

    if (!hs_invert3(A, it->A_inv)) return 0;
    hs_matvec3(it->A_inv, rhs, it->s);

    double p = (it->s[0] + it->s[1] + it->s[2]) / 3.0;
    it->p_c = st->p_c0 + 2.0 * st->H * dlambda_cap * (p + prm->a);
    for (int k = 0; k < n_act; ++k)
        if (surf[k].type == HS_CAP) hs_cap_function(prm, surf[k].w, it->s, it->p_c, it->n[k], NULL);

    for (int k = 0; k < n_act; ++k)
        it->f[k] = hs_surface_function(st, &surf[k], it->s, it->gamma_p, it->p_c, NULL, NULL, NULL);
    return 1;
}

static double hs_merit(const double f[], const double scale[], int n_act)
{
    double sum = 0.0;
    for (int k = 0; k < n_act; ++k) sum += (f[k] / scale[k]) * (f[k] / scale[k]);
    return sum;
}

static int hs_is_converged(const double f[], const double scale[], int n_act, double tol)
{
    for (int k = 0; k < n_act; ++k)
        if (fabs(f[k]) > tol * scale[k]) return 0;
    return 1;
}

/*
 * Newton iteration on the plastic multipliers of a fixed set of active surfaces, such that all
 * active yield functions vanish (consistency conditions of Eqs. 25-26 and 34). Returns 1 on
 * convergence with the multipliers in dlambda and the returned state in it.
 */
static int hs_solve_active_set(const HSStep* st, const HSSurface* surf, int n_act,
                               const double s_tr[3], const double scale[], double dlambda[],
                               HSIterate* it)
{
    const HSParams* prm = st->prm;
    HSIterate trial;
    double lambda_try[HS_MAX_ACTIVE];

    for (int k = 0; k < n_act; ++k) dlambda[k] = 0.0;
    if (!hs_evaluate_iterate(st, surf, n_act, s_tr, dlambda, it)) return 0;
    double merit = hs_merit(it->f, scale, n_act);

    for (int iter = 0; iter < HS_MAX_LOCAL_ITER; ++iter)
    {
        if (hs_is_converged(it->f, scale, n_act, HS_YIELD_TOL)) return 1;

        /* derivatives of the stress and of p_c w.r.t. the multipliers:
         * d sigma / d dLambda_l = -(I + sum_{cap k} dLambda_k D N_k)^-1 D n_l */
        double ds_dl[HS_MAX_ACTIVE][3], dpc_dl[HS_MAX_ACTIVE];
        double p = (it->s[0] + it->s[1] + it->s[2]) / 3.0;
        double dlambda_cap = 0.0;
        for (int k = 0; k < n_act; ++k)
            if (surf[k].type == HS_CAP) dlambda_cap += dlambda[k];

        for (int l = 0; l < n_act; ++l)
        {
            double Dn[3];
            hs_matvec3(st->D, it->n[l], Dn);
            hs_matvec3(it->A_inv, Dn, ds_dl[l]);
            for (int r = 0; r < 3; ++r) ds_dl[l][r] = -ds_dl[l][r];

            dpc_dl[l] = 2.0 * st->H *
                        (((surf[l].type == HS_CAP) ? (p + prm->a) : 0.0) +
                         dlambda_cap * (ds_dl[l][0] + ds_dl[l][1] + ds_dl[l][2]) / 3.0);
        }

        /* Jacobian d f_k / d dLambda_l */
        double J[HS_MAX_ACTIVE * HS_MAX_ACTIVE], rhs[HS_MAX_ACTIVE], step[HS_MAX_ACTIVE];
        for (int k = 0; k < n_act; ++k)
        {
            double df_ds[3], df_dgamma, df_dpc;
            hs_surface_function(st, &surf[k], it->s, it->gamma_p, it->p_c, df_ds, &df_dgamma,
                                &df_dpc);
            for (int l = 0; l < n_act; ++l)
            {
                double value = df_ds[0] * ds_dl[l][0] + df_ds[1] * ds_dl[l][1] +
                               df_ds[2] * ds_dl[l][2] + df_dpc * dpc_dl[l];
                /* d gamma_p / d dLambda_l = 1 for the shear surfaces */
                if (surf[l].type == HS_CONE || surf[l].type == HS_MC) value += df_dgamma;
                J[k * HS_MAX_ACTIVE + l] = value;
            }
            rhs[k] = -it->f[k];
        }
        if (!hs_solve_linear(n_act, J, rhs, step)) return 0;

        /* backtracking line search on the scaled residual */
        int accepted = 0;
        double t = 1.0;
        for (int ls = 0; ls < HS_MAX_LINE_SEARCH && !accepted; ++ls, t *= 0.5)
        {
            for (int k = 0; k < n_act; ++k) lambda_try[k] = dlambda[k] + t * step[k];
            if (!hs_evaluate_iterate(st, surf, n_act, s_tr, lambda_try, &trial)) continue;
            double merit_try = hs_merit(trial.f, scale, n_act);
            if (merit_try < merit)
            {
                accepted = 1;
                merit = merit_try;
                *it = trial;
                for (int k = 0; k < n_act; ++k) dlambda[k] = lambda_try[k];
            }
        }
        if (!accepted) break; /* no further decrease: stagnation at round-off level */
    }
    return hs_is_converged(it->f, scale, n_act, HS_YIELD_TOL_ACCEPT);
}

/*
 * Return mapping of the principal trial stress s_tr with an active set strategy:
 *  - shear: the hardening surface f13 or, when the result violates q <= qf, the Mohr-Coulomb
 *    surface (end of Sec. 3).
 *  - cap: when the trial stress or the returned stress lies outside the cap.
 *  - triaxial corners: when the principal stress order is lost, the return is made to the
 *    compression (sigma_2 = sigma_3) or extension (sigma_1 = sigma_2) corner. Both the shear
 *    surfaces (f13 and f12 or f23, Koiter, Eq. 26) and q~ of the cap have a corner there. As more
 *    surfaces meet than there are independent conditions, the stress and the hardening variables
 *    are unique but the split of the multipliers within each pair is not. The corner is therefore
 *    solved with the averaged flow of each pair plus a free antisymmetric plastic strain mu that
 *    enforces the equality of the principal stresses. It is admissible when mu can be carried by
 *    the pairs with non-negative multipliers.
 * Surfaces with a negative multiplier are removed from the active set.
 *
 * Returns 1 on success, with *plastic = 0 for an elastic step (s = s_tr). The active surfaces of
 * the returned stress (including the corner equality) are stored in active / n_active.
 */
static int hs_return_mapping(const HSStep* st, const double s_tr[3], double s[3], double* gamma_p,
                             double* p_c, int* plastic, int* on_failure,
                             HSSurface active[HS_MAX_ACTIVE], int* n_active)
{
    enum
    {
        SHEAR_NONE,
        SHEAR_CONE,
        SHEAR_MC
    };

    const HSParams* prm = st->prm;
    double w_ordered[3];
    hs_cap_weights(prm, HS_CORNER_NONE, w_ordered);

    /* fixed residual scales for the convergence checks */
    double stress_scale = fmax(fabs(s_tr[0]), fabs(s_tr[2])) + prm->a + 1.0e-3 * prm->p_ref;
    double scale_cone = prm->k_f / prm->Rf * stress_scale;
    double scale_mc = stress_scale;
    double scale_cap = (stress_scale + fabs(st->p_c0) + prm->a) * (stress_scale + fabs(st->p_c0) + prm->a);
    double order_tol = 1.0e-10 * stress_scale;

    int shear = SHEAR_NONE;
    HSCorner corner = HS_CORNER_NONE;
    int cap = 0;

    if (hs_cone_function(st, s_tr, st->gamma_p0, NULL, NULL) > HS_YIELD_TOL * scale_cone)
        shear = SHEAR_CONE;
    else if (hs_mc_function(prm, s_tr, NULL) > HS_YIELD_TOL * scale_mc)
        shear = SHEAR_MC;
    if (prm->use_cap &&
        hs_cap_function(prm, w_ordered, s_tr, st->p_c0, NULL, NULL) > HS_YIELD_TOL * scale_cap)
        cap = 1;

    for (int pass = 0; pass < HS_MAX_ACTIVE_SET_ITER; ++pass)
    {
        HSSurface surf[HS_MAX_ACTIVE];
        double scale[HS_MAX_ACTIVE], dlambda[HS_MAX_ACTIVE];
        HSIterate it;
        int n_act = 0;

        if (shear != SHEAR_NONE)
        {
            surf[n_act].type = (shear == SHEAR_CONE) ? HS_CONE : HS_MC;
            surf[n_act].corner = corner;
            scale[n_act++] = (shear == SHEAR_CONE) ? scale_cone : scale_mc;
        }
        if (cap)
        {
            surf[n_act].type = HS_CAP;
            surf[n_act].corner = corner;
            hs_cap_weights(prm, corner, surf[n_act].w);
            scale[n_act++] = scale_cap;
        }
        if (corner != HS_CORNER_NONE && n_act > 0)
        {
            surf[n_act].type = HS_CORNER;
            surf[n_act].corner = corner;
            scale[n_act++] = stress_scale;
        }

        if (!hs_solve_active_set(st, surf, n_act, s_tr, scale, dlambda, &it)) return 0;

        /* 1. surfaces with a negative multiplier are not active */
        int changed = 0;
        for (int k = 0; k < n_act; ++k)
        {
            if (surf[k].type == HS_CORNER || dlambda[k] >= 0.0) continue;
            changed = 1;
            if (surf[k].type == HS_CAP)
                cap = 0;
            else
                shear = SHEAR_NONE;
        }
        if (changed) continue;

        /* 2. at a corner, the antisymmetric plastic strain mu must be carried by the pairs:
         *    |mu| <= (1 +- sin(psi))/2 dLambda_shear + 2 |q~| / M^2 (2 alpha - 1 | 2 - alpha) dLambda_cap */
        if (corner != HS_CORNER_NONE && n_act > 0)
        {
            int compression = (corner == HS_CORNER_COMPRESSION);
            double capacity = 0.0, mu = 0.0;
            for (int k = 0; k < n_act; ++k)
            {
                if (surf[k].type == HS_CORNER)
                    mu = dlambda[k];
                else if (surf[k].type == HS_CAP)
                {
                    double q_tilde = surf[k].w[0] * it.s[0] + surf[k].w[1] * it.s[1] + surf[k].w[2] * it.s[2];
                    capacity += 2.0 * fabs(q_tilde) / (prm->M_cap * prm->M_cap) *
                                (compression ? 2.0 * prm->alpha - 1.0 : 2.0 - prm->alpha) * dlambda[k];
                }
                else
                {
                    double sin_psi = st->sin_psi_m;
                    capacity += 0.5 * (compression ? 1.0 + sin_psi : 1.0 - sin_psi) * dlambda[k];
                }
            }
            if (fabs(mu) > capacity * (1.0 + 1.0e-8) + SMALL_VALUE) return 0;
        }

        /* 3. principal stress order lost: return to the triaxial corner */
        if (n_act > 0 && corner == HS_CORNER_NONE)
        {
            if (it.s[1] < it.s[2] - order_tol)
            {
                corner = HS_CORNER_COMPRESSION;
                changed = 1;
            }
            else if (it.s[0] < it.s[1] - order_tol)
            {
                corner = HS_CORNER_EXTENSION;
                changed = 1;
            }
        }

        /* 4. failure criterion q <= qf, otherwise return to the Mohr-Coulomb surface (Sec. 3) */
        if (shear == SHEAR_CONE && hs_mc_function(prm, it.s, NULL) > HS_YIELD_TOL * scale_mc)
        {
            shear = SHEAR_MC;
            changed = 1;
        }

        /* 5. surfaces that are violated by the returned stress */
        if (shear == SHEAR_NONE)
        {
            if (hs_cone_function(st, it.s, it.gamma_p, NULL, NULL) > HS_YIELD_TOL * scale_cone)
            {
                shear = SHEAR_CONE;
                changed = 1;
            }
            else if (hs_mc_function(prm, it.s, NULL) > HS_YIELD_TOL * scale_mc)
            {
                shear = SHEAR_MC;
                changed = 1;
            }
        }
        if (prm->use_cap && !cap &&
            hs_cap_function(prm, w_ordered, it.s, it.p_c, NULL, NULL) > HS_YIELD_TOL * scale_cap)
        {
            cap = 1;
            changed = 1;
        }

        if (!changed)
        {
            for (int r = 0; r < 3; ++r) s[r] = it.s[r];
            *gamma_p = it.gamma_p;
            *p_c = it.p_c;
            *plastic = (n_act > 0) ? 1 : 0;
            *on_failure = (shear == SHEAR_MC) ? 1 : 0;
            *n_active = n_act;
            for (int k = 0; k < n_act; ++k) active[k] = surf[k];
            return 1;
        }
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Single step and integration with automatic sub-stepping             */
/* ------------------------------------------------------------------ */

/*
 * Return to the apex of the Mohr-Coulomb pyramid, sigma_1 = sigma_2 = sigma_3 = -a, for trial
 * stresses in tension that cannot be returned to its faces or edges. The paper does not treat
 * tension; this keeps the stress integration defined there. gamma_p accumulates the plastic shear
 * strain eps_1^p - eps_3^p.
 */
static void hs_apex_return(const HSStep* st, const double s_tr[3], double s[3], double* gamma_p,
                           double* p_c)
{
    double G = st->Eur / (2.0 * (1.0 + st->prm->nu));
    for (int r = 0; r < 3; ++r) s[r] = -st->prm->a;
    *gamma_p = st->gamma_p0 + (s_tr[0] - s_tr[2]) / (2.0 * G);
    *p_c = st->p_c0;
}

/* Principal values v (tensor with principal directions Q) in Voigt notation with engineering
 * shear components, the form of strain-like vectors and of stress gradients. */
static void hs_principal_to_voigt(const double v[3], double Q[3][3], double out[VOIGTSIZE_3D])
{
    for (int i = 0; i < VOIGTSIZE_3D; ++i) out[i] = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        out[XX] += v[k] * Q[0][k] * Q[0][k];
        out[YY] += v[k] * Q[1][k] * Q[1][k];
        out[ZZ] += v[k] * Q[2][k] * Q[2][k];
        out[XY] += 2.0 * v[k] * Q[0][k] * Q[1][k];
        out[YZ] += 2.0 * v[k] * Q[1][k] * Q[2][k];
        out[XZ] += 2.0 * v[k] * Q[0][k] * Q[2][k];
    }
}

/*
 * Elasto-plastic (continuum) tangent for the given active surfaces of the returned stress (Koiter):
 *   D_ep = De - De N (G^T De N + H)^-1 G^T De
 * with N the flow directions, G the yield function gradients and H the hardening moduli
 * H_kl = -(df_k/dgamma_p dgamma_p/dLambda_l + df_k/dp_c dp_c/dLambda_l). At a triaxial corner the
 * gradients of the shear and cap surfaces are averaged over both equal principal stresses, like
 * the flow directions, so that the tangent does not depend on the arbitrary principal directions
 * within the corner. The corner equality, when present, removes the stiffness of the corner mode:
 * the difference of the two equal principal stresses and the shear in their plane stay zero.
 * Returns 0 (with ddsdde = De) if the Koiter matrix is singular.
 */
static int hs_koiter_tangent(const HSStep* st, const HSSurface* active, int n_active,
                             const double s[3], double gamma_p, double p_c, double Q[3][3],
                             double ddsdde[VOIGTSIZE_3D * VOIGTSIZE_3D])
{
    const HSParams* prm = st->prm;
    double g6[HS_MAX_ACTIVE][VOIGTSIZE_3D], De_n[HS_MAX_ACTIVE][VOIGTSIZE_3D];
    double De_g[HS_MAX_ACTIVE][VOIGTSIZE_3D];
    double df_dgamma[HS_MAX_ACTIVE], df_dpc[HS_MAX_ACTIVE];
    double M[HS_MAX_ACTIVE * HS_MAX_ACTIVE], M_inv[HS_MAX_ACTIVE * HS_MAX_ACTIVE];
    double De[VOIGTSIZE_3D * VOIGTSIZE_3D];

    calculate_elastic_stiffness_matrix_3d(st->Eur, prm->nu, De);
    copy_array(De, VOIGTSIZE_3D * VOIGTSIZE_3D, ddsdde);
    if (n_active == 0) return 1;

    for (int k = 0; k < n_active; ++k)
    {
        double n[3], g[3], n6[VOIGTSIZE_3D];
        hs_surface_function(st, &active[k], s, gamma_p, p_c, g, &df_dgamma[k], &df_dpc[k]);
        if (active[k].type == HS_CAP)
            for (int r = 0; r < 3; ++r) n[r] = g[r]; /* associated flow, Eq. 30 */
        else
            hs_shear_flow(st, &active[k], n);

        if (active[k].type != HS_CORNER && active[k].corner == HS_CORNER_COMPRESSION)
            g[1] = g[2] = 0.5 * (g[1] + g[2]);
        else if (active[k].type != HS_CORNER && active[k].corner == HS_CORNER_EXTENSION)
            g[0] = g[1] = 0.5 * (g[0] + g[1]);

        hs_principal_to_voigt(n, Q, n6);
        hs_principal_to_voigt(g, Q, g6[k]);
        matrix_vector_multiply(De, n6, VOIGTSIZE_3D, De_n[k]);
        matrix_vector_multiply(De, g6[k], VOIGTSIZE_3D, De_g[k]);
    }

    /* M = G^T De N + H, inverted column by column */
    double p = (s[0] + s[1] + s[2]) / 3.0;
    for (int k = 0; k < n_active; ++k)
    {
        for (int l = 0; l < n_active; ++l)
        {
            double value = vector_dot_product(g6[k], De_n[l], VOIGTSIZE_3D);
            if (active[l].type == HS_CAP)
                value -= df_dpc[k] * 2.0 * st->H * (p + prm->a);
            else if (active[l].type != HS_CORNER)
                value -= df_dgamma[k];
            M[k * HS_MAX_ACTIVE + l] = value;
        }
    }
    for (int col = 0; col < n_active; ++col)
    {
        double J[HS_MAX_ACTIVE * HS_MAX_ACTIVE], e[HS_MAX_ACTIVE], x[HS_MAX_ACTIVE];
        for (int k = 0; k < HS_MAX_ACTIVE * HS_MAX_ACTIVE; ++k) J[k] = M[k];
        for (int k = 0; k < n_active; ++k) e[k] = (k == col) ? 1.0 : 0.0;
        if (!hs_solve_linear(n_active, J, e, x)) return 0;
        for (int k = 0; k < n_active; ++k) M_inv[k * HS_MAX_ACTIVE + col] = x[k];
    }

    for (int i = 0; i < VOIGTSIZE_3D; ++i)
        for (int j = 0; j < VOIGTSIZE_3D; ++j)
            for (int k = 0; k < n_active; ++k)
                for (int l = 0; l < n_active; ++l)
                    ddsdde[i * VOIGTSIZE_3D + j] -= De_n[k][i] * M_inv[k * HS_MAX_ACTIVE + l] * De_g[l][j];

    /* shear in the plane of the two equal principal stresses at a corner: G (s_a - s_b) /
     * (s_a^tr - s_b^tr) with s_a = s_b, i.e. no stiffness */
    for (int k = 0; k < n_active; ++k)
    {
        if (active[k].type != HS_CORNER) continue;
        int a = (active[k].corner == HS_CORNER_COMPRESSION) ? 1 : 0, b = a + 1;
        double S[VOIGTSIZE_3D];
        S[XX] = Q[0][a] * Q[0][b];
        S[YY] = Q[1][a] * Q[1][b];
        S[ZZ] = Q[2][a] * Q[2][b];
        S[XY] = 0.5 * (Q[0][a] * Q[1][b] + Q[0][b] * Q[1][a]);
        S[YZ] = 0.5 * (Q[1][a] * Q[2][b] + Q[1][b] * Q[2][a]);
        S[XZ] = 0.5 * (Q[0][a] * Q[2][b] + Q[0][b] * Q[2][a]);
        double G = st->Eur / (2.0 * (1.0 + prm->nu));
        for (int i = 0; i < VOIGTSIZE_3D; ++i)
            for (int j = 0; j < VOIGTSIZE_3D; ++j) ddsdde[i * VOIGTSIZE_3D + j] -= 4.0 * G * S[i] * S[j];
    }
    return 1;
}

/*
 * Tangent returned to the host (DDSDDE). Away from the triaxial corners this is the Koiter tangent,
 * which matches the derivative of the stress integration. At a corner the exact tangent has no
 * stiffness in the corner mode, which can make the global stiffness matrix singular (e.g. a 3D K0
 * state with equal horizontal stresses). A fraction HS_CORNER_STIFFNESS_FRACTION of the stiffness
 * of that mode is therefore kept, by blending with the tangent without the corner equality.
 */
static void hs_elastoplastic_tangent(const HSStep* st, const HSSurface* active, int n_active,
                                     const double s[3], double gamma_p, double p_c, double Q[3][3],
                                     double ddsdde[VOIGTSIZE_3D * VOIGTSIZE_3D])
{
    HSSurface smooth[HS_MAX_ACTIVE];
    int n_smooth = 0;
    for (int k = 0; k < n_active; ++k)
        if (active[k].type != HS_CORNER) smooth[n_smooth++] = active[k];

    if (!hs_koiter_tangent(st, active, n_active, s, gamma_p, p_c, Q, ddsdde) || n_smooth == n_active)
    {
        if (n_smooth != n_active) hs_koiter_tangent(st, smooth, n_smooth, s, gamma_p, p_c, Q, ddsdde);
        return;
    }

    double D_smooth[VOIGTSIZE_3D * VOIGTSIZE_3D];
    hs_koiter_tangent(st, smooth, n_smooth, s, gamma_p, p_c, Q, D_smooth);
    for (int i = 0; i < VOIGTSIZE_3D * VOIGTSIZE_3D; ++i)
        ddsdde[i] = (1.0 - HS_CORNER_STIFFNESS_FRACTION) * ddsdde[i] + HS_CORNER_STIFFNESS_FRACTION * D_smooth[i];
}

/*
 * One (sub-)step: elastic predictor (Eq. 20) and return mapping on the principal trial stresses.
 * Elasticity and yield surfaces are isotropic, so the returned stress keeps the principal
 * directions of the trial stress. Returns 1 on success; the stress and the state variables are
 * then updated in place and, if ddsdde is given, the elasto-plastic tangent at the returned stress
 * is stored in it.
 */
static int hs_single_step(double stress[VOIGTSIZE_3D], double* gamma_p, double* p_c,
                          int* at_failure, const double deps[VOIGTSIZE_3D], double void_ratio,
                          const HSParams* prm, double* ddsdde)
{
    HSStep st;
    HSSurface active[HS_MAX_ACTIVE];
    double s0[3], s_tr[3], s[3];
    double Q0[3][3], Q[3][3];
    double Ce[VOIGTSIZE_3D * VOIGTSIZE_3D];
    double delta_sigma[VOIGTSIZE_3D];
    double stress_trial[VOIGTSIZE_3D];
    int plastic = 0, on_failure = 0, n_active = 0;

    /* stiffness and dilatancy at the start of the step (Euler explicit, Sec. 3) */
    calculate_principal_system(stress, s0, Q0);
    double factor = hs_stiffness_factor(prm, s0[2]);
    st.prm = prm;
    st.Eur = prm->Eur_ref * factor;
    st.Ei = prm->Ei_ref * factor;
    hs_principal_elastic_matrix(st.Eur, prm->nu, st.D);
    st.sin_psi_m = hs_sin_psi_mobilised(prm, s0, void_ratio);
    st.H = prm->use_cap ? hs_cap_modulus(prm, *p_c) : 0.0;
    st.gamma_p0 = *gamma_p;
    st.p_c0 = *p_c;

    /* elastic predictor in the global frame, Eq. 20 */
    calculate_elastic_stiffness_matrix_3d(st.Eur, prm->nu, Ce);
    matrix_vector_multiply(Ce, deps, VOIGTSIZE_3D, delta_sigma);
    add_vectors(stress, delta_sigma, VOIGTSIZE_3D, stress_trial);

    /* return mapping in principal stress space */
    calculate_principal_system(stress_trial, s_tr, Q);
    int ok = hs_return_mapping(&st, s_tr, s, gamma_p, p_c, &plastic, &on_failure, active, &n_active);

    /* a trial stress in tension beyond the apex has no valid return to the faces or edges */
    int at_apex = 0;
    double apex_tol = 1.0e-10 * (fabs(s_tr[0]) + prm->a + prm->p_ref);
    if (!ok || (plastic && s[2] + prm->a < -apex_tol))
    {
        if (s_tr[2] + prm->a >= 0.0) return 0;
        hs_apex_return(&st, s_tr, s, gamma_p, p_c);
        plastic = 1;
        on_failure = 1;
        n_active = 0;
        at_apex = 1;
    }

    if (plastic)
        calculate_stress_from_principal_system(s, Q, stress);
    else
        copy_array(stress_trial, VOIGTSIZE_3D, stress);
    *at_failure = on_failure;

    if (ddsdde)
    {
        calculate_elastic_stiffness_matrix_3d(st.Eur, prm->nu, ddsdde);
        //hs_elastoplastic_tangent(&st, active, plastic ? n_active : 0, s, *gamma_p, *p_c, Q, ddsdde);
        /* at the apex the stress does not change under further loading: keep only a fraction of
         * the elastic stiffness, as for the corner mode, to keep the global matrix regular */
        if (at_apex)
            for (int i = 0; i < VOIGTSIZE_3D * VOIGTSIZE_3D; ++i) ddsdde[i] *= HS_CORNER_STIFFNESS_FRACTION;
    }
    return 1;
}

/*
 * Integrates the strain increment with sub-steps. On success the stress and state variables are
 * updated and ddsdde holds the elasto-plastic tangent of the last sub-step.
 *
 * The sub-steps have a fixed size (HS_SUBSTEP_STRESS_FRACTION of sigma_3 + a in terms of the
 * elastic stress increment) and the first sub-step takes the remaining fraction. The integrated
 * stress is then a continuous function of the strain increment. With a rounded-up number of equal
 * sub-steps, the explicit integration error jumps whenever the number changes, and those jumps put
 * a floor under the residual of the global Newton iteration.
 */
static int hs_integrate(double stress[VOIGTSIZE_3D], double* gamma_p, double* p_c,
                        int* at_failure, const double dstrain[VOIGTSIZE_3D], double eps_v0,
                        const HSParams* prm, double* ddsdde)
{
    double Ce[VOIGTSIZE_3D * VOIGTSIZE_3D];
    double delta_sigma[VOIGTSIZE_3D];

    /* (real) number of sub-steps from the size of the elastic stress increment */
    double sigma_3 = hs_minor_principal_stress(stress);
    double Eur = prm->Eur_ref * hs_stiffness_factor(prm, sigma_3);
    calculate_elastic_stiffness_matrix_3d(Eur, prm->nu, Ce);
    matrix_vector_multiply(Ce, dstrain, VOIGTSIZE_3D, delta_sigma);

    double increment = fmax(hs_q(delta_sigma), fabs(hs_mean_stress(delta_sigma)));
    double stress_scale = fmax(sigma_3 + prm->a, HS_MIN_STRESS_RATIO * (prm->p_ref + prm->a));
    double n_real = increment / (HS_SUBSTEP_STRESS_FRACTION * stress_scale);
    if (!(n_real >= 1.0)) n_real = 1.0;
    if (n_real > HS_MAX_SUBSTEPS) n_real = HS_MAX_SUBSTEPS;

    double deps_v = dstrain[XX] + dstrain[YY] + dstrain[ZZ];

    for (int trial = 0; trial < HS_MAX_SUBSTEP_TRIALS; ++trial)
    {
        double sigma[VOIGTSIZE_3D];
        double deps_step[VOIGTSIZE_3D];
        copy_array(stress, VOIGTSIZE_3D, sigma);
        double gp = *gamma_p;
        double pc = *p_c;
        int fail_flag = *at_failure;

        /* n_full sub-steps of size 1 / n_real, preceded by the remainder */
        int n_full = (int)floor(n_real);
        double remainder = n_real - (double)n_full;
        int n_sub = n_full + ((remainder > 0.0) ? 1 : 0);
        double fraction_done = 0.0;

        int ok = 1;
        for (int step = 0; step < n_sub && ok; ++step)
        {
            double fraction = ((step == 0 && remainder > 0.0) ? remainder : 1.0) / n_real;
            for (int i = 0; i < VOIGTSIZE_3D; ++i) deps_step[i] = dstrain[i] * fraction;

            /* void ratio at the start of the sub-step, Eq. 39 (eps_v = 0 at e = e0) */
            double eps_v = eps_v0 + deps_v * fraction_done;
            double void_ratio = (1.0 + prm->e0) * exp(-eps_v) - 1.0;
            ok = hs_single_step(sigma, &gp, &pc, &fail_flag, deps_step, void_ratio, prm,
                                (step == n_sub - 1) ? ddsdde : NULL);
            fraction_done += fraction;
        }

        if (ok)
        {
            copy_array(sigma, VOIGTSIZE_3D, stress);
            *gamma_p = gp;
            *p_c = pc;
            *at_failure = fail_flag;
            return 1;
        }

        if (n_real >= HS_MAX_SUBSTEPS) break;
        n_real = (2.0 * n_real > HS_MAX_SUBSTEPS) ? HS_MAX_SUBSTEPS : 2.0 * n_real; /* refine, retry */
    }

    return 0;
}

/* ------------------------------------------------------------------ */
/* Parameters and initial state                                        */
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
    if (props[4] <= 0.0 || props[4] >= 90.0) { fprintf(stderr, "UMAT Error: phi in (0, 90).\n"); n_errors++; }
    if (props[5] > props[4]) { fprintf(stderr, "UMAT Error: psi must not exceed phi.\n"); n_errors++; }
    if (props[6] <= 0.0) { fprintf(stderr, "UMAT Error: p_ref must be positive.\n"); n_errors++; }
    if (props[7] <= 0.0 || props[7] >= 1.0) { fprintf(stderr, "UMAT Error: Rf in (0, 1).\n"); n_errors++; }
    else if (props[1] <= 2.0 * props[0] / (2.0 - props[7]))
    {
        fprintf(stderr, "UMAT Error: Eur_ref must exceed Ei_ref = 2 E50_ref / (2 - Rf).\n");
        n_errors++;
    }
    if (props[8] < 0.0 || props[8] >= 0.5) { fprintf(stderr, "UMAT Error: nu in [0, 0.5).\n"); n_errors++; }
    if (props[9] < 0.0) { fprintf(stderr, "UMAT Error: M_cap must be non-negative (0 disables the cap).\n"); n_errors++; }
    if (props[9] > 0.0 && props[10] <= 1.0) { fprintf(stderr, "UMAT Error: K_ratio must be > 1.\n"); n_errors++; }
    if (nprops > 12 && props[12] > 0.0 && props[11] <= 0.0)
    {
        fprintf(stderr, "UMAT Error: e0 must be positive when the dilatancy cut-off is used.\n");
        n_errors++;
    }
    return (n_errors > 0) ? 1 : 0;
}

static void hs_set_parameters(int nprops, const double* props, HSParams* prm)
{
    prm->E50_ref = props[0];
    prm->Eur_ref = props[1];
    prm->m = props[2];
    prm->c = props[3];
    prm->phi = props[4] * PI / 180.0;
    prm->psi = props[5] * PI / 180.0;
    prm->p_ref = props[6];
    prm->Rf = props[7];
    prm->nu = props[8];
    prm->M_cap = props[9];
    prm->K_ratio = props[10];
    prm->e0 = (nprops > 11) ? props[11] : 0.0;
    prm->e_cv = (nprops > 12) ? props[12] : 0.0;

    prm->sin_phi = sin(prm->phi);
    prm->cos_phi = cos(prm->phi);
    prm->sin_psi = sin(prm->psi);
    prm->sin_phi_cv = (prm->sin_phi - prm->sin_psi) / (1.0 - prm->sin_phi * prm->sin_psi);
    prm->a = prm->c * prm->cos_phi / prm->sin_phi;
    prm->k_f = 2.0 * prm->sin_phi / (1.0 - prm->sin_phi);
    prm->alpha = (3.0 + prm->sin_phi) / (3.0 - prm->sin_phi);

    /* Hyperbola (Eq. 1) with E50 the secant stiffness at q = qf / 2 (Sec. 2.1). */
    prm->Ei_ref = 2.0 * prm->E50_ref / (2.0 - prm->Rf);

    prm->use_cap = (prm->M_cap > 0.0) ? 1 : 0;
    prm->use_cutoff = (prm->e_cv > 0.0) ? 1 : 0;

    /* H = Ks Kc / (Ks - Kc) = Ks / (K_ratio - 1), Eq. 32 */
    prm->H_ref = 0.0;
    if (prm->use_cap)
    {
        double Ks_ref = prm->Eur_ref / (3.0 * (1.0 - 2.0 * prm->nu));
        prm->H_ref = Ks_ref / (prm->K_ratio - 1.0);
    }
}

/*
 * First call: initialise the hardening parameters such that the initial stress lies on the shear
 * hardening surface (Eq. 8) and on the cap (Eq. 27), i.e. a normally consolidated state.
 */
static void hs_initialise_state(const HSParams* prm, const double stress[VOIGTSIZE_3D],
                                double* gamma_p, double* p_c)
{
    double s[3];
    double Q[3][3];
    calculate_principal_system(stress, s, Q);

    double p = (s[0] + s[1] + s[2]) / 3.0;
    double q_tilde = s[0] + (prm->alpha - 1.0) * s[1] - prm->alpha * s[2];
    double q_term = prm->use_cap ? q_tilde / prm->M_cap : 0.0;
    double p_c0 = sqrt(q_term * q_term + (p + prm->a) * (p + prm->a)) - prm->a;
    double p_c_min = 1.0e-2 * prm->p_ref;
    *p_c = (p_c0 > p_c_min) ? p_c0 : p_c_min;

    double qf = prm->k_f * (s[2] + prm->a);
    if (qf > 0.0)
    {
        double factor = hs_stiffness_factor(prm, s[2]);
        double q = s[0] - s[2];
        if (q > qf) q = qf;
        double gamma_0 = hs_hyperbolic_gamma_p(q, qf / prm->Rf, prm->Ei_ref * factor,
                                               prm->Eur_ref * factor);
        if (gamma_0 > *gamma_p) *gamma_p = gamma_0;
    }
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
    (void)TIME; (void)DTIME; (void)TEMP; (void)DTEMP; (void)PREDEF;
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
    hs_set_parameters(*NPROPS, PROPS, &prm);

    /* --- read state --- */
    double gamma_p = STATEV[0];
    double p_c = STATEV[1];
    int at_failure = (STATEV[2] > 0.5) ? 1 : 0;

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
    double eps_v0 = -(STRAN[XX] + STRAN[YY] + STRAN[ZZ]);

    if (p_c <= 0.0) hs_initialise_state(&prm, stress, &gamma_p, &p_c);

    /* --- integrate the stress --- */
    int converged = hs_integrate(stress, &gamma_p, &p_c, &at_failure, dstrain, eps_v0, &prm, DDSDDE);

    if (!converged)
    {
//        fprintf(stderr, "UMAT Warning: Hardening Soil return mapping did not converge; "
//                        "requesting a smaller time increment.\n");
//#ifdef HS_DEBUG
//        fprintf(stderr, "[hs] stress = [%g %g %g %g %g %g], gamma_p = %g, p_c = %g\n", STRESS[XX],
//                STRESS[YY], STRESS[ZZ], STRESS[XY], STRESS[YZ], STRESS[XZ], STATEV[0], STATEV[1]);
//#endif
        /* keep the stress and state at the start of the increment, with the elastic tangent */
        double factor = hs_stiffness_factor(&prm, hs_minor_principal_stress(stress));
        calculate_elastic_stiffness_matrix_3d(prm.Eur_ref * factor, prm.nu, DDSDDE);
        if (PNEWDT) *PNEWDT = 0.5;
        return;
    }

    /* --- write results --- */
    for (int i = 0; i < VOIGTSIZE_3D; ++i) STRESS[i] = -stress[i];
    STATEV[0] = gamma_p;
    STATEV[1] = p_c;
    STATEV[2] = (double)at_failure;

    if (SPD) *SPD = 0.0;
    return;
}
