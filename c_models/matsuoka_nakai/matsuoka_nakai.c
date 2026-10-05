/**
 * @file matsuoka_nakai.c
 * @brief UMAT for a linear elastic, perfectly plastic model with the Matsuoka-Nakai yield surface
 *        \cite Matsuoka_1974 in the formulation of \cite Lagioia_2016.
 *
 * Model
 * -----
 *  - Isotropic linear elasticity (E, nu).
 *  - Yield function (tension positive, see yield_surfaces/matsuoka_nakai_surface.h)
 *        f = M p - K + J alpha cos(acos(beta sin(3 theta)) / 3)
 *    which is the Matsuoka-Nakai criterion I1 I2 / I3 = 9 + 8 tan^2(phi) of the stresses shifted
 *    by c cot(phi). It passes through all six corners of the Mohr-Coulomb pyramid with the same c
 *    and phi; for phi = 0 it is the von Mises cylinder through the corners of the Tresca hexagon.
 *  - Non-associated flow with the same function with the dilation angle psi as plastic potential.
 *
 * Integration
 * -----------
 *  - Implicit (backward Euler) return mapping: Newton iteration on the stress and the plastic
 *    multiplier, with the second derivative of the plastic potential and a line search.
 *  - DDSDDE is the consistent (algorithmic) tangent of the return mapping.
 *  - Return to the apex of the cone, p = c cot(phi), when the plastic strain of that return is a
 *    flow direction of the plastic potential at its apex (beyond the apex for psi > 0; for
 *    psi < 0, on the compressive side, only if there is no return to the cone), and for all trial
 *    stresses beyond the apex when psi <= 0, which cannot be returned to the cone. With psi < 0 the
 *    return mapping can have no solution at all; in the smallest sub-steps the stress is then
 *    returned to the apex (see integrate_step). The stress does not change under further loading
 *    at the apex; a fraction MN_APEX_STIFFNESS_FRACTION of the elastic stiffness is kept in DDSDDE
 *    so that the global stiffness matrix stays regular.
 *  - If the Newton iteration does not converge, the strain increment is divided into sub-steps.
 *    If that fails as well, the stress is not updated and PNEWDT = 0.5 is requested.
 *
 * Conventions
 * -----------
 *  - Voigt ordering: [xx, yy, zz, xy, yz, xz] (see globals.h), engineering shear strains.
 *  - Tension positive (Abaqus convention), like the other models in this library.
 *
 * Material properties (PROPS)
 * ---------------------------
 *   [0]  E   - Young's modulus                       [stress]
 *   [1]  nu  - Poisson's ratio                       [-]
 *   [2]  c   - cohesion                              [stress]
 *   [3]  phi - friction angle, 0 <= phi < 90         [degrees]
 *   [4]  psi - dilation angle, -90 < psi < 90        [degrees]
 *
 * State variables (STATEV)
 * ------------------------
 *   [0]  state of the increment: 0 elastic, 1 plastic (cone), 2 plastic (apex)
 *
 * Energies
 * --------
 *   SSE is updated with the elastic strain energy increment, SPD with the plastic dissipation
 *   STRESS : dEps_p of the increment.
 */

#include <math.h>
#include <stdio.h>

#include "../elastic_laws/hookes_law.h"
#include "../globals.h"
#include "../stress_utils.h"
#include "../utils.h"
#include "../yield_surfaces/matsuoka_nakai_surface.h"

// Define necessary calling conventions and export macros (adjust for your compiler/system)
// For MSVC on Windows:
#if defined(_WIN32) || defined(_WIN64)
#define UMAT_EXPORT __declspec(dllexport)
#define UMAT_CALLCONV __stdcall  // Abaqus often uses stdcall
#else
// For GCC/Clang on Linux/macOS (usually no special decoration needed)
#define UMAT_EXPORT
#define UMAT_CALLCONV
#endif

#define MN_MAX_ITERATIONS 50    // maximum number of Newton iterations of the return mapping
#define MN_MAX_LINE_SEARCH 12   // maximum number of step halvings in the line search
#define MN_MAX_SUBSTEP_LEVEL 6  // at most 2^MN_MAX_SUBSTEP_LEVEL sub-steps

// Tolerance on the yield function and the residual of the return mapping, relative to the stress
// level |p| + J + K of the trial stress
#define MN_YIELD_TOL 1.0e-10

// Fraction of the elastic stiffness kept in DDSDDE at the apex, where the exact tangent is zero
#ifndef MN_APEX_STIFFNESS_FRACTION
#define MN_APEX_STIFFNESS_FRACTION 1.0e-2
#endif

// state of an increment, stored in STATEV[0]
#define MN_STATE_ELASTIC 0
#define MN_STATE_PLASTIC 1
#define MN_STATE_APEX 2

#define VOIGTSIZE_3D_SQ (VOIGTSIZE_3D * VOIGTSIZE_3D)

/**
 * @brief Elastic stiffness, compliance and yield surface constants.
 */
typedef struct
{
    double Ce[VOIGTSIZE_3D_SQ];      ///< elastic stiffness matrix (row-major)
    double Ce_inv[VOIGTSIZE_3D_SQ];  ///< elastic compliance matrix (row-major)
    double bulk_modulus;             ///< elastic bulk modulus
    double shear_modulus;            ///< elastic shear modulus
    MatsuokaNakaiConstants yield;      ///< constants of the yield function (phi, c)
    MatsuokaNakaiConstants potential;  ///< constants of the plastic potential (psi)
    int has_apex;                      ///< 1 if the cone has an apex (phi > 0)
    double p_apex;                     ///< mean stress at the apex, c cot(phi) (tension positive)
} MNModel;

/**
 * @brief Function to check the validity of material properties.
 *
 * @param[in]  NPROPS Number of properties.
 * @param[in]  PROPS Pointer to the array of material properties.
 * @return 0 if properties are valid, 1 otherwise.
 */
int check_properties(const int NPROPS, const double* PROPS);

/**
 * @brief Elastic compliance matrix (inverse of the Hooke stiffness) with engineering shear
 * strains.
 */
static void calculate_elastic_compliance_matrix_3d(const double E, const double nu,
                                                   double compliance[VOIGTSIZE_3D_SQ])
{
    const double G = E / (2.0 * (1.0 + nu));
    for (int i = 0; i < VOIGTSIZE_3D_SQ; ++i) compliance[i] = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j) compliance[i * VOIGTSIZE_3D + j] = (i == j) ? 1.0 / E : -nu / E;
        compliance[(i + 3) * VOIGTSIZE_3D + (i + 3)] = 1.0 / G;
    }
}

/**
 * @brief Inverts a 6x6 matrix by Gauss-Jordan elimination with partial pivoting.
 *
 * @return 1 on success, 0 if the matrix is singular.
 */
static int invert_matrix_6x6(const double matrix[VOIGTSIZE_3D_SQ],
                             double inverse[VOIGTSIZE_3D_SQ])
{
    const int n = VOIGTSIZE_3D;
    double a[VOIGTSIZE_3D_SQ];
    double max_entry = 0.0;
    for (int i = 0; i < n * n; ++i)
    {
        a[i] = matrix[i];
        inverse[i] = (i % (n + 1) == 0) ? 1.0 : 0.0;
        if (fabs(a[i]) > max_entry) max_entry = fabs(a[i]);
    }
    if (max_entry == 0.0) return 0;

    for (int col = 0; col < n; ++col)
    {
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (fabs(a[r * n + col]) > fabs(a[piv * n + col])) piv = r;
        if (fabs(a[piv * n + col]) < 1.0e-14 * max_entry) return 0;

        if (piv != col)
        {
            for (int k = 0; k < n; ++k)
            {
                double tmp = a[col * n + k];
                a[col * n + k] = a[piv * n + k];
                a[piv * n + k] = tmp;
                tmp = inverse[col * n + k];
                inverse[col * n + k] = inverse[piv * n + k];
                inverse[piv * n + k] = tmp;
            }
        }

        const double pivot = a[col * n + col];
        for (int k = 0; k < n; ++k)
        {
            a[col * n + k] /= pivot;
            inverse[col * n + k] /= pivot;
        }
        for (int r = 0; r < n; ++r)
        {
            if (r == col) continue;
            const double factor = a[r * n + col];
            if (factor == 0.0) continue;
            for (int k = 0; k < n; ++k)
            {
                a[r * n + k] -= factor * a[col * n + k];
                inverse[r * n + k] -= factor * inverse[col * n + k];
            }
        }
    }
    return 1;
}

/**
 * @brief Quantities of the return mapping evaluated at a stress state.
 */
typedef struct
{
    double p;                        ///< mean stress
    double J;                        ///< sqrt(J2)
    double f;                        ///< yield function
    double grad_f[VOIGTSIZE_3D];     ///< gradient of the yield function
    double grad_g[VOIGTSIZE_3D];     ///< gradient of the plastic potential
    double hess_g[VOIGTSIZE_3D_SQ];  ///< second derivative of the plastic potential
    double residual[VOIGTSIZE_3D];   ///< stress - stress_trial + delta_gamma Ce grad_g
} MNEvaluation;

static void evaluate(const MNModel* model, const double stress[VOIGTSIZE_3D],
                     const double stress_trial[VOIGTSIZE_3D], const double delta_gamma,
                     MNEvaluation* ev)
{
    double theta, j2, j3, s_dev[VOIGTSIZE_3D];
    calculate_stress_invariants_3d(stress, &ev->p, &ev->J, &theta, &j2, &j3, s_dev);
    calculate_yield_function(ev->p, theta, ev->J, model->yield, &ev->f);
    calculate_yield_derivatives(s_dev, j2, j3, model->yield, ev->grad_f, NULL);
    calculate_yield_derivatives(s_dev, j2, j3, model->potential, ev->grad_g, ev->hess_g);

    double Ce_grad_g[VOIGTSIZE_3D];
    matrix_vector_multiply(model->Ce, ev->grad_g, VOIGTSIZE_3D, Ce_grad_g);
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
        ev->residual[i] = stress[i] - stress_trial[i] + delta_gamma * Ce_grad_g[i];
}

/**
 * @brief Scaled squared norm of the residual of the return mapping (residual stress and yield
 * function).
 */
static double merit_function(const MNEvaluation* ev, const double scale)
{
    double merit = ev->f * ev->f + vector_dot_product(ev->residual, ev->residual, VOIGTSIZE_3D);
    return merit / (scale * scale);
}

/**
 * @brief Inverse of Ce^-1 + delta_gamma d2g/dsigma2, the matrix of the linearised return mapping.
 */
static int calculate_algorithmic_stiffness(const MNModel* model, const MNEvaluation* ev,
                                           const double delta_gamma,
                                           double Xi[VOIGTSIZE_3D_SQ])
{
    double A[VOIGTSIZE_3D_SQ];
    for (int i = 0; i < VOIGTSIZE_3D_SQ; ++i) A[i] = model->Ce_inv[i] + delta_gamma * ev->hess_g[i];
    return invert_matrix_6x6(A, Xi);
}

/**
 * @brief Implicit return of a trial stress to the smooth part of the yield surface.
 *
 * Solves  Ce^-1 (sigma - sigma_trial) + delta_gamma dg/dsigma(sigma) = 0  and  f(sigma) = 0  by
 * Newton iteration with a line search on the residual.
 *
 * @param[in]  model Model constants.
 * @param[in]  stress_trial Elastic trial stress.
 * @param[in]  scale Stress level used for the tolerances.
 * @param[out] stress Returned stress.
 * @param[out] ddsdde Consistent tangent.
 * @return 1 on convergence, 0 otherwise.
 */
static int return_to_cone(const MNModel* model, const double stress_trial[VOIGTSIZE_3D],
                          const double scale, double stress[VOIGTSIZE_3D],
                          double ddsdde[VOIGTSIZE_3D_SQ])
{
    double delta_gamma = 0.0;
    double Xi[VOIGTSIZE_3D_SQ];
    MNEvaluation ev;

    copy_array(stress_trial, VOIGTSIZE_3D, stress);
    evaluate(model, stress, stress_trial, delta_gamma, &ev);
    double merit = merit_function(&ev, scale);

    int converged = 0;
    for (int iter = 0; iter < MN_MAX_ITERATIONS; ++iter)
    {
        if (merit <= MN_YIELD_TOL * MN_YIELD_TOL)
        {
            converged = 1;
            break;
        }

        // Newton direction
        if (!calculate_algorithmic_stiffness(model, &ev, delta_gamma, Xi)) return 0;
        double r_eps[VOIGTSIZE_3D], Xi_r[VOIGTSIZE_3D], Xi_g[VOIGTSIZE_3D];
        matrix_vector_multiply(model->Ce_inv, ev.residual, VOIGTSIZE_3D, r_eps);
        matrix_vector_multiply(Xi, r_eps, VOIGTSIZE_3D, Xi_r);
        matrix_vector_multiply(Xi, ev.grad_g, VOIGTSIZE_3D, Xi_g);
        const double denom = vector_dot_product(ev.grad_f, Xi_g, VOIGTSIZE_3D);
        if (!(denom > 0.0)) return 0;

        const double d_gamma = (ev.f - vector_dot_product(ev.grad_f, Xi_r, VOIGTSIZE_3D)) / denom;
        double d_stress[VOIGTSIZE_3D];
        for (int i = 0; i < VOIGTSIZE_3D; ++i) d_stress[i] = -(Xi_r[i] + d_gamma * Xi_g[i]);

        // line search on the merit function; states on (or across) the hydrostatic axis are not
        // admissible, the cone is not differentiable there
        double step = 1.0;
        int accepted = 0;
        for (int ls = 0; ls <= MN_MAX_LINE_SEARCH; ++ls, step *= 0.5)
        {
            double stress_new[VOIGTSIZE_3D];
            for (int i = 0; i < VOIGTSIZE_3D; ++i) stress_new[i] = stress[i] + step * d_stress[i];
            const double delta_gamma_new = delta_gamma + step * d_gamma;

            MNEvaluation ev_new;
            evaluate(model, stress_new, stress_trial, delta_gamma_new, &ev_new);
            const double merit_new = merit_function(&ev_new, scale);
            if (!(ev_new.J > ZERO_TOL * scale) || !isfinite(merit_new)) continue;

            if (merit_new <= (1.0 - 1.0e-4 * step) * merit || ls == MN_MAX_LINE_SEARCH)
            {
                copy_array(stress_new, VOIGTSIZE_3D, stress);
                delta_gamma = delta_gamma_new;
                ev = ev_new;
                merit = merit_new;
                accepted = 1;
                break;
            }
        }
        if (!accepted) return 0;
    }
    if (!converged && !(merit <= MN_YIELD_TOL * MN_YIELD_TOL)) return 0;
    if (delta_gamma < 0.0) return 0;

    // consistent tangent: Xi - (Xi grad_g) (Xi grad_f)^T / (grad_f . Xi grad_g), Xi symmetric
    if (!calculate_algorithmic_stiffness(model, &ev, delta_gamma, Xi)) return 0;
    double Xi_f[VOIGTSIZE_3D], Xi_g[VOIGTSIZE_3D];
    matrix_vector_multiply(Xi, ev.grad_f, VOIGTSIZE_3D, Xi_f);
    matrix_vector_multiply(Xi, ev.grad_g, VOIGTSIZE_3D, Xi_g);
    const double denom = vector_dot_product(ev.grad_f, Xi_g, VOIGTSIZE_3D);
    if (!(denom > 0.0)) return 0;
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
        for (int j = 0; j < VOIGTSIZE_3D; ++j)
            ddsdde[i * VOIGTSIZE_3D + j] = Xi[i * VOIGTSIZE_3D + j] - Xi_g[i] * Xi_f[j] / denom;

    return 1;
}

/**
 * @brief Whether the return to the apex solves the return mapping.
 *
 * The plastic strain of the return to the apex, Ce^-1 (sigma_trial - sigma_apex), must be a plastic
 * flow direction of the potential at its apex, lambda (M_g / 3 I + n) with lambda >= 0 and
 * h*(n) <= 1 (calculate_deviatoric_support_function):
 *     lambda = (p_trial - p_apex) / (K M_g) >= 0   and   h*(s_trial) <= 2 G lambda.
 * For M_g > 0 (dilatancy) this region lies beyond the apex and does not overlap with the trial
 * stresses that can be returned to the cone. For M_g < 0 (contraction) it lies on the compressive
 * side of the apex and overlaps with them.
 */
static int apex_return_is_admissible(const MNModel* model, const double p_trial,
                                     const double j2_trial, const double j3_trial,
                                     const double scale)
{
    const double M_g = model->potential.M;
    if (M_g == 0.0) return 0;

    const double lambda = (p_trial - model->p_apex) / (model->bulk_modulus * M_g);
    if (!(lambda > 0.0)) return 0;
    return calculate_deviatoric_support_function(j2_trial, j3_trial, model->potential) <=
           2.0 * model->shear_modulus * lambda + MN_YIELD_TOL * scale;
}

/**
 * @brief Return to the apex: the stress does not change under further loading, a fraction
 * MN_APEX_STIFFNESS_FRACTION of the elastic stiffness is kept as tangent.
 */
static void return_to_apex(const MNModel* model, double stress[VOIGTSIZE_3D],
                           double ddsdde[VOIGTSIZE_3D_SQ], int* state)
{
    for (int i = 0; i < VOIGTSIZE_3D; ++i) stress[i] = (i < 3) ? model->p_apex : 0.0;
    for (int i = 0; i < VOIGTSIZE_3D_SQ; ++i) ddsdde[i] = MN_APEX_STIFFNESS_FRACTION * model->Ce[i];
    *state = MN_STATE_APEX;
}

/**
 * @brief Elastic predictor and return mapping for one (sub-)step.
 *
 * The return goes to the apex when that is the solution and otherwise to the cone, except:
 *  - beyond the apex a potential with M_g <= 0 cannot bring the stress back to the cone; these
 *    trial stresses are returned to the apex, although the volumetric plastic strain then does not
 *    follow the flow rule;
 *  - with contractant flow (M_g < 0) the cone and the apex can both be solutions; the return to the
 *    cone is preferred, as the solution that is continuous for small increments. For trial stresses
 *    far outside the cone in the deviatoric direction neither is a solution (also not at the apex
 *    itself under continued shearing). With apex_fallback set, these are returned to the apex,
 *    where contractant flow leads anyway, although the plastic strain then does not follow the
 *    flow rule.
 *
 * @return 1 on success, 0 if the return mapping fails.
 */
static int integrate_step(const MNModel* model, const double stress_n[VOIGTSIZE_3D],
                          const double dstrain[VOIGTSIZE_3D], double stress[VOIGTSIZE_3D],
                          double ddsdde[VOIGTSIZE_3D_SQ], int* state, const int apex_fallback)
{
    // elastic predictor
    double stress_trial[VOIGTSIZE_3D], delta_stress[VOIGTSIZE_3D];
    matrix_vector_multiply(model->Ce, dstrain, VOIGTSIZE_3D, delta_stress);
    add_vectors(stress_n, delta_stress, VOIGTSIZE_3D, stress_trial);

    double p_trial, J_trial, theta_trial, j2_trial, j3_trial, s_dev[VOIGTSIZE_3D];
    calculate_stress_invariants_3d(stress_trial, &p_trial, &J_trial, &theta_trial, &j2_trial,
                                   &j3_trial, s_dev);
    double f_trial;
    calculate_yield_function(p_trial, theta_trial, J_trial, model->yield, &f_trial);

    double scale = fabs(p_trial) + J_trial + fabs(model->yield.K);
    if (!(scale > 0.0)) scale = 1.0;

    if (f_trial <= MN_YIELD_TOL * scale)
    {
        copy_array(stress_trial, VOIGTSIZE_3D, stress);
        copy_array(model->Ce, VOIGTSIZE_3D_SQ, ddsdde);
        *state = MN_STATE_ELASTIC;
        return 1;
    }

    const double M_g = model->potential.M;
    if (model->has_apex)
    {
        const int beyond_apex = p_trial - model->p_apex >= -MN_YIELD_TOL * scale;
        if ((M_g <= 0.0 && beyond_apex) ||
            (M_g > 0.0 && apex_return_is_admissible(model, p_trial, j2_trial, j3_trial, scale)))
        {
            return_to_apex(model, stress, ddsdde, state);
            return 1;
        }
    }

    if (return_to_cone(model, stress_trial, scale, stress, ddsdde))
    {
        *state = MN_STATE_PLASTIC;
        return 1;
    }

    if (model->has_apex && M_g < 0.0 &&
        (apex_fallback || apex_return_is_admissible(model, p_trial, j2_trial, j3_trial, scale)))
    {
        return_to_apex(model, stress, ddsdde, state);
        return 1;
    }
    return 0;
}

/**
 * @brief Integrates the strain increment, with sub-steps if the return mapping of the full
 * increment does not converge. DDSDDE is the tangent of the last (sub-)step. The apex fallback of
 * contractant flow (see integrate_step) is only used with the smallest sub-steps.
 *
 * @return 1 on success, 0 otherwise.
 */
static int integrate(const MNModel* model, const double stress_n[VOIGTSIZE_3D],
                     const double dstrain[VOIGTSIZE_3D], double stress[VOIGTSIZE_3D],
                     double ddsdde[VOIGTSIZE_3D_SQ], int* state)
{
    if (integrate_step(model, stress_n, dstrain, stress, ddsdde, state, 0)) return 1;

    for (int level = 1; level <= MN_MAX_SUBSTEP_LEVEL; ++level)
    {
        const int n_substeps = 1 << level;
        double dstrain_sub[VOIGTSIZE_3D], stress_sub[VOIGTSIZE_3D];
        vector_scalar_multiply(dstrain, 1.0 / n_substeps, VOIGTSIZE_3D, dstrain_sub);
        copy_array(stress_n, VOIGTSIZE_3D, stress_sub);

        const int apex_fallback = level == MN_MAX_SUBSTEP_LEVEL;
        int ok = 1;
        for (int k = 0; k < n_substeps && ok; ++k)
        {
            ok = integrate_step(model, stress_sub, dstrain_sub, stress, ddsdde, state, apex_fallback);
            copy_array(stress, VOIGTSIZE_3D, stress_sub);
        }
        if (ok) return 1;
    }
    return 0;
}

// Define the UMAT function signature expected by the FEA software.
//       Check your specific FEA software documentation for exact C interface requirements if
//       available. Some systems might require all arguments to be pointers, even scalars.

UMAT_EXPORT void UMAT_CALLCONV umat(
    // Outputs (to be updated by the subroutine)
    double* STRESS,  // Stress tensor at end of increment (NTENS components)
    double* STATEV,  // State variables at end of increment (NSTATV components)
    double* DDSDDE,  // Jacobian matrix (NTENS * NTENS components)
    double* SSE,     // Specific elastic strain energy
    double* SPD,     // Plastic dissipation
    double* SCD,     // Creep dissipation
    double* RPL,     // Volumetric heat generation
    double* DDSDDT,  // Stress rate dependency on temperature (NTENS components)
    double* DRPLDE,  // Derivative of RPL wrt strain (NTENS components)
    double* DRPLDT,  // Derivative of RPL wrt temperature
    // Inputs (provided by the FEA software)
    double* STRAN,   // Total strain at start of increment (NTENS components)
    double* DSTRAN,  // Increment in total strain (NTENS components)
    double* TIME,    // Step time [0] and total time [1]
    double* DTIME,   // Time increment
    double* TEMP,    // Temperature at start of increment
    double* DTEMP,   // Increment in temperature
    double* PREDEF,  // Predefined field variables at start (NPREDFIELD components)
    double* DPRED,   // Increment in predefined field variables (NPREDFIELD components)
    char* CMNAME,    // Material name (passed typically as CHARACTER*80 in Fortran)
    int* NDI,        // Number of direct stress components (e.g., 3 for 3D)
    int* NSHR,       // Number of shear stress components (e.g., 3 for 3D)
    int* NTENS,      // Total number of stress components (NDI + NSHR)
    int* NSTATV,     // Number of state variables
    double* PROPS,   // User-defined material properties (NPROPS components)
    int* NPROPS,     // Number of properties
    double* COORDS,  // Coordinates of the integration point (3 components)
    double* DROT,    // Rotation increment matrix (3x3 = 9 components)
    double* PNEWDT,  // Suggested new time increment size (can be modified)
    double* CELENT,  // Characteristic element length
    double* DFGRD0,  // Deformation gradient at start (9 components)
    double* DFGRD1,  // Deformation gradient at end (9 components)
    int* NOEL,       // Element number
    int* NPT,        // Integration point number
    int* LAYER,      // Layer number (for shells/beams)
    int* KSPT,       // Section point number
    int* KSTEP,      // Step number
    int* KINC        // Increment number
    // Note: Size of CMNAME requires careful handling between C and Fortran
)
{
    // avoid unused variable warnings
    (void)KINC;
    (void)KSTEP;
    (void)KSPT;
    (void)LAYER;
    (void)NPT;
    (void)NOEL;
    (void)DFGRD0;
    (void)DFGRD1;
    (void)CELENT;
    (void)DROT;
    (void)COORDS;
    (void)CMNAME;
    (void)DPRED;
    (void)PREDEF;
    (void)DTEMP;
    (void)TEMP;
    (void)DTIME;
    (void)TIME;
    (void)STRAN;
    (void)DRPLDT;
    (void)DRPLDE;
    (void)DDSDDT;
    (void)RPL;
    (void)SCD;

    // Check Inputs
    if (*NTENS != VOIGTSIZE_3D || *NDI != 3 || *NSHR != 3)
    {
        // Handle error - This UMAT is specifically for 3D
        // For simplicity, we'll print an error and potentially stop (though stopping is usually
        // bad)
        fprintf(stderr, "UMAT Error: NTENS != 6. This UMAT requires 3D elements.\n");
        // exit(1); // Avoid exiting in production code if possible
        return;  // Or try to handle gracefully
    }
    if (check_properties(*NPROPS, PROPS))
    {
        fprintf(stderr, "UMAT Error: invalid material properties.\n");
        return;
    }

    if (*NSTATV < 1)
    {
        fprintf(stderr, "UMAT Error: NSTATV < 1. Requires at least 1 state variable.\n");
        return;
    }

    // Material Properties
    const double E_mod = PROPS[0];    // Young's Modulus
    const double nu = PROPS[1];       // Poisson's Ratio
    const double c = PROPS[2];        // Cohesion
    const double phi_deg = PROPS[3];  // Friction angle
    const double psi_deg = PROPS[4];  // Dilation angle

    // Convert angles to radians
    const double phi_rad = phi_deg * PI / 180.0;
    const double psi_rad = psi_deg * PI / 180.0;

    MNModel model;
    calculate_elastic_stiffness_matrix_3d(E_mod, nu, model.Ce);
    calculate_elastic_compliance_matrix_3d(E_mod, nu, model.Ce_inv);
    model.bulk_modulus = E_mod / (3.0 * (1.0 - 2.0 * nu));
    model.shear_modulus = E_mod / (2.0 * (1.0 + nu));
    model.yield = calculate_matsuoka_nakai_constants(phi_rad, c);
    model.potential = calculate_matsuoka_nakai_constants(psi_rad, c);
    model.has_apex = model.yield.M > 0.0;
    model.p_apex = model.has_apex ? model.yield.K / model.yield.M : 0.0;

    double stress[VOIGTSIZE_3D];
    int state = MN_STATE_ELASTIC;
    if (!integrate(&model, STRESS, DSTRAN, stress, DDSDDE, &state))
    {
        fprintf(stderr, "UMAT Warning: Matsuoka-Nakai return mapping did not converge; "
                        "requesting a smaller time increment.\n");
        // keep the stress at the start of the increment, with the elastic tangent
        copy_array(model.Ce, VOIGTSIZE_3D_SQ, DDSDDE);
        if (PNEWDT) *PNEWDT = 0.5;
        return;
    }

    // elastic strain increment Ce^-1 (stress - STRESS) and plastic strain increment
    double delta_stress[VOIGTSIZE_3D], dEps_el[VOIGTSIZE_3D], dEps_p[VOIGTSIZE_3D];
    double mean_stress[VOIGTSIZE_3D];
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        delta_stress[i] = stress[i] - STRESS[i];
        mean_stress[i] = 0.5 * (stress[i] + STRESS[i]);
    }
    matrix_vector_multiply(model.Ce_inv, delta_stress, VOIGTSIZE_3D, dEps_el);
    for (int i = 0; i < VOIGTSIZE_3D; ++i) dEps_p[i] = DSTRAN[i] - dEps_el[i];

    // elastic strain energy (exact for linear elasticity) and plastic dissipation
    if (SSE) *SSE += vector_dot_product(mean_stress, dEps_el, VOIGTSIZE_3D);
    if (SPD) *SPD += vector_dot_product(stress, dEps_p, VOIGTSIZE_3D);

    copy_array(stress, VOIGTSIZE_3D, STRESS);
    STATEV[0] = (double)state;

    return;
}

int check_properties(const int NPROPS, const double* PROPS)
{
    if (NPROPS < 5)
    {
        fprintf(stderr, "UMAT Error: NPROPS < 5. Requires E, nu,c, phi and psi.\n");
        return 1;
    }

    int n_errors = 0;

    if (PROPS[0] <= 0.0)
    {
        fprintf(stderr, "UMAT Error: Young's Modulus must be positive.\n");
        n_errors++;
    }
    if (PROPS[1] < 0.0 || PROPS[1] >= 0.5)
    {
        fprintf(stderr, "UMAT Error: Poisson's Ratio must be between [0.0, 0.5).\n");
        n_errors++;
    }
    if (PROPS[2] < 0.0)
    {
        fprintf(stderr, "UMAT Error: Cohesion must be non-negative.\n");
        n_errors++;
    }
    if (PROPS[3] < 0.0 || PROPS[3] >= 90.0)
    {
        fprintf(stderr, "UMAT Error: Friction angle must be between [0.0, 90.0).\n");
        n_errors++;
    }
    if (PROPS[4] <= -90.0 || PROPS[4] >= 90.0)
    {
        fprintf(stderr, "UMAT Error: Dilation angle must be between (-90.0, 90.0).\n");
        n_errors++;
    }

    if (n_errors > 0)
    {
        return 1;
    }
    return 0;  // All checks passed
}
