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
#define UMAT_CALLCONV   // Abaqus often uses stdcall
#else
// For GCC/Clang on Linux/macOS (usually no special decoration needed)
#define UMAT_EXPORT
#define UMAT_CALLCONV
#endif

// Define the UMAT function signature expected by the FEA software.
//       Check your specific FEA software documentation for exact C interface requirements if
//       available. Some systems might require all arguments to be pointers, even scalars.

UMAT_EXPORT void UMAT_CALLCONV umat(
    // Outputs (to be updated by the subroutine)
    double* STRESS,  // Stress tensor at end of increment (NTENS components)
    double* STATEV,  // State variables at end of increment (NSTATV components)
    double* DDSDDE,  // Jacobian matrix (NTENS * NTENS components)
    double* SSE,     // Elastic strain energy density
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
    (void)PNEWDT;
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
    (void)SSE;

    // --- 0. Initialization and Material Properties ---
    int i;

    // --- 0. Check Inputs ---
    if (*NTENS != VOIGTSIZE_3D || *NDI != 3 || *NSHR != 3)
    {
        // Handle error - This UMAT is specifically for 3D
        // For simplicity, we'll print an error and potentially stop (though stopping is usually
        // bad)
        fprintf(stderr, "UMAT Error: NTENS != 6. This UMAT requires 3D elements.\n");
        // exit(1); // Avoid exiting in production code if possible
        return;  // Or try to handle gracefully
    }
    if (*NPROPS < 5)
    {
        fprintf(stderr, "UMAT Error: NPROPS < 5. Requires E, nu,c, phi and psi.\n");
        return;
    }
    if (*NSTATV < 1)
    {
        fprintf(stderr, "UMAT Error: NSTATV < 1. Requires at least 1 state variable.\n");
        return;
    }

    // --- 1. Material Properties ---
    double E_mod = PROPS[0];    // Young's Modulus
    double nu = PROPS[1];       // Poisson's Ratio
    double c = PROPS[2];        // Cohesion (not used in this example)
    double phi_deg = PROPS[3];  // Friction angle (not used in this example)
    double psi_deg = PROPS[4];  // Dilation angle (not used in this example)

    // Convert angles to radians
    double phi_rad = phi_deg * PI / 180.0;
    double psi_rad = psi_deg * PI / 180.0;

    // Local arrays
    double stress_trial[VOIGTSIZE_3D];
    double delta_stress[VOIGTSIZE_3D];
    double Ce_matrix[VOIGTSIZE_3D * VOIGTSIZE_3D];  // Elastic stiffness (6x6 as 1D row-major)
    double s_dev[VOIGTSIZE_3D];                     // Deviatoric stress tensor (Voigt)
    double grad_f[VOIGTSIZE_3D];                    // Gradient of yield function df/dsigma (A_vec)
    double grad_g[VOIGTSIZE_3D];     // Gradient of plastic potential dg/dsigma (g_vec, flow vector)
    double Ce_grad_g[VOIGTSIZE_3D];  // Ce * grad_g
    double Ce_grad_f[VOIGTSIZE_3D];  // Ce * grad_f
    double dEps_p[VOIGTSIZE_3D];     // Plastic strain increment
    double dp_dsig[VOIGTSIZE_3D];
    double dJ_dsig[VOIGTSIZE_3D];
    double dtheta_dsig[VOIGTSIZE_3D];

    // --- 1. Calculate Elastic Stiffness Matrix ---
    calculate_elastic_stiffness_matrix_3d(E_mod, nu, Ce_matrix);

    // Initialize Jacobian DDSDDE to elastic matrix (default assumption)
    copy_array(Ce_matrix, VOIGTSIZE_3D * VOIGTSIZE_3D, DDSDDE);

    // --- 2. Elastic Predictor Step ---
    // stress_trial = STRESS_n + Ce * DSTRAN
    matrix_vector_multiply(Ce_matrix, DSTRAN, VOIGTSIZE_3D, delta_stress);

    add_vectors(STRESS, delta_stress, VOIGTSIZE_3D,
                stress_trial);  // stress_trial = STRESS + delta_stress

    // calculate invariants
    double p_trial, J_trial, theta_trial, j2_trial, j3_trial;
    calculate_stress_invariants_3d(stress_trial, &p_trial, &J_trial, &theta_trial, &j2_trial,
                                   &j3_trial, s_dev);  // s_dev also calculated here

    // matsuoka nakai constants
    double alpha = 0;
    double beta = 0;
    double gamma = 0;
    double K = 0;
    double M = 0;

    // yield function
    double f_trial = 0.0;

    // Calculate yield function value
    calculate_matsuoka_nakai_constants(phi_rad, c, &alpha, &beta, &gamma, &K, &M);
    calculate_yield_function(p_trial, theta_trial, J_trial, alpha, beta, gamma, K, M, &f_trial);

    // yield function greater than zero, calculate plastic correction
    if (f_trial > ZERO_TOL)
    {
        // gradient yield function
        calculate_stress_invariants_derivatives_3d(J_trial, s_dev, j2_trial, j3_trial, dp_dsig,
                                                   dJ_dsig, dtheta_dsig);

        double mats_nak_constants[5] = {phi_rad, M, alpha, beta, gamma};

        calculate_yield_gradient(theta_trial, J_trial, mats_nak_constants, dp_dsig, dJ_dsig,
                                 dtheta_dsig, grad_f);

        // gradient potential function, g, it is required to recalculate the constants using psi
        calculate_matsuoka_nakai_constants(psi_rad, c, &alpha, &beta, &gamma, &K, &M);
        double mats_nak_constants_plastic_potential[5] = {psi_rad, M, alpha, beta, gamma};

        calculate_yield_gradient(theta_trial, J_trial, mats_nak_constants_plastic_potential,
                                 dp_dsig, dJ_dsig, dtheta_dsig, grad_g);

        // Calculate terms needed for delta_gamma and Jacobian
        matrix_vector_multiply(Ce_matrix, grad_g, VOIGTSIZE_3D, Ce_grad_g);  // Ce * g
        matrix_vector_multiply(Ce_matrix, grad_f, VOIGTSIZE_3D, Ce_grad_f);  // Ce * f

        // Calculate denominator for delta_gamma (and Jacobian)
        // Assumes perfect plasticity (Hardening modulus H=0)
        // Denom = A_vec : Ce : g_vec = grad_f . (Ce * grad_g)
        double denom = vector_dot_product(grad_f, Ce_grad_g, VOIGTSIZE_3D);
        double delta_gamma = f_trial / denom;

        // --- Update Stress ---
        // STRESS_{n+1} = stress_trial - delta_gamma * Ce * g_vec
        for (i = 0; i < VOIGTSIZE_3D; ++i)
        {
            STRESS[i] = stress_trial[i] - delta_gamma * Ce_grad_g[i];
        }

        // dEps_p = delta_gamma * g_vec
        vector_scalar_multiply(grad_g, delta_gamma, VOIGTSIZE_3D, dEps_p);

        double dSpd = vector_dot_product(STRESS, dEps_p, VOIGTSIZE_3D);
        *SPD += dSpd;

        // Calculate outer product of Ce_grad_g and Ce_grad_f
        double Ce_grad_g_outer_Ce_grad_f[VOIGTSIZE_3D * VOIGTSIZE_3D];
        vector_outer_product(Ce_grad_g, Ce_grad_f, VOIGTSIZE_3D, Ce_grad_g_outer_Ce_grad_f);

        // --- Calculate Consistent Tangent Modulus (Jacobian DDSDDE) ---
        // DDSDDE = Ce - (Ce * g) outer_prod (Ce * f)^T / denom
        if (fabs(denom) > ZERO_TOL)
        {  // Redundant check, but safe
            for (i = 0; i < VOIGTSIZE_3D * VOIGTSIZE_3D; ++i)
            {
                // Accessing 1D array DDSDDE with row-major logic
                DDSDDE[i] = Ce_matrix[i] - Ce_grad_g_outer_Ce_grad_f[i] / denom;
            }
        }

        STATEV[0] = 1.0;  // Indicate that yield surface was reached
    }
    else
    {
        // --- Elastic Step ---
        // No yield, treat as elastic
        // STRESS_{n+1} = stress_trial
        // DDSDDE remains Ce_matrix (already set)

        copy_array(stress_trial, VOIGTSIZE_3D, STRESS);  // Copy trial stress to output stress

        STATEV[0] = 0.0;  // No plastic strain increment
    }

    *SCD = 0.0;  // No creep
    return;
}
