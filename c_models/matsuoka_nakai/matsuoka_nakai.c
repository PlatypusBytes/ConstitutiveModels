#include <stdio.h>
#include <math.h>

#include "../utils.h"
#include "../yield_surfaces/matsuoka_nakai_surface.h"
#include "../stress_utils.h"
#include "../globals.h"

//#include <stdlib.h> // For exit() if needed, though usually avoided in UMAT

// Define necessary calling conventions and export macros (adjust for your compiler/system)
// For MSVC on Windows:
#if defined(_WIN32) || defined(_WIN64)
#define UMAT_EXPORT __declspec(dllexport)
#define UMAT_CALLCONV __stdcall // Abaqus often uses stdcall
#else
// For GCC/Clang on Linux/macOS (usually no special decoration needed)
#define UMAT_EXPORT
#define UMAT_CALLCONV
#endif


void calculate_elastic_stiffness(double E, double nu, double* DDSDDE, int NTENS);


// Define the UMAT function signature expected by the FEA software.
//       Check your specific FEA software documentation for exact C interface requirements if available.
//       Some systems might require all arguments to be pointers, even scalars.

UMAT_EXPORT void UMAT_CALLCONV umat(
    // Outputs (to be updated by the subroutine)
    double* STRESS,      // Stress tensor at end of increment (NTENS components)
    double* STATEV,      // State variables at end of increment (NSTATV components)
    double* DDSDDE,      // Jacobian matrix (NTENS * NTENS components)
    double* SSE,         // Elastic strain energy density
    double* SPD,         // Plastic dissipation
    double* SCD,         // Creep dissipation
    double* RPL,         // Volumetric heat generation
    double* DDSDDT,      // Stress rate dependency on temperature (NTENS components)
    double* DRPLDE,      // Derivative of RPL wrt strain (NTENS components)
    double* DRPLDT,      // Derivative of RPL wrt temperature
    // Inputs (provided by the FEA software)
    double* STRAN,       // Total strain at start of increment (NTENS components)
    double* DSTRAN,      // Increment in total strain (NTENS components)
    double* TIME,        // Step time [0] and total time [1]
    double* DTIME,       // Time increment
    double* TEMP,        // Temperature at start of increment
    double* DTEMP,       // Increment in temperature
    double* PREDEF,      // Predefined field variables at start (NPREDFIELD components)
    double* DPRED,       // Increment in predefined field variables (NPREDFIELD components)
    char* CMNAME,      // Material name (passed typically as CHARACTER*80 in Fortran)
    int* NDI,         // Number of direct stress components (e.g., 3 for 3D)
    int* NSHR,        // Number of shear stress components (e.g., 3 for 3D)
    int* NTENS,       // Total number of stress components (NDI + NSHR)
    int* NSTATV,      // Number of state variables
    double* PROPS,       // User-defined material properties (NPROPS components)
    int* NPROPS,      // Number of properties
    double* COORDS,      // Coordinates of the integration point (3 components)
    double* DROT,        // Rotation increment matrix (3x3 = 9 components)
    double* PNEWDT,      // Suggested new time increment size (can be modified)
    double* CELENT,      // Characteristic element length
    double* DFGRD0,      // Deformation gradient at start (9 components)
    double* DFGRD1,      // Deformation gradient at end (9 components)
    int* NOEL,        // Element number
    int* NPT,         // Integration point number
    int* LAYER,       // Layer number (for shells/beams)
    int* KSPT,        // Section point number
    int* KSTEP,       // Step number
    int* KINC         // Increment number
    // Note: Size of CMNAME requires careful handling between C and Fortran
) {

    // --- 0. Initialization and Material Properties ---
    int i, j;
    int n_dim_dir = *NDI;
    int n_dim_shr = *NSHR;
    int n_tensor = *NTENS; // Should be 6 for 3D

    // --- 0. Check Inputs ---
    if (*NTENS != 6 || *NDI != 3 || *NSHR != 3) {
        // Handle error - This UMAT is specifically for 3D
        // For simplicity, we'll print an error and potentially stop (though stopping is usually bad)
        fprintf(stderr, "UMAT Error: NTENS != 6. This UMAT requires 3D elements.\n");
        // exit(1); // Avoid exiting in production code if possible
        return; // Or try to handle gracefully
    }
    if (*NPROPS < 3) {
        fprintf(stderr, "UMAT Error: NPROPS < 3. Requires E, nu, tension_threshold.\n");
        return;
    }
    if (*NSTATV < 1) {
        fprintf(stderr, "UMAT Error: NSTATV < 1. Requires at least 1 state variable.\n");
        return;
    }


    // --- 1. Material Properties ---
    double E_mod = PROPS[0];                // Young's Modulus
    double nu = PROPS[1];               // Poisson's Ratio
    double c = PROPS[2];                // Cohesion (not used in this example)
    double phi_deg = PROPS[3];              // Friction angle (not used in this example)
    double psi_deg = PROPS[4];               // Dilation angle (not used in this example)

    // Convert angles to radians
    double phi_rad = phi_deg * PI / 180.0;
    double psi_rad = psi_deg * PI / 180.0;
    double sin_phi = sin(phi_rad);
    double cos_phi = cos(phi_rad);
    double sin_psi = sin(psi_rad);

    // State variables
    double peeq_n = STATEV[0]; // Equivalent plastic strain at start of increment

        // Local arrays
    double stress_trial[6];
    double Ce_matrix[36]; // Elastic stiffness (6x6 as 1D row-major)
    double s_dev[6];      // Deviatoric stress tensor (Voigt)
    double grad_f[6];     // Gradient of yield function df/dsigma (A_vec)
    double grad_g[6];     // Gradient of plastic potential dg/dsigma (g_vec, flow vector)
    double Ce_grad_g[6];  // Ce * grad_g
    double Ce_grad_f[6];  // Ce * grad_f
    double dEps_p[6];     // Plastic strain increment
    double dp_dsig[6];
    double dJ_dsig[6];
    double dtheta_dsig[6];

    // --- 1. Calculate Elastic Stiffness Matrix ---
    calculate_elastic_stiffness(E_mod, nu, Ce_matrix, n_tensor);

    // Initialize Jacobian DDSDDE to elastic matrix (default assumption)
    for (i = 0; i < n_tensor * n_tensor; ++i) {
        DDSDDE[i] = Ce_matrix[i];
    }

    // --- 2. Elastic Predictor Step ---
    // stress_trial = STRESS_n + Ce * DSTRAN
    matrix_vector_multiply(Ce_matrix, DSTRAN, n_tensor, Ce_grad_g); // Use Ce_grad_g as temporary storage for Ce*DSTRAN
    for (i = 0; i < n_tensor; ++i) {
        stress_trial[i] = STRESS[i] + Ce_grad_g[i];
    }

    // calculate invariants
    double p_trial, J_trial, theta_trial, j2_trial, j3_trial;
    calculate_stress_invariants_3d(stress_trial,&p_trial, &J_trial, &theta_trial, &j2_trial, &j3_trial, s_dev); // s_dev also calculated here

    // matsuoka nakai constants
    double alpha =0;
    double beta =0;
    double gamma =0;
    double K =0;
    double M =0;

    // yield function
    double f_trial = 0.0;

    // Calculate yield function value
    calculate_matsuoka_nakai_constants(phi_rad, c, &alpha, &beta, &gamma, &K, &M);
    calculate_yield_function(p_trial, theta_trial,J_trial, c, phi_rad,alpha,beta,gamma,K,M, &f_trial);

    // yield function greater than zero, calculate plastic correction
    if (f_trial > ZERO_TOL)
    {

        // gradient yield function
        calculate_stress_invariants_derivatives_3d(J_trial, s_dev, j2_trial, j3_trial,dp_dsig,  dJ_dsig,  dtheta_dsig);


        double mats_nak_constants[5] = {phi_rad, M, alpha, beta, gamma};

        calculate_yield_gradient(theta_trial, J_trial, mats_nak_constants, dp_dsig, dJ_dsig, dtheta_dsig, grad_f);

        // gradient potential function, g, it is required to recalculate the constants using psi
        calculate_matsuoka_nakai_constants(psi_rad, c, &alpha, &beta, &gamma, &K, &M);
        double mats_nak_constants_plastic_potential[5] ={psi_rad, M, alpha, beta, gamma};

        calculate_yield_gradient(theta_trial, J_trial, mats_nak_constants_plastic_potential, dp_dsig, dJ_dsig, dtheta_dsig, grad_g);


        // Calculate terms needed for delta_gamma and Jacobian
        matrix_vector_multiply(Ce_matrix, grad_g, n_tensor, Ce_grad_g); // Ce * g
        matrix_vector_multiply(Ce_matrix, grad_f, n_tensor, Ce_grad_f); // Ce * f


        // Calculate denominator for delta_gamma (and Jacobian)
        // Assumes perfect plasticity (Hardening modulus H=0)
        // Denom = A_vec : Ce : g_vec = grad_f . (Ce * grad_g)
        double denom = vector_dot_product(grad_f, Ce_grad_g, n_tensor);
        double delta_gamma = f_trial / denom;

        // --- Update Stress ---
        // STRESS_{n+1} = stress_trial - delta_gamma * Ce * g_vec
        for (i = 0; i < n_tensor; ++i) {
            STRESS[i] = stress_trial[i] - delta_gamma * Ce_grad_g[i];
        }

        // --- Update State Variables (PEEQ) ---
        // dEps_p = delta_gamma * g_vec
        for(i=0; i<n_tensor; ++i) {
            dEps_p[i] = delta_gamma * grad_g[i];
        }
        double dSpd = vector_dot_product(STRESS, dEps_p, n_tensor);
        *SPD += dSpd;

        // --- Calculate Consistent Tangent Modulus (Jacobian DDSDDE) ---
        // DDSDDE = Ce - (Ce * g) cross_prod (Ce * f)^T / denom
        if (fabs(denom) > ZERO_TOL) { // Redundant check, but safe
            for (i = 0; i < n_tensor; ++i) {
                for (j = 0; j < n_tensor; ++j) {
                    // Accessing 1D array DDSDDE with row-major logic
                    DDSDDE[i * n_tensor + j] = Ce_matrix[i * n_tensor + j] - (Ce_grad_g[i] * Ce_grad_f[j]) / denom ;
                }
            }
        }
    }
    else{

        // --- Elastic Step ---
        // No yield, treat as elastic
        // STRESS_{n+1} = stress_trial
        // DDSDDE remains Ce_matrix (already set)
        // STATEV[0] remains peeq_n


        for (i = 0; i < n_tensor; ++i) {
            STRESS[i] = stress_trial[i]; // Use trial stress
        }
        // DDSDDE remains Ce_matrix
        // STATEV[0] remains peeq_n
        }

     *SCD = 0.0; // No creep
    return;

}

void calculate_elastic_stiffness(double E, double nu, double* DDSDDE, int NTENS) {
    if (NTENS != 6) return; // Only for 3D

    double G = E / (2.0 * (1.0 + nu));      // Shear modulus
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); // Lame's first param

    double factor = lambda + 2.0 * G;

    // Initialize to zero
    for (int i = 0; i < NTENS * NTENS; ++i) {
        DDSDDE[i] = 0.0;
    }

    // Populate using row-major indexing: DDSDDE[row * NTENS + col]
    // Normal Stresses
    DDSDDE[0 * NTENS + 0] = factor; // C_1111 (row 0, col 0)
    DDSDDE[0 * NTENS + 1] = lambda; // C_1122 (row 0, col 1)
    DDSDDE[0 * NTENS + 2] = lambda; // C_1133 (row 0, col 2)

    DDSDDE[1 * NTENS + 0] = lambda; // C_2211 (row 1, col 0)
    DDSDDE[1 * NTENS + 1] = factor; // C_2222 (row 1, col 1)
    DDSDDE[1 * NTENS + 2] = lambda; // C_2233 (row 1, col 2)

    DDSDDE[2 * NTENS + 0] = lambda; // C_3311 (row 2, col 0)
    DDSDDE[2 * NTENS + 1] = lambda; // C_3322 (row 2, col 1)
    DDSDDE[2 * NTENS + 2] = factor; // C_3333 (row 2, col 2)

    // Shear Stresses (using engineering shear strain convention gamma = 2*epsilon_shear)
    // The stiffness term C_ij = G for i != j in engineering strain notation
    DDSDDE[3 * NTENS + 3] = G;      // C_1212 (row 3, col 3)
    DDSDDE[4 * NTENS + 4] = G;      // C_2323 (row 4, col 4)
    DDSDDE[5 * NTENS + 5] = G;      // C_1313 (row 5, col 5)
}
