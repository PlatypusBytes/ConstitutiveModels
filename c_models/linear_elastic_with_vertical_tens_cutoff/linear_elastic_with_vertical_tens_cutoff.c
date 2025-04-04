#include <stdio.h>
#include <math.h>
#include <stdlib.h> // For exit() if needed, though usually avoided in UMAT

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

// Define the UMAT function signature expected by the FEA software.
// NOTE: This is a C adaptation. Abaqus *officially* expects Fortran.
//       Argument names match Fortran standard for clarity.
//       Pointers are used for all arrays/outputs.
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
    // --- 0. Check Inputs ---
    if (*NTENS != 6 || *NDI != 3 || *NSHR != 3) {
        // Handle error - This UMAT is specifically for 3D
        // In Abaqus, you might write to MSG file or call utility routines.
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
    double E = PROPS[0];                // Young's Modulus
    double nu = PROPS[1];               // Poisson's Ratio
    double tension_threshold = PROPS[2]; // Tensile strength threshold in Y

    // --- 2. Calculate Elastic Constants ---
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); // Lame's first parameter
    double G = E / (2.0 * (1.0 + nu));                     // Shear modulus (Lame's second parameter)
    double K = E / (3.0 * (1.0 - 2.0 * nu));                 // Bulk modulus (optional, for checks)

    // --- 3. Calculate Elastic Stiffness Matrix (DDSDDE_elastic) ---
    // Stored as a 1D array (NTENS * NTENS). Assuming column-major storage like Fortran.
    // C = | C11 C12 C13  0   0   0  |   Indices (0-based C array):
    //     | C12 C22 C23  0   0   0  |   Row 0: 0, 1, 2, 3, 4, 5
    //     | C13 C23 C33  0   0   0  |   Row 1: 6, 7, 8, 9, 10, 11
    //     |  0   0   0  C44  0   0  |   Row 2: 12, 13, 14, 15, 16, 17
    //     |  0   0   0   0  C55  0  |   ... etc ...
    //     |  0   0   0   0   0  C66 |
    // C11 = C22 = C33 = lambda + 2G
    // C12 = C13 = C23 = lambda
    // C44 = C55 = C66 = G

    double DDSDDE_elastic[36]; // NTENS * NTENS = 6 * 6 = 36
    for (int i = 0; i < 36; ++i) DDSDDE_elastic[i] = 0.0; // Initialize

    double C11 = lambda + 2.0 * G;
    double C12 = lambda;
    double C44 = G;

    // Diagonal terms
    DDSDDE_elastic[0] = C11; // C(1,1) index 0
    DDSDDE_elastic[7] = C11; // C(2,2) index 7 = 1*6 + 1
    DDSDDE_elastic[14] = C11; // C(3,3) index 14 = 2*6 + 2
    DDSDDE_elastic[21] = C44; // C(4,4) index 21 = 3*6 + 3
    DDSDDE_elastic[28] = C44; // C(5,5) index 28 = 4*6 + 4
    DDSDDE_elastic[35] = C44; // C(6,6) index 35 = 5*6 + 5

    // Off-diagonal terms (symmetric)
    DDSDDE_elastic[1] = C12; // C(1,2) index 1 = 0*6 + 1
    DDSDDE_elastic[6] = C12; // C(2,1) index 6 = 1*6 + 0

    DDSDDE_elastic[2] = C12; // C(1,3) index 2 = 0*6 + 2
    DDSDDE_elastic[12] = C12; // C(3,1) index 12 = 2*6 + 0

    DDSDDE_elastic[8] = C12; // C(2,3) index 8 = 1*6 + 2
    DDSDDE_elastic[13] = C12; // C(3,2) index 13 = 2*6 + 1


    // --- 4. Calculate Strain at End of Increment ---
    double strain_end[6];
    for (int i = 0; i < 6; ++i) {
        strain_end[i] = STRAN[i] + DSTRAN[i];
    }

    // --- 5. Calculate Elastic Trial Stress ---
    // stress_trial = DDSDDE_elastic * strain_end
    // Matrix-vector multiplication (remembering column-major storage for DDSDDE)
    double stress_trial[6];
    for (int i = 0; i < 6; ++i) {
        stress_trial[i] = 0.0;
        for (int j = 0; j < 6; ++j) {
            // DDSDDE_elastic[row + col*NROWS] = DDSDDE_elastic[i + j*6]
            stress_trial[i] += DDSDDE_elastic[i + j * 6] * strain_end[j];
        }
    }

    // --- 6. Apply Tension Cutoff Logic ---
    // Check vertical stress component (sigma_yy, index 1)
    double sigma_yy_trial = stress_trial[1];
    int tension_cutoff_active = 0; // Flag: 0 = elastic, 1 = cutoff

    if (sigma_yy_trial > tension_threshold) {
        tension_cutoff_active = 1;
    }

    // --- 7. Update Stress, Jacobian, and State Variables ---
    if (tension_cutoff_active) {
        // Tension cutoff is active: Zero out all stresses and the Jacobian
        for (int i = 0; i < 6; ++i) {
            STRESS[i] = 0.0;
        }
        for (int i = 0; i < 36; ++i) {
            // Set Jacobian to zero (or a very small number for numerical stability if needed)
            DDSDDE[i] = 0.0; // Or maybe 1.0e-12 * DDSDDE_elastic[i] ? Test this.
        }
        STATEV[0] = 1.0; // Indicate cutoff state
    }
    else {
        // Elastic behavior: Use the trial stress and elastic Jacobian
        for (int i = 0; i < 6; ++i) {
            STRESS[i] = stress_trial[i];
        }
        for (int i = 0; i < 36; ++i) {
            DDSDDE[i] = DDSDDE_elastic[i];
        }
        STATEV[0] = 0.0; // Indicate elastic state
    }

    //// --- 8. Update Optional Outputs (Set to Zero for Simplicity) ---
    //*SSE = 0.0; // Could calculate 0.5 * stress * elastic_strain if needed
    //*SPD = 0.0;
    //*SCD = 0.0;
    //*RPL = 0.0;
    //// DDSDDT, DRPLDE, DRPLDT are often zero unless thermal coupling is active
    //for (int i = 0; i < *NTENS; ++i) DDSDDT[i] = 0.0;
    //for (int i = 0; i < *NTENS; ++i) DRPLDE[i] = 0.0;
    //*DRPLDT = 0.0;

    // PNEWDT can be reduced if convergence is difficult, e.g., *PNEWDT = 0.5 * (*PNEWDT);
    // No change here means *PNEWDT remains as suggested by Abaqus.

    // --- 9. Return ---
    return;
}