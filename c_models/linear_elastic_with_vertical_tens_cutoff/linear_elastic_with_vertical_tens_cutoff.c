#include <math.h>
#include <stdio.h>

#include "../elastic_laws/hookes_law.h"
#include "../globals.h"
#include "../utils.h"

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
    (void)SCD;
    (void)SPD;

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
    if (*NPROPS < 4)
    {
        fprintf(
            stderr,
            "UMAT Error: NPROPS < 3. Requires E, nu, tension_threshold and normal_axis_index\n");
        return;
    }
    if (*NSTATV < 1)
    {
        fprintf(stderr, "UMAT Error: NSTATV < 1. Requires at least 1 state variable.\n");
        return;
    }

    // Material Properties ---
    double E = PROPS[0];                  // Young's Modulus
    double nu = PROPS[1];                 // Poisson's Ratio
    double tension_threshold = PROPS[2];  // Tensile strength threshold in Y
    int normal_axis_index =
        (int)round(PROPS[3]);  // Normal axis for tension cutoff (1 for Y, 2 for Z)

    // Calculate Elastic Stiffness Matrix (DDSDDE_elastic) ---
    double DDSDDE_elastic[VOIGTSIZE_3D * VOIGTSIZE_3D];
    calculate_elastic_stiffness_matrix_3d(E, nu, DDSDDE_elastic);

    // Calculate Elastic Trial Stress ---
    // delta_stress = D * DSTRAN
    double delta_stress[VOIGTSIZE_3D];
    matrix_vector_multiply(DDSDDE_elastic, DSTRAN, VOIGTSIZE_3D, delta_stress);

    // stress_trial = STRESS + delta_stress
    double stress_trial[VOIGTSIZE_3D];
    add_vectors(STRESS, delta_stress, VOIGTSIZE_3D, stress_trial);

    // Apply Tension Cutoff Logic ---
    // Check normal stress component
    const double sigma_normal_trial = stress_trial[normal_axis_index];

    if (sigma_normal_trial > tension_threshold)
    {
        //  Update Stress, Jacobian, and State Variables ---
        // Tension cutoff is active: set vertical stress to the threshold
        STRESS[normal_axis_index] = tension_threshold;

        // set stiffness matrix terms to a small value for numerical stability
        // DDSDDE = DDSDDE_elastic * small_value
        double small_value = 1.0e-16;
        vector_scalar_multiply(DDSDDE_elastic, small_value, VOIGTSIZE_3D * VOIGTSIZE_3D, DDSDDE);

        STATEV[0] = 1.0;  // Indicate cutoff state
    }
    else
    {
        // Elastic behavior: Use the trial stress and elastic Jacobian

        // Copy the trial stress to the output stress
        copy_array(stress_trial, VOIGTSIZE_3D, STRESS);

        // Copy the elastic stiffness matrix to the output Jacobian
        copy_array(DDSDDE_elastic, VOIGTSIZE_3D * VOIGTSIZE_3D, DDSDDE);

        STATEV[0] = 0.0;  // Indicate elastic state
    }

    return;
}
