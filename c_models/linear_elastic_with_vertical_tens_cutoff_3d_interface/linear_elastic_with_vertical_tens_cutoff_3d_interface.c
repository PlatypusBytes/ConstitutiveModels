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

/**
 * @brief Function to check the validity of material properties.
 *
 * @param[in]  NPROPS Number of properties.
 * @param[in]  PROPS Pointer to the array of material properties.
 * @return 0 if properties are valid, 1 otherwise.
 */
int check_properties(const int NPROPS, const double* PROPS);

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
    (void)SCD;

    // Check Inputs
    if (*NTENS != VOIGTSIZE_3D_INTERFACE || *NDI != 1 || *NSHR != 2)
    {
        // Handle error - This UMAT is specifically for 3D interface elements
        fprintf(stderr, "UMAT Error: NTENS != 3. This UMAT requires 3D interface elements.\n");
        // exit(1); // Avoid exiting in production code if possible
        return;  // Or try to handle gracefully
    }
    if (check_properties(*NPROPS, PROPS))
    {
        fprintf(stderr, "UMAT Error: invalid material properties.\n");
        return;
    }

    // state variable is 0 no tension cut-off, 1 tension cut-off
    if (*NSTATV < 1)
    {
        fprintf(stderr, "UMAT Error: NSTATV < 1. Requires at least 1 state variable.\n");
        return;
    }

    // Material Properties
    double E = PROPS[0];                  // Young's Modulus
    double nu = PROPS[1];                 // Poisson's Ratio
    double tension_threshold = PROPS[2];  // Tensile strength threshold in Y

    double DDSDDE_elastic[VOIGTSIZE_3D_INTERFACE * VOIGTSIZE_3D_INTERFACE];
    double delta_stress[VOIGTSIZE_3D_INTERFACE];
    double stress_trial[VOIGTSIZE_3D_INTERFACE];
    double elastic_delta_strain_vector[VOIGTSIZE_3D_INTERFACE];

    int normal_axis_index = 0;

    // Calculate Elastic Stiffness Matrix (DDSDDE_elastic)
    calculate_elastic_stiffness_matrix_3d_interface(E, nu, DDSDDE_elastic);

    // Calculate Elastic Trial Stress
    // delta_stress = D * DSTRAN
    matrix_vector_multiply(DDSDDE_elastic, DSTRAN, VOIGTSIZE_3D_INTERFACE, delta_stress);


    // stress_trial = STRESS + delta_stress
    add_vectors(STRESS, delta_stress, VOIGTSIZE_3D_INTERFACE, stress_trial);

    // Apply Tension Cutoff Logic
    // Check normal stress component
    const double sigma_normal_trial = stress_trial[normal_axis_index];

    if (sigma_normal_trial > tension_threshold)
    {
        // elastic normal stiffness
        const double elastic_normal_stiffness =
            DDSDDE_elastic[normal_axis_index * VOIGTSIZE_3D_INTERFACE + normal_axis_index];

        // calculate elastic and plastic strains
        const double dEps_el = (tension_threshold - STRESS[normal_axis_index]) / elastic_normal_stiffness;
        const double dEps_pl = DSTRAN[normal_axis_index] - dEps_el;

        // Copy the trial stress to the output stress
        copy_array(DSTRAN, VOIGTSIZE_3D_INTERFACE, elastic_delta_strain_vector);
        elastic_delta_strain_vector[normal_axis_index] = dEps_el;

        // calculate elastic and plastic strain energy
        *SSE += vector_dot_product(STRESS, elastic_delta_strain_vector,
                                   VOIGTSIZE_3D_INTERFACE);  // Specific elastic strain energy
        *SPD += tension_threshold * dEps_pl;       // Specific plastic strain dissipation

        //  Update Stress, Jacobian, and State Variables
        // Tension cutoff is active: set vertical stress to the threshold
        STRESS[normal_axis_index] = tension_threshold;

        // set stiffness matrix terms to a small value for numerical stability
        // DDSDDE = DDSDDE_elastic * small_value

        vector_scalar_multiply(DDSDDE_elastic, SMALL_VALUE, VOIGTSIZE_3D_INTERFACE * VOIGTSIZE_3D_INTERFACE, DDSDDE);

        STATEV[0] = 1.0;  // Indicate cutoff state
    }
    else
    {
        // Elastic behavior: Use the trial stress and elastic Jacobian

        // Copy the trial stress to the output stress
        copy_array(stress_trial, VOIGTSIZE_3D_INTERFACE, STRESS);

        // Copy the elastic stiffness matrix to the output Jacobian
        copy_array(DDSDDE_elastic, VOIGTSIZE_3D_INTERFACE * VOIGTSIZE_3D_INTERFACE, DDSDDE);

        *SSE += vector_dot_product(STRESS, DSTRAN, VOIGTSIZE_3D_INTERFACE);  // Specific elastic strain energy

        STATEV[0] = 0.0;  // Indicate elastic state
    }

    return;
}

int check_properties(const int NPROPS, const double* PROPS)
{
    if (NPROPS < 3)
    {
        fprintf(
            stderr,
            "UMAT Error: NPROPS < 3. Requires E, nu, tension_threshold\n");
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
        fprintf(stderr, "UMAT Error: tension_threshold must be non-negative.\n");
        n_errors++;
    }

    if (n_errors > 0)
    {
        return 1;
    }
    return 0;  // All checks passed
}
