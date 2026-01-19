#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // For exit()

#include "../elastic_laws/hookes_law.h"
#include "../viscoelastic_laws/viscoelastic_kelvin_voigt.h"
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
    // Outputs
    double* STRESS,  // Stress tensor at end of increment (NTENS components)
    double* STATEV,  // State variables at end of increment (NSTATV components)
    double* DDSDDE,  // Jacobian matrix (NTENS * NTENS components)
    double* SSE,     // Elastic strain energy density
    double* SPD,     // Viscous dissipation (formerly plastic)
    double* SCD,     // Creep dissipation (unused here, using SPD for total dissipation)
    double* RPL,     // Volumetric heat generation (unused)
    double* DDSDDT,  // Stress rate dependency on temperature (unused)
    double* DRPLDE,  // Derivative of RPL wrt strain (unused)
    double* DRPLDT,  // Derivative of RPL wrt temperature (unused)
    // Inputs
    const double* STRAN,   // Total strain at START of increment (NTENS components)
    const double* DSTRAN,  // Increment in total strain (NTENS components)
    const double* TIME,    // Step time [0] and total time [1]
    const double* DTIME,   // Time increment
    const double* TEMP,    // Temperature at start of increment (unused)
    const double* DTEMP,   // Increment in temperature (unused)
    const double* PREDEF,  // Predefined field variables (unused)
    const double* DPRED,   // Increment in predefined field variables (unused)
    char* CMNAME,          // Material name (unused)
    const int* NDI,        // Number of direct stress components (3)
    const int* NSHR,       // Number of shear stress components (3)
    const int* NTENS,      // Total number of stress components (6)
    const int* NSTATV,     // Number of state variables
    const double* PROPS,   // User-defined material properties (NPROPS components)
    const int* NPROPS,     // Number of properties
    const double* COORDS,  // Coordinates of the integration point (unused)
    const double* DROT,    // Rotation increment matrix (unused)
    double* PNEWDT,        // Suggested new time increment size
    const double* CELENT,  // Characteristic element length (unused)
    const double* DFGRD0,  // Deformation gradient at start (unused)
    const double* DFGRD1,  // Deformation gradient at end (unused)
    const int* NOEL,       // Element number (unused)
    const int* NPT,        // Integration point number (unused)
    const int* LAYER,      // Layer number (unused)
    const int* KSPT,       // Section point number (unused)
    const int* KSTEP,      // Step number (unused)
    const int* KINC        // Increment number (unused)
)
{
    // Avoid unused variable warnings
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
    (void)TIME;
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

    // Material Properties
    // PROPS[0] and PROPS[1] are E and nu (elastic)
    // PROPS[2] and PROPS[3] are Ev and nu_v (viscous)
    const double E = PROPS[0];     // Elastic Young's Modulus
    const double nu = PROPS[1];    // Poisson's Ratio
    const double eta = PROPS[2];    // Viscosity coefficent (porportional to the strain rate)

    double DDSDDE_elastic[VOIGTSIZE_3D * VOIGTSIZE_3D];
    double delta_elastic_stress[VOIGTSIZE_3D];
    double delta_viscous_stress[VOIGTSIZE_3D];
    double updated_stress[VOIGTSIZE_3D];
    double elastic_delta_strain_vector[VOIGTSIZE_3D];

    // Calculate Elastic Stiffness Matrix (DDSDDE_elastic)
    calculate_elastic_stiffness_matrix_3d(E, nu, DDSDDE_elastic);

    // Calculate Elastic Stress
    // delta_stress = D * DSTRAN
    matrix_vector_multiply(DDSDDE_elastic, DSTRAN, VOIGTSIZE_3D, delta_elastic_stress);

    // Viscous Stress Update (Kelvin-Voigt Model)
    // delta_viscous_stress; += eta * (DSTRAN / DTIME)
    vector_scalar_multiply(DSTRAN, eta / *DTIME, VOIGTSIZE_3D, delta_viscous_stress);

    // add viscous contribution
    // updated_stress = delta_stress
    add_vectors(delta_elastic_stress, delta_viscous_stress, VOIGTSIZE_3D, updated_stress);


    // Copy the trial stress to the output stress
    copy_array(updated_stress, VOIGTSIZE_3D, STRESS);

    // Copy the elastic stiffness matrix to the output Jacobian
    copy_array(DDSDDE_elastic, VOIGTSIZE_3D * VOIGTSIZE_3D, DDSDDE);

    *SSE += vector_dot_product(STRESS, DSTRAN, VOIGTSIZE_3D);  // Specific elastic strain energy

    return;
}

int check_properties(const int NPROPS, const double* PROPS)
{
    // Now requires 3 properties: E, nu, eta
    if (NPROPS < 3)
    {
        fprintf(stderr, "UMAT Error: NPROPS < 3. Requires E, nu, and eta.\n");
        return 1;
    }

    int n_errors = 0;

    // Elastic checks
    if (PROPS[0] <= 0.0)
    { /* E > 0 */
        n_errors++;
    }
    if (PROPS[1] < 0.0 || PROPS[1] >= 0.5)
    { /* nu in [0, 0.5) */
        n_errors++;
    }

    // Viscous checks (eta must be greater than or equal to zero)
    if (PROPS[2] < 0.0)  // eta can be zero (pure elastic) or positive
    {
        fprintf(stderr, "UMAT Error: Viscosity coefficient (eta) must be non-negative.\n");
        n_errors++;
    }

    if (n_errors > 0)
    {
        fprintf(stderr, "UMAT Error: %d property errors found.\n", n_errors);
        return 1;
    }
    return 0;
}