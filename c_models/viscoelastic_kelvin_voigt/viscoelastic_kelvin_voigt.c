#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // For exit()

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
//
// This is the implementation for the Kelvin-Voigt viscoelastic model for 3D problems.
// The model is formulated in total stress and total strain, with viscous effects
// incorporated via strain rate dependence.
// Kelvin-Voigt model: σ = E·ε + η·(dε/dt)
//
// State Variables (STATEV):
//   STATEV[0-5]: Elastic stress components (σ_xx, σ_yy, σ_zz, σ_xy, σ_xz, σ_yz)
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
    (void)STRAN;
    (void)NSTATV;
    (void) KINC;
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

    if (*NSTATV != 6)
    {
        fprintf(stderr, "UMAT Error: NSTATV != 6. Requires exactly 6 state variables,\n"
            "             corresponding to the elastic stress components.\n");
        return;
    }

    // Material Properties
    const double E = PROPS[0];     // Elastic Young's Modulus
    const double nu = PROPS[1];    // Poisson's Ratio
    const double eta = PROPS[2];   // Viscosity coefficient

    double DDSDDE_elastic[VOIGTSIZE_3D * VOIGTSIZE_3D];
    double elastic_stress_increment[VOIGTSIZE_3D];
    double new_elastic_stress[VOIGTSIZE_3D];
    double viscous_stress_increment[VOIGTSIZE_3D];
    double delta_sigma[VOIGTSIZE_3D];
    double elastic_stress_old[VOIGTSIZE_3D];

    // Retrieve old elastic stress from STATEV
    for (int i=0; i<VOIGTSIZE_3D; i++) {
        elastic_stress_old[i] = STATEV[i];
    }

    // Calculate Elastic Stiffness Matrix (DDSDDE_elastic)
    calculate_elastic_stiffness_matrix_3d(E, nu, DDSDDE_elastic);

    // Calculate Elastic Stress Increment
    // elastic_stress_increment = D * DSTRAN
    matrix_vector_multiply(DDSDDE_elastic, DSTRAN, VOIGTSIZE_3D, elastic_stress_increment);

    // Calculate the New Elastic Stress
    // new_elastic_stress = elastic_stress_old + (D * DSTRAN)
    add_vectors(elastic_stress_old, elastic_stress_increment, VOIGTSIZE_3D, new_elastic_stress);

    // Calculate Viscous Stress Increment
    // viscous_stress_increment = eta / DTIME * DSTRAN
    vector_scalar_multiply(DSTRAN, eta / (*DTIME), VOIGTSIZE_3D, viscous_stress_increment);

    // Total Stress
    add_vectors(new_elastic_stress, viscous_stress_increment, VOIGTSIZE_3D, STRESS);

    // Consistent tangent: DDSDDE_elastic + eta/DTIME * I
    copy_array(DDSDDE_elastic, VOIGTSIZE_3D * VOIGTSIZE_3D, DDSDDE);
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        DDSDDE[i * VOIGTSIZE_3D + i] += eta / (*DTIME);
    }

    // update STATEV with New Elastic Stress
    copy_array(new_elastic_stress, VOIGTSIZE_3D, STATEV);

    // Energy Calculations
    *SSE += vector_dot_product(elastic_stress_increment, DSTRAN, VOIGTSIZE_3D);
    *SPD += vector_dot_product(viscous_stress_increment, DSTRAN, VOIGTSIZE_3D);

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