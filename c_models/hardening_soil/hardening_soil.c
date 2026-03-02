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

double calculate_stress_dependency_factor(const double c, const double phi_rad, const double p_ref, const double sigma_3, const double m);
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
    // --- 0. Initialization and Material Properties ---
    int i;
    int n_dim_dir = *NDI;
    int n_dim_shr = *NSHR;

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
    double c = props[0];  // Cohesion (not used in this example)
    double phi_deg = props[1];  // Friction angle (not used in this example)
    double psi_deg = props[2];  // Dilation angle (not used in this example)
    double tension_threshold = props[3];  // Tensile strength threshold
    double E_50_ref = props[4];  // Reference secant stiffness in drained triaxial test
    double E_oed_ref = props[5];   // Reference tangent stiffness for primary loading in oedometer test
    double E_ur_ref = props[6];  // Reference unload reload stiffness
    double m = props[7];  // power of stress level dependency of stiffness
    double nu = props[8];  // Poisson's Ratio
    double p_ref = props[9];  // reference pressure for stiffness
    double K0nc = props[10];  // K0 value for normal consolidation
    double Rf = props[11];  // failure ratio qf/qa


    // Convert angles to radians
    double phi_rad = phi_deg * PI / 180.0;
    double psi_rad = psi_deg * PI / 180.0;
    double sin_phi = sin(phi_rad);
    double cos_phi = cos(phi_rad);
    double sin_psi = sin(psi_rad);

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

    double principle_stresses[3];

    calculate_principle_stresses(STRESS, principle_stresses);

    double sigma_1 = principle_stresses[0];
    double sigma_2 = principle_stresses[1];
    double sigma_3 = principle_stresses[2];

    double stress_dependency_factor = calculate_stress_dependency_factor(c, phi_rad, p_ref, sigma_3, m);
    double E_50 = E_50_ref * stress_dependency_factor; // confining stiffness for primary loading
    double E_oed = E_oed_ref * stress_dependency_factor;
    double E_ur = E_ur_ref * stress_dependency_factor; // unloading/reloading stiffness


    // initialization
    // run mohr coulomb

    // --- 1. Calculate Elastic Stiffness Matrix ---
    calculate_elastic_stiffness_matrix_3d(E_mod, nu, Ce_matrix);

    // Initialize Jacobian DDSDDE to elastic matrix (default assumption)
    copy_array(Ce_matrix, VOIGTSIZE_3D * VOIGTSIZE_3D, DDSDDE);

    // --- 2. Elastic Predictor Step ---
    // stress_trial = STRESS_n + Ce * DSTRAN
    matrix_vector_multiply(Ce_matrix, DSTRAN, VOIGTSIZE_3D, delta_stress);

    add_vectors(STRESS, delta_stress, VOIGTSIZE_3D,
                stress_trial);  // stress_trial = STRESS + delta_stress

    double mean_stress = (sigma_1 + sigma_2 + sigma_3) / 3.0;
    double q = sigma_1 - sigma_3;

        // ultimate deviatoric stress
    double q_f = calculate_ultimate_deviatoric_stress(c, phi_rad, mean_stress);
    double q_a = Rf * q_f; // hardening parameter, approximated as a fraction of the ultimate stress

    double plastic_strain_1 = STATEV[0]; // Assuming plastic_strain_1 is stored in the first state variable
    double plastic_strain_2 = STATEV[1]; // Assuming plastic_strain_2 is stored in the second state variable
    double plastic_strain_3 = STATEV[2]; // Assuming plastic_strain_3 is stored in the third state variable

    double gamma_ps = calculate_shear_hardening_parameter(plastic_strain_1, plastic_strain_2, plastic_strain_3)

    f = shear_yield_function(q_a,  q, E_50, E_ur, gamma_ps)
    res = f;

    int i = 0;
    while (i < max_iter && abs(res) > ZERO_TOL)
    {
        // calculate cone parameters
        parameters = calculate_cone_parameters();

        double dgdsigma[VOIGTSIZE_3D];

        derivative_shear_yield_function(const double q_a, const double q, const double E_50, const double E_ur, const double gamma_ps, dgdsigma)

        dgdsigma = calculate_cone_dgdsigma();
        // normalize


        dfdh = calculate_dfdh_cone();
        dhdlamda = calculate_dhdlambda();

        dfdsigma = calculate_dfdsigma();

        matrix_vector_multiply(D, dgdsigma, VOIGTSIZE_3D, Ddgdsigma);

        double dfdsigmaDdgdsigma =  vector_dot_product(dfdsigma, Ddgdsigma, VOIGTSIZE_3D);

        dfdlambda = dfdsigmaDdgdsigma - dfdh * dhdlamda;

        // implicit
        d_lambda = d_lambda + res /(dfdlambda); // + viscosity / (d_t + SMALL_VALUE));

        d_eps_p = d_lambda * dgdsigma;
        d_eps_p_vol = calculate_volumetric_strain_3d(d_eps_p);
        eps_p_vol = eps_p_vol_state + d_eps_p_vol;

        princ_stress = princ_stress_elast - d_lambda * Ddgdsigma;
        // resort princ stress
        sort(princ_stress, VOIGTSIZE_3D);

        f = cone_yield_function();
        res = calculate_residual(f);

        i++;

    }



//    // yield function greater than zero, calculate plastic correction
//    if (f_trial > ZERO_TOL)
//    {
//        // gradient yield function
//        calculate_stress_invariants_derivatives_3d(J_trial, s_dev, j2_trial, j3_trial, dp_dsig,
//                                                   dJ_dsig, dtheta_dsig);
//
//        double mats_nak_constants[5] = {phi_rad, M, alpha, beta, gamma};
//
//        calculate_yield_gradient(theta_trial, J_trial, mats_nak_constants, dp_dsig, dJ_dsig,
//                                 dtheta_dsig, grad_f);
//
//        // gradient potential function, g, it is required to recalculate the constants using psi
//        calculate_matsuoka_nakai_constants(psi_rad, c, &alpha, &beta, &gamma, &K, &M);
//        double mats_nak_constants_plastic_potential[5] = {psi_rad, M, alpha, beta, gamma};
//
//        calculate_yield_gradient(theta_trial, J_trial, mats_nak_constants_plastic_potential,
//                                 dp_dsig, dJ_dsig, dtheta_dsig, grad_g);
//
//        // Calculate terms needed for delta_gamma and Jacobian
//        matrix_vector_multiply(Ce_matrix, grad_g, VOIGTSIZE_3D, Ce_grad_g);  // Ce * g
//        matrix_vector_multiply(Ce_matrix, grad_f, VOIGTSIZE_3D, Ce_grad_f);  // Ce * f
//
//        // Calculate denominator for delta_gamma (and Jacobian)
//        // Assumes perfect plasticity (Hardening modulus H=0)
//        // Denom = A_vec : Ce : g_vec = grad_f . (Ce * grad_g)
//        double denom = vector_dot_product(grad_f, Ce_grad_g, VOIGTSIZE_3D);
//        double delta_gamma = f_trial / denom;
//
//        // --- Update Stress ---
//        // STRESS_{n+1} = stress_trial - delta_gamma * Ce * g_vec
//        for (i = 0; i < VOIGTSIZE_3D; ++i)
//        {
//            STRESS[i] = stress_trial[i] - delta_gamma * Ce_grad_g[i];
//        }
//
//        // dEps_p = delta_gamma * g_vec
//        vector_scalar_multiply(grad_g, delta_gamma, VOIGTSIZE_3D, dEps_p);
//
//        double dSpd = vector_dot_product(STRESS, dEps_p, VOIGTSIZE_3D);
//        *SPD += dSpd;
//
//        // Calculate outer product of Ce_grad_g and Ce_grad_f
//        const double Ce_grad_g_outer_Ce_grad_f[VOIGTSIZE_3D * VOIGTSIZE_3D];
//        vector_outer_product(Ce_grad_g, Ce_grad_f, VOIGTSIZE_3D, Ce_grad_g_outer_Ce_grad_f);
//
//        // --- Calculate Consistent Tangent Modulus (Jacobian DDSDDE) ---
//        // DDSDDE = Ce - (Ce * g) outer_prod (Ce * f)^T / denom
//        if (fabs(denom) > ZERO_TOL)
//        {  // Redundant check, but safe
//            for (i = 0; i < VOIGTSIZE_3D * VOIGTSIZE_3D; ++i)
//            {
//                // Accessing 1D array DDSDDE with row-major logic
//                DDSDDE[i] = Ce_matrix[i] - Ce_grad_g_outer_Ce_grad_f[i] / denom;
//            }
//        }
//
//        STATEV[0] = 1.0;  // Indicate that yield surface was reached
//    }
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


double calculate_stress_dependency_factor(const double c, const double phi_rad, const double p_ref, const double sigma_3, const double m)
{
    const double c_cos_phi = c * cos(phi_rad);
    double factor = pow((c_cos_phi - sigma_3*sin(phi_rad))/(c_cos_phi + p_ref * sin(phi_rad)),m);
    return factor;
}


double calculate_ultimate_deviatoric_stress(const double c, const double phi_rad, const double mean_stress)
{
    // follows from mohr coulomb yield criterion
    double q_f = 6 * sin(phi_rad) / (3- sin(phi_rad)) * (mean_stress +  c /tan(phi_rad));

    return q_f;
}

void calculate_stress_on_yield_surface()
{
    double f_tension = get_tension_yield_value(tension_cutoff, princ_stress)

    double f_cone = get_cone_yield_value(cone_cutoff, princ_stress)

//    double f_complimentary = get_complementary_yield_value(complementary_cutoff, princ_stress)

    double f_cap = cap_yield_function(q_special, M, p, pc)


     if (f_cone > ZERO_TOL)
     {
         if (f_tension > ZERO_TOL)
         {
             // first cone hardening
             int is_converged = calculate_cone_hardening();

             if (is_converged)
             {

                 double pc_state_var = pc; // todo change this

                 double f_r_cap = cap_yield_function(q_special, M, p, pc_state_var)

                 if (f_r_cap > ZERO_TOL)
                 {
                      if (f_cap > ZERO_TOL)
                      {

                           is_converged = calculate_cone_and_cap_hardening();

                      }
                      else {
                        is_converged = calculate_cap_hardening();
                      }
                      if (is_converged)
                      {
                          double f_r_tension = get_tension_yield_value(tension_cutoff, princ_stress)
                          if (f_r_tension > ZERO_TOL)
                          {
                             is_converged = 0;
                             int is_tension = 1;
                          }

                      }
                 }
                 return
             }
             if (is_tension)
             {
                 double c_cot_phi = c / tan(phi_rad);
                 if (c_cot_phi > ZERO_TOL && c_cot_phi > tension_cutoff)
                 {
                    is_converged = calculate_cone_hardening_with_tension();
                    if (is_converged)
                    {
                        return
                    }
                    else
                    {
                        is_converged = tension_criterium();
                        return
                    }
                 }
                 else
                 {
                     is_converged = tension_criterium();
                     return
                 }
             }
             if (!is_converged && !is_tension)
             {
                fprintf(stderr, "UMAT Error: (no tension or convergence, maybe close to c_cot_phi?).\n");
             }





         }
         else
         {
            is_converged = calculate_cone_hardening();

            if (!is_converged && !is_tension)
            {
                fprintf(stderr, "UMAT Error: (no tension or convergence, maybe close to c_cot_phi?).\n");
            }

            // check cap later
         }
     }
     else
     {
        if(f_cap > ZERO_TOL)
        {
            is_converged = calculate_cap_hardening();

            if (is_converged)
            {
               double f_r_tension = get_tension_yield_value(tension_cutoff, princ_stress)
               if(f_r_tension > ZERO_TOL)
               {
                  is_converged = 0;
                  int is_tension = 1;
                  return
               }

            }

            if (is_converged)
            {
                double f_r_cone = get_cone_yield_value(cone_cutoff, princ_stress)

                if (f_r_cone > ZERO_TOL)
                {
                    is_converged = calculate_cone_and_cap_hardening();
                    if (is_converged)
                    {
                        return;
                    }
                    else
                    {
                        is_converged = cone_hardening();
                        return;
                    }
                }
                else
                {
                    return;
                }
            }
            else
            {
                return;
            }

        }
        else if (f_tension > ZERO_TOL)
        {
            is_converged = tension_criterium();
            return;
        }

     }
     else{
         // elastic
         is_converged = 1;
         is_tension = 0;
         return;

     }

     if (!is_converged)
     {
        if (is_tension)
         {
             is_converged = tension_criterium();
             return;
         }
     }
     else{
         // check cap

         f_r_cap = cap_yield_function(q_special, M, p, pc_state_var);
         if (f_r_cap > ZERO_TOL)
         {
             if (f_cap > ZERO_TOL)
             {
                 if ( f_cone > ZERO_TOL)
                 {
                    is_converged = calculate_cone_and_cap_hardening();

                    if ( is_converged)
                    {
                        f_r_tension = get_tension_yield_value(tension_cutoff, princ_stress)
                        if (f_r_tension > ZERO_TOL)
                        {
                            is_converged = 0;
                            is_tension = 1;
                            return
                        }
                    }
                    else
                    {
                        is_converged = calculate_cap_hardening();
                    }
                 }
             }
             else
             {
                is_converged = calculate_cap_hardening();
             }
         }

     }



}