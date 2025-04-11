#include <stdio.h>
#include <math.h>

#include "../utils.h"
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

#define ZERO_TOL 1.0e-10 // Tolerance for zero checks (q, denominators)
#define PI 3.14159265358979323846

// Define tensor component mapping (Voigt, 0-based)
#define XX 0
#define YY 1
#define ZZ 2
#define XY 3
#define YZ 4
#define XZ 5

void calculate_elastic_stiffness(double E, double nu, double* DDSDDE, int NTENS);

void calculate_stress_invariants(const double* stress, int NDI, int NSHR, int NTENS,
                                 double* p, double* q, double* theta,
                                 double* j2, double* j3, double* s_dev);

void calculate_yield_function(double p, double theta,double q,
                              double c, double phi_rad, double alpha, double beta, double gamma, double K, double M, double* f)    ;

void calculate_stress_gradient(const double* s_dev, double p, double q, double theta,
                               double j2, double j3,
                               double angle_rad, // phi for df/dsigma, psi for dg/dsigma
                               double M, double alpha, double beta, double gamma,
                               int NTENS, double* grad); // Output gradient vector

void calculate_matsuoka_nakai_constants(double *alpha, double *beta, double *gamma, double *K, double *M, double phi_rad, double c);

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
    double p_trial, q_trial, theta_trial, j2_trial, j3_trial;
    calculate_stress_invariants(stress_trial, n_dim_dir, n_dim_shr, n_tensor,
                                &p_trial, &q_trial, &theta_trial,
                                &j2_trial, &j3_trial, s_dev); // s_dev also calculated here


    // matsuoka nakai constants
    double alpha =0;
    double beta =0;
    double gamma =0;
    double K =0;
    double M =0;

    // yield function
    double f_trial = 0.0;

    // Calculate yield function value
    calculate_matsuoka_nakai_constants(&alpha, &beta, &gamma, &K, &M, phi_rad, c);
    calculate_yield_function(p_trial, theta_trial,q_trial, c, phi_rad,alpha,beta,gamma,K,M, &f_trial);

    // yield function greater than zero, calculate plastic correction
    if (f_trial > ZERO_TOL)
    {

        // gradient yield function
        calculate_stress_gradient(s_dev, p_trial, q_trial, theta_trial,j2_trial, j3_trial,phi_rad, // phi for df/dsigma, psi for dg/dsigma
                                   M, alpha, beta, gamma,
                                   n_tensor, grad_f);

        // gradient potential function, g, it is required to recalculate the constants using psi
        calculate_matsuoka_nakai_constants(&alpha, &beta, &gamma, &K, &M, psi_rad, c);
        calculate_stress_gradient(s_dev, p_trial, q_trial, theta_trial,j2_trial, j3_trial,psi_rad, // phi for df/dsigma, psi for dg/dsigma
                               M, alpha, beta, gamma, n_tensor, grad_g);

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

void calculate_stress_invariants(const double* stress, int NDI, int NSHR, int NTENS,
                                 double* p, double* q, double* theta,
                                 double* j2, double* j3, double* s_dev) // s_dev is output
{
    if (NTENS != 6) return;

    // Mean stress (pressure p = -trace(sigma)/3, but often p = trace(sigma)/3 in geomech)
    *p = (stress[XX] + stress[YY] + stress[ZZ]) / 3.0;

    // Deviatoric stress tensor s = sigma - p * I
    s_dev[XX] = stress[XX] - *p;
    s_dev[YY] = stress[YY] - *p;
    s_dev[ZZ] = stress[ZZ] - *p;
    s_dev[XY] = stress[XY];
    s_dev[YZ] = stress[YZ];
    s_dev[XZ] = stress[XZ];

    // Second invariant of deviatoric stress J2 = 0.5 * s:s
    // J2 = 0.5 * (sxx^2 + syy^2 + szz^2) + sxy^2 + syz^2 + sxz^2
     *j2 = 0.5 * (s_dev[XX] * s_dev[XX] + s_dev[YY] * s_dev[YY] + s_dev[ZZ] * s_dev[ZZ])
         + (s_dev[XY] * s_dev[XY] + s_dev[YZ] * s_dev[YZ] + s_dev[XZ] * s_dev[XZ]);

    // this is not the same as the von Mises stress, todo check this
    *q = sqrt((*j2));

    // Third invariant of deviatoric stress J3 = det(s)
    // J3 = sxx*syy*szz + 2*sxy*syz*sxz - sxx*syz^2 - syy*sxz^2 - szz*sxy^2

    calculate_determinant_voight_vector(s_dev, j3);

    // Lode angle theta
    // sin(3*theta) = - (J3 / 2) * (3 / J2)^(3/2) = - (3 * sqrt(3) / 2) * (J3 / q^3) if q!=0
    // or sin(3*theta) = - sqrt(27/8) * J3 / J2^(3/2)

    //double arg = (3.0 * sqrt(3.0) / 2.0) * (*j3) / pow(*q, 3.0);
    if (*q > ZERO_TOL) { // Avoid division by zero if q is very small
        //double arg = - (3.0 * sqrt(3.0) / 2.0) * (*j3) / pow(*q, 3.0);

        double arg = ( sqrt(27.0) / 2.0) * (*j3) / pow(*j2, 1.5);

//
        // Clamp argument to asin range [-1, 1] due to potential numerical inaccuracies
        if (arg > 1.0) arg = 1.0;
        if (arg < -1.0) arg = -1.0;

        *theta = asin(arg) / 3.0; // Result is in [-pi/6, pi/6]
    } else {
         *theta = 0.0; // Set to convention if q is zero
    }

    if (*theta < -(PI / 6.0) + ZERO_TOL){
    *theta = -(PI / 6.0) + ZERO_TOL;
    }
    else if (*theta > (PI / 6.0) - ZERO_TOL){
    *theta = (PI / 6.0) - ZERO_TOL;
    }

     // Ensure theta is within [-pi/6, pi/6] range
     // (asin should already handle this, but double check if implementing manually)
}


void calculate_yield_function(double p, double theta,double q,
                              double c, double phi_rad, double alpha, double beta, double gamma, double K, double M, double* f)
{

//        fprintf(stderr, "calculate_yield_function: p = %f, theta = %f, q = %f, c = %f, phi_rad = %f, alpha = %f, beta = %f, gamma = %f, K = %f, M = %f\n",
//                p, theta, q, c, phi_rad, alpha, beta, gamma, K, M);

//        self.yield_residual = (-self.K + self.M * self.stress_utils.mean_stress) + J * self.alpha * np.cos(
//            np.acos(self.beta * np.sin(3 * self.stress_utils.lode_angle)) / 3 - self.gamma * np.pi / 6)


        *f  = (-K + M * p) + q * alpha * cos(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0);

//        fprintf(stderr, "calculate_yield_function inside: f = %f\n", *f);

        // in paper:
        //*f  = -(K + M * p) + q * alpha * cos(acos(beta * sin(3 * theta)) / 3 - gamma * PI / 6);

 }



 void calculate_stress_gradient(const double* s_dev, double p, double q, double theta,
                               double j2, double j3,
                               double angle_rad, // phi for df/dsigma, psi for dg/dsigma
                               double M, double alpha, double beta, double gamma,
                               int NTENS, double* grad) // Output gradient vector
{
    if (NTENS != 6) return;

    // Computes df/dsigma or dg/dsigma based on the angle provided
    // Uses chain rule: grad = (dF/dp)*(dp/dsig) + (dF/dq)*(dq/dsig) + (dF/dtheta)*(dtheta/dsig)
    // where F is the yield function f or potential g, using angle_rad.

    // Initialize gradient to zero
    for (int i = 0; i < NTENS; ++i) grad[i] = 0.0;


    // --- Components for Chain Rule ---
    double sin_angle = sin(angle_rad);
    double cos_angle = cos(angle_rad); // Needed if using the c*cos(angle) form

    // 1. Derivatives of F (yield/potential function) w.r.t. invariants
    double C1 = alpha * cos(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0);

     double dF_dp = M;
//    double dF_dp = -M;
     double dF_dq = C1;

    // ai generated, check
    double dC1_dtheta = alpha * beta * cos(3 * theta) * sin(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0) / sqrt(1.0 - beta*beta * pow(sin(3.0 * theta),2.0));
    double dF_dtheta = q * dC1_dtheta;


    // 2. Derivatives of invariants w.r.t. sigma (in Voigt)
    // dp/dsigma = [1/3, 1/3, 1/3, 0, 0, 0]^T
    double dp_dsig[6] = {1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0, 0.0, 0.0};

    double dq_dsig[6];
    if (q < ZERO_TOL) {
        // Handle case where q is very small or zero
        for (int i = 0; i < NTENS; ++i) {
            dq_dsig[i] = 0.0; // Set to zero or some other value as needed
        }

    }
    else {
            dq_dsig[0] =s_dev[XX] / (2.0*q);
            dq_dsig[1] =s_dev[YY] / (2.0*q);
            dq_dsig[2] =s_dev[ZZ] / (2.0*q);
            dq_dsig[3] =s_dev[XY] / q;
            dq_dsig[4] =s_dev[YZ] / q;
            dq_dsig[5] =s_dev[XZ] / q;

    }


    // dtheta/dsigma = (dtheta/dJ2)*(dJ2/dsig) + (dtheta/dJ3)*(dJ3/dsig)
    // This is the complex part, requires dJ2/dsig, dJ3/dsig, dtheta/dJ2, dtheta/dJ3

    // dJ2/dsigma = s_dev (Voigt)
    //const double* dJ2_dsig = s_dev; // Pointer assignment is fine
    double dJ2_dsig[6] = {s_dev[XX], s_dev[YY], s_dev[ZZ], 2.0*s_dev[XY], 2.0*s_dev[YZ], 2.0*s_dev[XZ]};


    double dJ3_dsig[6];


    double term1 = s_dev[XX]*s_dev[XX] + s_dev[XY]*s_dev[XY] + s_dev[XZ]*s_dev[XZ] - 2.0*j2/3.0; // s_xx^2 + s_xy^2 + s_xz^2 - 2/3 * J2
    double term2 = s_dev[YY]*s_dev[YY] + s_dev[XY]*s_dev[XY] + s_dev[YZ]*s_dev[YZ] - 2.0*j2/3.0; // s_yy^2 + s_xy^2 + s_yz^2 - 2/3 * J2
    double term3 = s_dev[ZZ]*s_dev[ZZ] + s_dev[YZ]*s_dev[YZ] + s_dev[XZ]*s_dev[XZ] - 2.0*j2/3.0; // s_zz^2 + s_yz^2 + s_xz^2 - 2/3 * J2
    double term4 = 2.0*((s_dev[XX] + s_dev[YY])*s_dev[XY] + s_dev[XZ]*s_dev[YZ]); // 2*((s_xx+s_yy)*s_xy + s_xz*s_yz)
    double term5 = 2.0*((s_dev[YY] + s_dev[ZZ])*s_dev[YZ] + s_dev[XY]*s_dev[XZ]); // 2*((s_yy+s_zz)*s_yz + s_xy*s_xz)
    double term6 = 2.0*((s_dev[ZZ] + s_dev[XX])*s_dev[XZ] + s_dev[XY]*s_dev[YZ]); // 2*((s_zz+s_xx)*s_xz + s_xy*s_yz)

    dJ3_dsig[XX] = term1;
    dJ3_dsig[YY] = term2;
    dJ3_dsig[ZZ] = term3;
    dJ3_dsig[XY] = term4;
    dJ3_dsig[YZ] = term5;
    dJ3_dsig[XZ] = term6;

    double dtheta_dJ2 = 0.0;
    double dtheta_dJ3 = 0.0;

    double j2_pow_3 = pow(j2, 3);

    double inner_term = 1.0- (27.0 * j3 * j3) / (4.0 * pow(j2, 3.0));

    if (inner_term < ZERO_TOL) {
        dtheta_dJ2 = 0.0;
        dtheta_dJ3= 0.0;
    }
    else {
        dtheta_dJ2 = -pow(3.0, 1.5) * j3 / (4.0* pow(j2,2.5) * sqrt(inner_term));
        dtheta_dJ3 = sqrt(3.0) / sqrt(4.0*pow(j2, 3.0) -27.0 * j3 * j3) ;
    }

    // dtheta/dsigma = dtheta_dJ2 * dJ2_dsig + dtheta_dJ3 * dJ3_dsig
    double dtheta_dsig[6];
    for(int i=0; i<NTENS; ++i) {
        dtheta_dsig[i] = dtheta_dJ2 * dJ2_dsig[i] + dtheta_dJ3 * dJ3_dsig[i];
    }

    // 3. Combine using chain rule: grad = dF_dp*dp_dsig + dF_dq*dq_dsig + dF_dtheta*dtheta_dsig
    for (int i = 0; i < NTENS; ++i) {
        grad[i] = dF_dp * dp_dsig[i] + dF_dq * dq_dsig[i] + dF_dtheta * dtheta_dsig[i];
    }
}


void calculate_matsuoka_nakai_constants(double *alpha, double *beta, double *gamma, double *K, double *M, double phi_rad, double c)
{
        //if phi_rad is zero, return tresca constants
        if (phi_rad < ZERO_TOL)
        {
            *M = 0;
            *K = 0;
            *alpha = 1 / (cos(PI / 6));
            *beta = 0.9999;
            *gamma = 1;
        }
        // else Matsuaoka Nakai constants
        else{

            *M = 1.0 / sqrt(3.0) * 6.0 * sin(phi_rad) / (3.0 - sin(phi_rad));
            *K = c / tan(phi_rad);
            double k_mn = (9.0 -pow(sin(phi_rad),2)) / (1.0 - pow(sin(phi_rad),2.0));
            double A1 = (k_mn - 3.0) / (k_mn - 9.0);
            double A2 = k_mn / (k_mn - 9.0);

            *alpha = 2.0 / sqrt(3.0) * sqrt(A1) * *M;
            *beta = A2 / (pow(A1,1.5));
            *gamma = 0.0;

        }
}