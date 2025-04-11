
#include <math.h>
#include "matsuoka_nakai_surface.h"

#define PI 3.14159265358979323846
#define ZERO_TOL 1.0e-12 // Tolerance for zero checks
#define VOIGTSIZE_3D 6

// Define tensor component mapping (Voigt, 0-based)
#define XX 0
#define YY 1
#define ZZ 2
#define XY 3
#define YZ 4
#define XZ 5



void calculate_yield_function(double p, double theta,double J,
                              double c, double phi_rad, double alpha, double beta, double gamma, double K, double M, double* f)
{

        *f  = (-K + M * p) + J * alpha * cos(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0);

        // in paper:
        //*f  = -(K + M * p) + J * alpha * cos(acos(beta * sin(3 * theta)) / 3 - gamma * PI / 6);

 }

void calculate_yield_gradient(double theta,
                               double J,
                               double* constants,
                               double dp_dsig[VOIGTSIZE_3D], double dJ_dsig[VOIGTSIZE_3D], double dtheta_dsig[VOIGTSIZE_3D], double grad[VOIGTSIZE_3D]) // Output gradient vector
{

    double angle_rad = constants[0]; // Angle in radians
    double M = constants[1]; // Matsuoka-Nakai constant
    double alpha = constants[2]; // Matsuoka-Nakai constant
    double beta = constants[3]; // Matsuoka-Nakai constant
    double gamma = constants[4]; // Matsuoka-Nakai constant

    // Computes df/dsigma or dg/dsigma based on the angle provided
    // Uses chain rule: grad = (dF/dp)*(dp/dsig) + (dF/dq)*(dq/dsig) + (dF/dtheta)*(dtheta/dsig)
    // where F is the yield function f or potential g, using angle_rad.

    // Initialize gradient to zero
    for (int i = 0; i < VOIGTSIZE_3D; ++i) grad[i] = 0.0;

    // --- Components for Chain Rule ---
    double sin_angle = sin(angle_rad);
    double cos_angle = cos(angle_rad); // Needed if using the c*cos(angle) form

    // 1. Derivatives of F (yield/potential function) w.r.t. invariants


    double dF_dp = M; // in paper : dF_dp = -M;
    double C1 = alpha * cos(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0);
    double dF_dJ = C1;

    double dC1_dtheta = alpha * beta * cos(3 * theta) * sin(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0) / sqrt(1.0 - beta*beta * pow(sin(3.0 * theta),2.0));
    double dF_dtheta = J * dC1_dtheta;

    // 3. Combine using chain rule: grad = dF_dp*dp_dsig + dF_dq*dq_dsig + dF_dtheta*dtheta_dsig
    for (int i = 0; i < VOIGTSIZE_3D; ++i) {
        grad[i] = dF_dp * dp_dsig[i] + dF_dJ * dJ_dsig[i] + dF_dtheta * dtheta_dsig[i];
    }
}




void calculate_matsuoka_nakai_constants(double angle_rad, double c, double *alpha, double *beta, double *gamma, double *K, double *M)
{
        //if angle_rad is zero, return tresca constants
        if (angle_rad < ZERO_TOL)
        {
            *M = 0;
            *K = 0;
            *alpha = 1 / (cos(PI / 6));
            *beta = 0.9999;
            *gamma = 1;
        }
        // else Matsuaoka Nakai constants
        else{

            *M = 1.0 / sqrt(3.0) * 6.0 * sin(angle_rad) / (3.0 - sin(angle_rad));
            *K = c / tan(angle_rad);
            double k_mn = (9.0 -pow(sin(angle_rad),2)) / (1.0 - pow(sin(angle_rad),2.0));
            double A1 = (k_mn - 3.0) / (k_mn - 9.0);
            double A2 = k_mn / (k_mn - 9.0);

            *alpha = 2.0 / sqrt(3.0) * sqrt(A1) * *M;
            *beta = A2 / (pow(A1,1.5));
            *gamma = 0.0;

        }
}