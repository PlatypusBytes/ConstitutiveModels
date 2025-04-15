
#include <math.h>
#include <stdio.h>

#include "../globals.h"
#include "matsuoka_nakai_surface.h"

void calculate_yield_function(const double p, const double theta, const double J, const double c,
                              const double phi_rad, const double alpha, const double beta,
                              const double gamma, const double K, const double M, double* f)
{
    *f = M * (-K + p) + J * alpha * cos(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0);

    // in paper:
    //*f  = -(K + M * p) + J * alpha * cos(acos(beta * sin(3 * theta)) / 3 - gamma * PI / 6);
}

void calculate_yield_gradient(const double theta, const double J, const double* constants,
                              const double dp_dsig[VOIGTSIZE_3D],
                              const double dJ_dsig[VOIGTSIZE_3D],
                              const double dtheta_dsig[VOIGTSIZE_3D],
                              double grad[VOIGTSIZE_3D])  // Output gradient vector
{
    const double angle_rad = constants[0];  // Angle in radians
    const double M = constants[1];          // Matsuoka-Nakai constant
    const double alpha = constants[2];      // Matsuoka-Nakai constant
    const double beta = constants[3];       // Matsuoka-Nakai constant
    const double gamma = constants[4];      // Matsuoka-Nakai constant

    // Computes df/dsigma or dg/dsigma based on the angle provided
    // Uses chain rule: grad = (dF/dp)*(dp/dsig) + (dF/dq)*(dq/dsig) + (dF/dtheta)*(dtheta/dsig)
    // where F is the yield function f or potential g, using angle_rad.

    // Initialize gradient to zero
    for (int i = 0; i < VOIGTSIZE_3D; ++i) grad[i] = 0.0;

    // --- Components for Chain Rule ---
    const double sin_angle = sin(angle_rad);
    const double cos_angle = cos(angle_rad);  // Needed if using the c*cos(angle) form

    // 1. Derivatives of F (yield/potential function) w.r.t. invariants

    const double dF_dp = M;  // in paper : dF_dp = -M;
    const double C1 = alpha * cos(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0);
    const double dF_dJ = C1;

    const double dC1_dtheta = alpha * beta * cos(3 * theta) *
                              sin(acos(beta * sin(3.0 * theta)) / 3.0 - gamma * PI / 6.0) /
                              sqrt(1.0 - beta * beta * pow(sin(3.0 * theta), 2.0));
    const double dF_dtheta = J * dC1_dtheta;

    // 3. Combine using chain rule: grad = dF_dp*dp_dsig + dF_dq*dq_dsig + dF_dtheta*dtheta_dsig
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        grad[i] = dF_dp * dp_dsig[i] + dF_dJ * dJ_dsig[i] + dF_dtheta * dtheta_dsig[i];
    }
}

void calculate_matsuoka_nakai_constants(const double angle_rad, const double c, double* alpha,
                                        double* beta, double* gamma, double* K, double* M)
{
    // if angle_rad is zero, return tresca constants
    if (angle_rad < ZERO_TOL)
    {
        *M = 0;
        *K = 0;
        *alpha = 1 / (cos(PI / 6));
        *beta = 0.9999;
        *gamma = 1;
    }
    // else Matsuaoka Nakai constants
    else
    {
        *M = 6.0 / sqrt(3.0) * sin(angle_rad) / (3.0 - sin(angle_rad));
        *K = c / tan(angle_rad);

        const double k_mn = (9.0 - pow(sin(angle_rad), 2)) / (1.0 - pow(sin(angle_rad), 2.0));
        const double A1 = (k_mn - 3.0) / (k_mn - 9.0);
        const double A2 = k_mn / (k_mn - 9.0);

        *alpha = 2.0 / sqrt(3.0) * sqrt(A1) * *M;

        *beta = A2 / (pow(A1, 1.5));
        *gamma = 0.0;
    }
}
