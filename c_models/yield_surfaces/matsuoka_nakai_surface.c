
#include <math.h>

#include "../globals.h"
#include "matsuoka_nakai_surface.h"

void calculate_yield_function(const double p, const double theta, const double J,
                              const MatsuokaNakaiConstants constants, double* f)
{
    *f = constants.M * (-constants.K + p) +
         J * constants.alpha *
             cos(acos(constants.beta * sin(3.0 * theta)) / 3.0 - constants.gamma * PI / 6.0);

    // in paper:
    //*f  = -(K + M * p) + J * alpha * cos(acos(beta * sin(3 * theta)) / 3 - gamma * PI / 6);
}

void calculate_yield_gradient(const double theta, const double J,
                              const MatsuokaNakaiConstants constants,
                              const double dp_dsig[VOIGTSIZE_3D],
                              const double dJ_dsig[VOIGTSIZE_3D],
                              const double dtheta_dsig[VOIGTSIZE_3D],
                              double grad[VOIGTSIZE_3D])  // Output gradient vector
{
    // Computes df/dsigma or dg/dsigma based on the angle provided
    // Uses chain rule: grad = (dF/dp)*(dp/dsig) + (dF/dq)*(dq/dsig) + (dF/dtheta)*(dtheta/dsig)
    // where F is the yield function f or potential g, using angle_rad.

    // Initialize gradient to zero
    for (int i = 0; i < VOIGTSIZE_3D; ++i) grad[i] = 0.0;

    // --- Components for Chain Rule ---

    // Derivatives of F (yield/potential function) w.r.t. invariants

    const double dF_dp = constants.M;  // in paper : dF_dp = -M;
    const double C1 = constants.alpha * cos(acos(constants.beta * sin(3.0 * theta)) / 3.0 -
                                            constants.gamma * PI / 6.0);
    const double dF_dJ = C1;

    const double dC1_dtheta =
        constants.alpha * constants.beta * cos(3.0 * theta) *
        sin(acos(constants.beta * sin(3.0 * theta)) / 3.0 - constants.gamma * PI / 6.0) /
        sqrt(1.0 - constants.beta * constants.beta * pow(sin(3.0 * theta), 2.0));
    const double dF_dtheta = J * dC1_dtheta;

    // Combine using chain rule: grad = dF_dp*dp_dsig + dF_dq*dq_dsig + dF_dtheta*dtheta_dsig
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        grad[i] = dF_dp * dp_dsig[i] + dF_dJ * dJ_dsig[i] + dF_dtheta * dtheta_dsig[i];
    }
}

MatsuokaNakaiConstants calculate_matsuoka_nakai_constants(const double angle_rad, const double c)
{
    double M;
    double K;
    double alpha;
    double beta;
    double gamma;

    // if angle_rad is zero, return tresca constants
    if (fabs(angle_rad) < ZERO_TOL)
    {
        M = 0.0;
        K = 0.0;
        alpha = 1.0 / (cos(PI / 6.0));
        beta = 0.9999;
        gamma = 1.0;
    }
    // else Matsuaoka Nakai constants
    else
    {
        M = 6.0 / sqrt(3.0) * sin(angle_rad) / (3.0 - sin(angle_rad));
        K = c / tan(angle_rad);

        const double k_mn = (9.0 - pow(sin(angle_rad), 2)) / (1.0 - pow(sin(angle_rad), 2.0));
        const double A1 = (k_mn - 3.0) / (k_mn - 9.0);
        const double A2 = k_mn / (k_mn - 9.0);

        alpha = 2.0 / sqrt(3.0) * sqrt(A1) * M;

        beta = A2 / (pow(A1, 1.5));
        gamma = 0.0;
    }

    MatsuokaNakaiConstants constants = {alpha, beta, gamma, K, M};

    return constants;
}
