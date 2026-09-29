#include "../globals.h"
#include "mohr_coulomb_surface.h"

void mohr_coulomb_principal_gradient(int i, int j, double sin_angle, double gradient[3])
{
    gradient[0] = gradient[1] = gradient[2] = 0.0;
    gradient[i] = 0.5 - 0.5 * sin_angle;
    gradient[j] = -0.5 - 0.5 * sin_angle;
}

double mohr_coulomb_principal_function(const double s[3], int i, int j, double sin_angle,
                                       double c_cos_angle, double gradient[3])
{
    if (gradient) mohr_coulomb_principal_gradient(i, j, sin_angle, gradient);
    return 0.5 * (s[i] - s[j]) - 0.5 * (s[i] + s[j]) * sin_angle - c_cos_angle;
}

double mohr_coulomb_mobilised_sin_phi(const double s[3], double a, double fallback)
{
    double denominator = s[0] + s[2] + 2.0 * a;
    if (denominator > ZERO_TOL) return (s[0] - s[2]) / denominator;
    return fallback;
}

double mohr_coulomb_failure_deviator_factor(double sin_phi)
{
    return 2.0 * sin_phi / (1.0 - sin_phi);
}
