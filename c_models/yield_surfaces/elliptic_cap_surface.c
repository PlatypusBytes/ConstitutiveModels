#include "elliptic_cap_surface.h"

double elliptic_cap_shape_factor(double sin_phi)
{
    return (3.0 + sin_phi) / (3.0 - sin_phi);
}

void elliptic_cap_weights(double alpha, PrincipalCorner corner, double w[3])
{
    w[0] = 1.0;
    w[1] = alpha - 1.0;
    w[2] = -alpha;
    if (corner == PRINCIPAL_CORNER_S2_EQ_S3)
    {
        w[1] = -0.5;
        w[2] = -0.5;
    }
    else if (corner == PRINCIPAL_CORNER_S1_EQ_S2)
    {
        w[0] = 0.5 * alpha;
        w[1] = 0.5 * alpha;
    }
}

double elliptic_cap_equivalent_deviator(const double w[3], const double s[3])
{
    return w[0] * s[0] + w[1] * s[1] + w[2] * s[2];
}

double elliptic_cap_yield_function(const double s[3], const double w[3], double p_c, double M,
                                   double a, double gradient[3], double* df_dpc)
{
    double M2 = M * M;
    double q_tilde = elliptic_cap_equivalent_deviator(w, s);
    double p = (s[0] + s[1] + s[2]) / 3.0;

    if (gradient)
    {
        double dp_term = 2.0 * (p + a) / 3.0;
        for (int r = 0; r < 3; ++r) gradient[r] = 2.0 * q_tilde / M2 * w[r] + dp_term;
    }
    if (df_dpc) *df_dpc = -2.0 * (p_c + a);

    return q_tilde * q_tilde / M2 + (p + a) * (p + a) - (p_c + a) * (p_c + a);
}
