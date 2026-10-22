#include "hyperbolic_shear_surface.h"

double hyperbolic_shear_yield_function(const double s[3], double gamma_p, double Ei, double Eur,
                                       double k_a, double a, double gradient[3], double* df_dgamma)
{
    double q = s[0] - s[2];
    double qa = k_a * (s[2] + a);
    double strain = 2.0 * q / Eur + gamma_p;

    if (gradient)
    {
        double df_dq = 2.0 / Ei * qa - 2.0 / Eur * (qa - q) + strain;
        double df_dqa = 2.0 / Ei * q - strain;
        gradient[0] = df_dq;
        gradient[1] = 0.0;
        gradient[2] = -df_dq + df_dqa * k_a;
    }
    if (df_dgamma) *df_dgamma = -(qa - q);

    return 2.0 / Ei * q * qa - strain * (qa - q);
}
