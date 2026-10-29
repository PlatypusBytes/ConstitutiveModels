#include "hyperbolic_shear_hardening.h"

double hyperbolic_plastic_shear_strain(double q, double qa, double Ei, double Eur)
{
    return 2.0 / Ei * q / (1.0 - q / qa) - 2.0 * q / Eur;
}

double hyperbolic_initial_stiffness(double E50, double Rf)
{
    return 2.0 * E50 / (2.0 - Rf);
}
