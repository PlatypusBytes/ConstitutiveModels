#include <math.h>

#include "logarithmic_elasticity.h"

double logarithmic_bulk_modulus(double p, double kappa, double p_min)
{
    return fmax(p, p_min) / kappa;
}

double logarithmic_mean_stress(double p0, double deps_v_e, double kappa, double p_min)
{
    return fmax(p0, p_min) * exp(deps_v_e / kappa);
}
