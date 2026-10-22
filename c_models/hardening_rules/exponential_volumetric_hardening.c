#include <math.h>

#include "exponential_volumetric_hardening.h"

double exponential_volumetric_hardening(double p_p0, double deps_v, double B)
{
    return p_p0 * exp(deps_v / B);
}
