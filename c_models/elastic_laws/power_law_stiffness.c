#include <math.h>

#include "power_law_stiffness.h"

double calculate_power_law_stiffness_factor(double sigma, double a, double p_ref, double m,
                                            double min_ratio)
{
    double ratio = (sigma + a) / (p_ref + a);
    if (ratio < min_ratio) ratio = min_ratio;
    return pow(ratio, m);
}
