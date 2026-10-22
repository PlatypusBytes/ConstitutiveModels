#include <math.h>

#include "isotache_creep.h"

double isotache_creep_rate(double p_eq, double p_p, double B, double C, double tau)
{
    if (p_eq <= 0.0 || p_p <= 0.0) return 0.0;
    return C / tau * pow(p_eq / p_p, B / C);
}

double isotache_creep_increment(double p_eq, double p_p, double B, double C, double tau,
                                double dt, double* d_increment_d_peq)
{
    if (d_increment_d_peq) *d_increment_d_peq = 0.0;
    if (dt <= 0.0 || p_eq <= 0.0 || p_p <= 0.0) return 0.0;

    // x = ln(dt / tau (p_eq / p_p)^(B / C)); deps_c = C ln(1 + e^x), with
    // d deps_c / d p_eq = B / p_eq e^x / (1 + e^x)
    double x = log(dt / tau) + B / C * log(p_eq / p_p);
    double increment, fraction;
    if (x > 0.0)
    {
        double e = exp(-x);
        increment = C * (x + log1p(e));
        fraction = 1.0 / (1.0 + e);
    }
    else
    {
        double e = exp(x);
        increment = C * log1p(e);
        fraction = e / (1.0 + e);
    }

    if (d_increment_d_peq) *d_increment_d_peq = B * fraction / p_eq;
    return increment;
}
