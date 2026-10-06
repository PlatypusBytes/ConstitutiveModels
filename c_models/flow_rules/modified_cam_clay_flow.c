#include <math.h>

#include "../globals.h"
#include "modified_cam_clay_flow.h"

#define MCC_RETURN_MAX_ITER 200
#define MCC_RETURN_TOL 1.0e-15

double modified_cam_clay_deviatoric_return(double q_trial, double p, double G, double deps_v,
                                           double M, double* dq_dp, double* dq_ddeps_v)
{
    double M2p2 = M * M * p * p;
    // no inelastic strain for deps_v <= 0; the derivatives are then the ones for deps_v -> 0+
    double c = 6.0 * G * fmax(deps_v, 0.0) * p;
    double q = q_trial;

    if (q_trial <= 0.0)
    {
        q = 0.0;
    }
    else if (deps_v > 0.0)
    {
        // h(q) = (q - q_trial) (q^2 - M^2 p^2) - 6 G deps_v p q, with h(0) >= 0 and h(q_max) <= 0:
        // safeguarded Newton iteration within the bracket [lo, hi]
        double q_max = fmin(q_trial, M * p);
        double lo = 0.0, hi = q_max;

        // initial guess from the linearisation for small q
        q = q_trial / (1.0 + c / M2p2);
        if (!(q > lo && q < hi)) q = 0.5 * (lo + hi);

        for (int iter = 0; iter < MCC_RETURN_MAX_ITER; ++iter)
        {
            double h = (q - q_trial) * (q * q - M2p2) - c * q;
            if (h == 0.0) break;
            if (h > 0.0)
                lo = q;
            else
                hi = q;

            double dh_dq = 3.0 * q * q - 2.0 * q_trial * q - M2p2 - c;
            double q_new = (dh_dq < 0.0) ? q - h / dh_dq : 0.5 * (lo + hi);
            if (!(q_new > lo && q_new < hi)) q_new = 0.5 * (lo + hi);

            double change = fabs(q_new - q);
            q = q_new;
            if (change <= MCC_RETURN_TOL * q_max) break;
        }
    }

    // implicit derivatives of h(q, p, deps_v) = 0
    if (dq_dp || dq_ddeps_v)
    {
        double dh_dq = 3.0 * q * q - 2.0 * q_trial * q - M2p2 - c;
        double dh_dp = -2.0 * M * M * p * (q - q_trial) - c / p * q;
        double dh_ddeps_v = -6.0 * G * p * q;
        int regular = (dh_dq < -SMALL_VALUE * (M2p2 + fabs(c))) ? 1 : 0;
        if (dq_dp) *dq_dp = regular ? -dh_dp / dh_dq : 0.0;
        if (dq_ddeps_v) *dq_ddeps_v = regular ? -dh_ddeps_v / dh_dq : 0.0;
    }
    return q;
}
