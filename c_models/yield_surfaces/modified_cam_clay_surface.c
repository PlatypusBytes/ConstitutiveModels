#include <math.h>

#include "../globals.h"
#include "../stress_utils.h"
#include "modified_cam_clay_surface.h"

double modified_cam_clay_equivalent_pressure(double p, double q, double M, double* dpeq_dp,
                                             double* dpeq_dq)
{
    double M2p = M * M * p;
    if (dpeq_dp) *dpeq_dp = 1.0 - q * q / (M2p * p);
    if (dpeq_dq) *dpeq_dq = 2.0 * q / M2p;
    return p + q * q / M2p;
}

double modified_cam_clay_equivalent_pressure_3d(const double stress[VOIGTSIZE_3D], double M,
                                                double gradient[VOIGTSIZE_3D], double* dpeq_dp)
{
    double p, J, theta, j2, j3, s_dev[VOIGTSIZE_3D];
    calculate_stress_invariants_3d(stress, &p, &J, &theta, &j2, &j3, s_dev);

    // q = sqrt(3 J2)
    double dpeq_dp_value, dpeq_dq;
    double p_eq =
        modified_cam_clay_equivalent_pressure(p, sqrt(3.0) * J, M, &dpeq_dp_value, &dpeq_dq);

    if (gradient)
    {
        double dp_dsig[VOIGTSIZE_3D], dJ_dsig[VOIGTSIZE_3D], dtheta_dsig[VOIGTSIZE_3D];
        calculate_stress_invariants_derivatives_3d(J, s_dev, j2, j3, dp_dsig, dJ_dsig,
                                                   dtheta_dsig);
        for (int i = 0; i < VOIGTSIZE_3D; ++i)
            gradient[i] = dpeq_dp_value * dp_dsig[i] + dpeq_dq * sqrt(3.0) * dJ_dsig[i];
    }
    if (dpeq_dp) *dpeq_dp = dpeq_dp_value;
    return p_eq;
}
