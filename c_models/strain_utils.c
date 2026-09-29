#include <math.h>

#include "globals.h"
#include "strain_utils.h"

double calculate_volumetric_strain(const double strain[VOIGTSIZE_3D])
{
    return strain[XX] + strain[YY] + strain[ZZ];
}

void calculate_strain_from_principal_system(const double principal_strain[3], double Q[3][3],
                                            double strain[VOIGTSIZE_3D])
{
    for (int i = 0; i < VOIGTSIZE_3D; ++i) strain[i] = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        strain[XX] += principal_strain[k] * Q[0][k] * Q[0][k];
        strain[YY] += principal_strain[k] * Q[1][k] * Q[1][k];
        strain[ZZ] += principal_strain[k] * Q[2][k] * Q[2][k];
        // engineering shear strains
        strain[XY] += 2.0 * principal_strain[k] * Q[0][k] * Q[1][k];
        strain[YZ] += 2.0 * principal_strain[k] * Q[1][k] * Q[2][k];
        strain[XZ] += 2.0 * principal_strain[k] * Q[0][k] * Q[2][k];
    }
}

double calculate_void_ratio(double e0, double eps_v)
{
    return (1.0 + e0) * exp(-eps_v) - 1.0;
}
