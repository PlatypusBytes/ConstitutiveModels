#include "globals.h"


// small strains only
double calculate_volumetric_strain_3d(double strains[VOIGTSIZE_3D])
{
    // Calculate the volumetric strain from the 3D strain tensor
    return strains[XX] + strains[YY] + strains[ZZ];
}

double calculate_deviatoric_strain_triaxial_3d(double sigma_1, double sigma_3)
{
    // Calculate the deviatoric strain from the 3D principle strain tensor
    return 2/3 * (sigma_1 - sigma_3);
}