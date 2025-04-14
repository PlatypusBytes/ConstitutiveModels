#include <math.h>
#include "globals.h"
#include "utils.h"
#include "stress_utils.h"


void calculate_stress_invariants_3d(const double stress[VOIGTSIZE_3D],
                                 double* p, double* J, double* theta,
                                 double* j2, double* j3, double s_dev[VOIGTSIZE_3D]) // s_dev is output
{
    // Mean stress pressure p = trace(sigma)/3
    *p = (stress[XX] + stress[YY] + stress[ZZ]) / 3.0;

    // Deviatoric stress tensor s = sigma - p * I
    s_dev[XX] = stress[XX] - *p;
    s_dev[YY] = stress[YY] - *p;
    s_dev[ZZ] = stress[ZZ] - *p;
    s_dev[XY] = stress[XY];
    s_dev[YZ] = stress[YZ];
    s_dev[XZ] = stress[XZ];

    // Second invariant of deviatoric stress J2 = 0.5 * s:s
    // J2 = 0.5 * (sxx^2 + syy^2 + szz^2) + sxy^2 + syz^2 + sxz^2
    *j2 = 0.5 * (s_dev[XX] * s_dev[XX] + s_dev[YY] * s_dev[YY] + s_dev[ZZ] * s_dev[ZZ])
         + (s_dev[XY] * s_dev[XY] + s_dev[YZ] * s_dev[YZ] + s_dev[XZ] * s_dev[XZ]);

    // for convenience, store sqrt(J2) in J
    *J = sqrt((*j2));

    // Third invariant of deviatoric stress J3 = det(s)
    calculate_determinant_voigt_vector_3d(s_dev, j3);

    // Lode angle theta
    if (*J > ZERO_TOL) { // Avoid division by zero if J is very small
        double arg = ( sqrt(27.0) / 2.0) * (*j3) / pow(*j2, 1.5);

        // Clamp argument to asin range [-1, 1] due to potential numerical inaccuracies
        if (arg > 1.0) arg = 1.0;
        if (arg < -1.0) arg = -1.0;

        *theta = asin(arg) / 3.0; // Result is in [-pi/6, pi/6]
    } else {
         *theta = 0.0; // Set to convention if q is zero
    }
    // Ensure theta is within [-pi/6, pi/6] range
    if (*theta < -(PI / 6.0) + ZERO_TOL){
    *theta = -(PI / 6.0) + ZERO_TOL;
    }
    else if (*theta > (PI / 6.0) - ZERO_TOL){
    *theta = (PI / 6.0) - ZERO_TOL;
    }

}

void calculate_stress_invariants_derivatives_3d(const double J, const double s_dev[VOIGTSIZE_3D], const double j2, const double j3,
                                                double dp_dsig[VOIGTSIZE_3D], double dJ_dsig[VOIGTSIZE_3D], double dtheta_dsig[VOIGTSIZE_3D])
{


    // 2. Derivatives of invariants w.r.t. sigma (in Voigt)
    // dp/dsigma = [1/3, 1/3, 1/3, 0, 0, 0]^T
    dp_dsig[0] = 1.0/3.0;
    dp_dsig[1] = 1.0/3.0;
    dp_dsig[2] = 1.0/3.0;
    dp_dsig[3] = 0.0;
    dp_dsig[4] = 0.0;
    dp_dsig[5] = 0.0;

    // calculate dJ/dsigma

    if (J < ZERO_TOL) {
        // Handle case where q is very small or zero
        for (int i = 0; i < VOIGTSIZE_3D; ++i) {
            dJ_dsig[i] = 0.0; // Set to zero or some other value as needed
        }
    }
    else {
            dJ_dsig[0] =s_dev[XX] / (2.0*J);
            dJ_dsig[1] =s_dev[YY] / (2.0*J);
            dJ_dsig[2] =s_dev[ZZ] / (2.0*J);
            dJ_dsig[3] =s_dev[XY] / J;
            dJ_dsig[4] =s_dev[YZ] / J;
            dJ_dsig[5] =s_dev[XZ] / J;
    }


    // dtheta/dsigma = (dtheta/dJ2)*(dJ2/dsig) + (dtheta/dJ3)*(dJ3/dsig)
    // This is the complex part, requires dJ2/dsig, dJ3/dsig, dtheta/dJ2, dtheta/dJ3


    double dJ2_dsig[VOIGTSIZE_3D] = {s_dev[XX], s_dev[YY], s_dev[ZZ], 2.0*s_dev[XY], 2.0*s_dev[YZ], 2.0*s_dev[XZ]};

    double dJ3_dsig[VOIGTSIZE_3D];
    dJ3_dsig[XX] = s_dev[XX]*s_dev[XX] + s_dev[XY]*s_dev[XY] + s_dev[XZ]*s_dev[XZ] - 2.0*j2/3.0; // s_xx^2 + s_xy^2 + s_xz^2 - 2/3 * J2
    dJ3_dsig[YY] = s_dev[YY]*s_dev[YY] + s_dev[XY]*s_dev[XY] + s_dev[YZ]*s_dev[YZ] - 2.0*j2/3.0; // s_yy^2 + s_xy^2 + s_yz^2 - 2/3 * J2
    dJ3_dsig[ZZ] = s_dev[ZZ]*s_dev[ZZ] + s_dev[YZ]*s_dev[YZ] + s_dev[XZ]*s_dev[XZ] - 2.0*j2/3.0; // s_zz^2 + s_yz^2 + s_xz^2 - 2/3 * J2
    dJ3_dsig[XY] = 2.0*((s_dev[XX] + s_dev[YY])*s_dev[XY] + s_dev[XZ]*s_dev[YZ]); // 2*((s_xx+s_yy)*s_xy + s_xz*s_yz)
    dJ3_dsig[YZ] = 2.0*((s_dev[YY] + s_dev[ZZ])*s_dev[YZ] + s_dev[XY]*s_dev[XZ]); // 2*((s_yy+s_zz)*s_yz + s_xy*s_xz)
    dJ3_dsig[XZ] = 2.0*((s_dev[ZZ] + s_dev[XX])*s_dev[XZ] + s_dev[XY]*s_dev[YZ]); // 2*((s_zz+s_xx)*s_xz + s_xy*s_yz)

    double dtheta_dJ2 = 0.0;
    double dtheta_dJ3 = 0.0;

    const double inner_term = 1.0- (27.0 * j3 * j3) / (4.0 * pow(j2, 3.0));

    // prevent sqrt zero
    if (inner_term < ZERO_TOL) {
        dtheta_dJ2 = 0.0;
        dtheta_dJ3= 0.0;
    }
    else {
        dtheta_dJ2 = -pow(3.0, 1.5) * j3 / (4.0* pow(j2,2.5) * sqrt(inner_term));
        dtheta_dJ3 = sqrt(3.0) / sqrt(4.0*pow(j2, 3.0) -27.0 * j3 * j3) ;
    }

    // dtheta/dsigma = dtheta_dJ2 * dJ2_dsig + dtheta_dJ3 * dJ3_dsig
    for(int i=0; i<VOIGTSIZE_3D; ++i) {
        dtheta_dsig[i] = dtheta_dJ2 * dJ2_dsig[i] + dtheta_dJ3 * dJ3_dsig[i];
    }

}