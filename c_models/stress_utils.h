#ifndef STRESS_UTILS_H
#define STRESS_UTILS_H

#include "globals.h"

/**
 * @brief calculates stress invariants for 3D stress tensor
 * @param stress 3D stress tensor (6 components in Voigt notation)
 * @param p Mean stress (pressure)
 * @param J square root of second deviatoric invariant of the stress tensor
 * @param theta Lode angle
 * @param j2 Second deviatoric invariant of the deviatoric stress tensor
 * @param j3 Third deviatoric invariant of the deviatoric stress tensor
 * @param s_dev Deviatoric stress tensor (output)
 */
void calculate_stress_invariants_3d(const double* stress, double* p, double* J, double* theta,
                                    double* j2, double* j3, double* s_dev);  // s_dev is output

/**
 * @brief calculates the derivatives of the stress invariants with respect to the stress tensor
 * @param J square root of second deviatoric invariant of the stress tensor
 * @param s_dev Deviatoric stress tensor (Voigt notation)
 * @param j2 Second deviatoric invariant of the deviatoric stress tensor
 * @param j3 Third deviatoric invariant of the deviatoric stress tensor
 * @param dp_dsig Derivative of mean stress with respect to stress tensor (output)
 * @param dJ_dsig Derivative of J with respect to stress tensor (output)
 * @param dtheta_dsig Derivative of theta with respect to stress tensor (output)
 */
void calculate_stress_invariants_derivatives_3d(const double J, const double* s_dev,
                                                const double j2, const double j3,
                                                double* dp_dsig, double* dJ_dsig,
                                                double* dtheta_dsig);


void calculate_principle_stresses_3d(const double stress[VOIGTSIZE_3D], double principle_stresses[3])

#endif  // STRESS_UTILS_H
