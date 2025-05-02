#pragma once

#include "globals.h"

/**
 * @brief calculates stress invariants for 3D stress tensor
 * @param[in]  stress 3D stress tensor (6 components in Voigt notation)
 * @param[out] p Mean stress (pressure)
 * @param[out] J square root of second deviatoric invariant of the stress tensor
 * @param[out] theta Lode angle
 * @param[out] j2 Second deviatoric invariant of the deviatoric stress tensor
 * @param[out] j3 Third deviatoric invariant of the deviatoric stress tensor
 * @param[out] s_dev Deviatoric stress tensor
 */
void calculate_stress_invariants_3d(const double stress[VOIGTSIZE_3D], double* p, double* J, double* theta,
                                    double* j2, double* j3, double s_dev[VOIGTSIZE_3D]);

/**
 * @brief calculates the derivatives of the stress invariants with respect to the stress tensor
 * @param[in]  J square root of second deviatoric invariant of the stress tensor
 * @param[in]  s_dev Deviatoric stress tensor (Voigt notation)
 * @param[in]  j2 Second deviatoric invariant of the deviatoric stress tensor
 * @param[in]  j3 Third deviatoric invariant of the deviatoric stress tensor
 * @param[out] dp_dsig Derivative of mean stress with respect to stress tensor (output)
 * @param[out] dJ_dsig Derivative of J with respect to stress tensor (output)
 * @param[out] dtheta_dsig Derivative of theta with respect to stress tensor (output)
 */
void calculate_stress_invariants_derivatives_3d(const double J, const double s_dev[VOIGTSIZE_3D],
                                                const double j2, const double j3,
                                                double dp_dsig[VOIGTSIZE_3D], double dJ_dsig[VOIGTSIZE_3D],
                                                double dtheta_dsig[VOIGTSIZE_3D]);

void calculate_principle_stresses_3d(const double stress[VOIGTSIZE_3D], double principle_stresses[3])



