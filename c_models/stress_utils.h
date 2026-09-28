#pragma once

#include "globals.h"


#define EIG_TOL 1.0e-10
#define JACOBI_MAX_ITER 50

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

/**
 * @brief calculates the principal stresses and principal directions (Jacobi eigensolver)
 * @param[in]  stress 3D stress tensor (6 components in Voigt notation)
 * @param[out] principal_stress principal stresses, sorted descending (s1 >= s2 >= s3)
 * @param[out] Q principal directions stored as columns, Q[:][i] belongs to principal_stress[i]
 */
void calculate_principal_system(const double stress[VOIGTSIZE_3D], double principal_stress[3], double Q[3][3]);

/**
 * @brief rotates principal stresses back to the global frame: stress = Q * diag(principal_stress) * Q^T
 * @param[in]  principal_stress principal stresses
 * @param[in]  Q principal directions stored as columns (as returned by calculate_principal_system)
 * @param[out] stress 3D stress tensor (6 components in Voigt notation)
 */
void calculate_stress_from_principal_system(const double principal_stress[3], double Q[3][3],
                                            double stress[VOIGTSIZE_3D]);



