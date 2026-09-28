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



/**
 * @brief Triaxial corners of isotropic yield surfaces in principal stress space, for principal
 * stresses ordered s1 >= s2 >= s3.
 */
typedef enum
{
    PRINCIPAL_CORNER_NONE,
    PRINCIPAL_CORNER_S2_EQ_S3, ///< s2 = s3 (triaxial compression when compression is positive)
    PRINCIPAL_CORNER_S1_EQ_S2  ///< s1 = s2 (triaxial extension when compression is positive)
} PrincipalCorner;

/**
 * @brief calculates the mean stress p = trace(stress) / 3
 * @param[in]  stress 3D stress tensor (6 components in Voigt notation)
 * @return mean stress
 */
double calculate_mean_stress(const double stress[VOIGTSIZE_3D]);

/**
 * @brief calculates the deviatoric stress s = stress - p I
 * @param[in]  stress 3D stress tensor (6 components in Voigt notation)
 * @param[in]  p mean stress of the stress tensor
 * @param[out] s_dev deviatoric stress tensor (6 components in Voigt notation)
 */
void calculate_deviatoric_stress(const double stress[VOIGTSIZE_3D], double p, double s_dev[VOIGTSIZE_3D]);

/**
 * @brief calculates the von Mises equivalent stress q = sqrt(3 J2), including the shear components
 * @param[in]  stress 3D stress tensor (6 components in Voigt notation)
 * @return von Mises equivalent stress
 */
double calculate_von_mises_stress(const double stress[VOIGTSIZE_3D]);

/**
 * @brief calculates the algebraically smallest principal stress
 * @param[in]  stress 3D stress tensor (6 components in Voigt notation)
 * @return smallest principal stress (the minor principal stress when compression is positive)
 */
double calculate_min_principal_stress(const double stress[VOIGTSIZE_3D]);
