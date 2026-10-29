#pragma once

#include "globals.h"

/**
 * @brief calculates the volumetric strain eps_v = trace(strain)
 * @param[in]  strain 3D strain tensor (6 components in Voigt notation)
 * @return volumetric strain
 */
double calculate_volumetric_strain(const double strain[VOIGTSIZE_3D]);

/**
 * @brief rotates principal values of a strain-like tensor back to the global frame:
 * strain = Q * diag(principal_strain) * Q^T, with engineering shear components (gamma = 2 eps).
 *
 * This is the strain counterpart of calculate_stress_from_principal_system and also applies to
 * stress gradients (e.g. yield function gradients df/dsigma) that are contracted with stresses.
 *
 * @param[in]  principal_strain principal values
 * @param[in]  Q principal directions stored as columns (as returned by calculate_principal_system)
 * @param[out] strain 3D strain tensor (6 components in Voigt notation, engineering shear strains)
 */
void calculate_strain_from_principal_system(const double principal_strain[3], double Q[3][3],
                                            double strain[VOIGTSIZE_3D]);

/**
 * @brief calculates the void ratio from the volumetric strain for finite volume changes:
 * e = (1 + e0) exp(-eps_v) - 1
 * @param[in]  e0 void ratio at eps_v = 0
 * @param[in]  eps_v volumetric strain, compression positive
 * @return void ratio
 */
double calculate_void_ratio(double e0, double eps_v);
