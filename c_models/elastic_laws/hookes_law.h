#pragma once

#include "../globals.h"

/**
 * @brief Function to calculate the elastic stiffness matrix for 3D isotropic materials using
 * Hooke's law.
 *
 * @param[in]  E Young's modulus of the material.
 * @param[in]  nu Poisson's ratio of the material.
 * @param[out] elastic_matrix Pointer to the output stiffness matrix (6x6) in row-major order.
 */
void calculate_elastic_stiffness_matrix_3d(double E, double nu,
                                           double elastic_matrix[VOIGTSIZE_3D * VOIGTSIZE_3D]);
