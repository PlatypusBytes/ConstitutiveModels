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


/**
 * @brief Function to calculate the elastic stiffness matrix for 2D interface materials using
 * Hooke's law. For 2D interface elements, the stiffness matrix is 2x2. With one component in the
 * normal direction and one in the shear direction.
 *
 * @param[in]  E Young's modulus of the material.
 * @param[in]  nu Poisson's ratio of the material.
 * @param[out] elastic_matrix Pointer to the output stiffness matrix (2x2) in row-major order.
 */
void calculate_elastic_stiffness_matrix_2d_interface(double E, double nu,
                                           double elastic_matrix[VOIGTSIZE_2D_INTERFACE * VOIGTSIZE_2D_INTERFACE]);

/**
 * @brief Function to calculate the elastic stiffness matrix for 3D interface materials using
 * Hooke's law. For 3D interface elements, the stiffness matrix is 3x3. With one component in the
 * normal direction and two in the shear direction.
 *
 * @param[in]  E Young's modulus of the material.
 * @param[in]  nu Poisson's ratio of the material.
 * @param[out] elastic_matrix Pointer to the output stiffness matrix (3x3) in row-major order.
 */
void calculate_elastic_stiffness_matrix_2d_interface(double E, double nu,
                                           double elastic_matrix[VOIGTSIZE_3D_INTERFACE * VOIGTSIZE_3D_INTERFACE]);