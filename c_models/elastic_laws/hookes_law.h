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
void calculate_elastic_stiffness_matrix_3d_interface(double E, double nu,
                                           double elastic_matrix[VOIGTSIZE_3D_INTERFACE * VOIGTSIZE_3D_INTERFACE]);
/**
 * @brief Function to calculate the isotropic elastic stiffness matrix acting on the principal
 * stresses and strains, d sigma_i = D_ij d eps_j.
 *
 * @param[in]  E Young's modulus of the material.
 * @param[in]  nu Poisson's ratio of the material.
 * @param[out] elastic_matrix Pointer to the output stiffness matrix (3x3) in row-major order.
 */
void calculate_elastic_stiffness_matrix_principal(double E, double nu, double elastic_matrix[9]);

/**
 * @brief Shear modulus G = E / (2 (1 + nu)).
 */
double calculate_shear_modulus(double E, double nu);

/**
 * @brief Bulk modulus K = E / (3 (1 - 2 nu)).
 */
double calculate_bulk_modulus(double E, double nu);

/**
 * @brief Shear modulus G = 3 K (1 - 2 nu) / (2 (1 + nu)) for a given bulk modulus and Poisson's
 * ratio, e.g. for a pressure dependent bulk modulus.
 */
double calculate_shear_modulus_from_bulk_modulus(double K, double nu);

/**
 * @brief Isotropic elastic stiffness matrix for 3D materials from the bulk and shear modulus,
 * which allows K and G to be evaluated independently (e.g. a pressure dependent K). With K = 0 the
 * matrix maps a strain on the deviatoric stress 2 G e.
 *
 * @param[in]  K Bulk modulus.
 * @param[in]  G Shear modulus.
 * @param[out] elastic_matrix Pointer to the output stiffness matrix (6x6) in row-major order,
 * engineering shear strains.
 */
void calculate_elastic_stiffness_matrix_3d_bulk_shear(
    double K, double G, double elastic_matrix[VOIGTSIZE_3D * VOIGTSIZE_3D]);

/**
 * @brief Isotropic elastic stiffness matrix acting on the principal stresses and strains, from the
 * bulk and shear modulus, d sigma_i = D_ij d eps_j.
 *
 * @param[in]  K Bulk modulus.
 * @param[in]  G Shear modulus.
 * @param[out] elastic_matrix Pointer to the output stiffness matrix (3x3) in row-major order.
 */
void calculate_elastic_stiffness_matrix_principal_bulk_shear(double K, double G,
                                                             double elastic_matrix[9]);
