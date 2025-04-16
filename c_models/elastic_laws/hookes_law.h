#ifndef HOOKES_LAW_H
#define HOOKES_LAW_H

/**
 * @brief Function to calculate the elastic stiffness matrix for 3D isotropic materials using
 * Hooke's law.
 *
 * @param E Young's modulus of the material.
 * @param nu Poisson's ratio of the material.
 * @param elastic_matrix Pointer to the output stiffness matrix (6x6) in row-major order.
 */
void calculate_elastic_stiffness_matrix_3d(double E, double nu, double* elastic_matrix);

#endif  // HOOKES_LAW_H
