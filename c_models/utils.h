#ifndef UTILS_H
#define UTILS_H

/**
 * @brief Calculates the determinant of a 3D tensor represented in Voigt notation.
 *
 * This function computes the determinant of a 3x3 symmetric tensor stored as a 6-component
 * vector in Voigt notation:
 * [xx, yy, zz, xy, yz, xz]
 *
 * @param[in]  vector      Pointer to an array of 6 doubles representing the tensor in Voigt
 * notation.
 * @return The determinant of the corresponding 3×3 tensor.
 */
double calculate_determinant_voigt_vector_3d(const double* vector);

/**
 * @brief Computes the dot product of two vectors.
 *
 * Calculates the scalar (dot) product of two vectors of length NTENS.
 *
 * @param[in] vec1 Pointer to the first input vector.
 * @param[in] vec2 Pointer to the second input vector.
 * @param[in] NTENS Number of components in each vector.
 * @return The dot product of vec1 and vec2.
 */
double vector_dot_product(const double* vec1, const double* vec2, int NTENS);

/**
 * @brief Multiplies a square matrix by a vector.
 *
 * Performs the multiplication of an NTENS×NTENS matrix with a vector of length NTENS.
 * The result is stored in the provided result array.
 *
 * @param[in]  matrix Pointer to the input matrix stored in row-major order (length NTENS×NTENS).
 * @param[in]  vector Pointer to the input vector (length NTENS).
 * @param[in]  NTENS  The dimension of the matrix and vector.
 * @param[out] result Pointer to the output vector (length NTENS) where the result is stored.
 */
void matrix_vector_multiply(const double* matrix, const double* vector, int NTENS, double* result);

#endif  // UTILS_H
