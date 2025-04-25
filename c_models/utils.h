#ifndef UTILS_H
#define UTILS_H

#include "globals.h"

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
double calculate_determinant_voigt_vector_3d(const double vector[VOIGTSIZE_3D]);

/**
 * @brief Computes the dot product of two vectors.
 *
 * Calculates the scalar (dot) product of two vectors of length length_vector.
 *
 * @param[in] vector_1 Pointer to the first input vector.
 * @param[in] vector_2 Pointer to the second input vector.
 * @param[in] length_vector Number of components in each vector.
 * @return The dot product of vector_1 and vector_2.
 */
double vector_dot_product(const double* vector_1, const double* vector_2, int length_vector);

/**
 * @brief Multiplies a square matrix by a vector.
 *
 * Performs the multiplication of an length_vector×length_vector matrix with a vector of length
 * length_vector. The result is stored in the provided result array.
 *
 * @param[in]  matrix Pointer to the input matrix stored in row-major order (length
 * length_vector×length_vector).
 * @param[in]  vector Pointer to the input vector (length length_vector).
 * @param[in]  length_vector  The dimension of the matrix and vector.
 * @param[out] result Pointer to the output vector (length length_vector) where the result is
 * stored.
 */
void matrix_vector_multiply(const double* matrix, const double* vector, int length_vector,
                            double* result);

/**
 * @brief Copies an array from source to destination.
 *
 * Copies an array of length length_vector from source to destination.
 *
 * @param[in]  source Pointer to the source array (length length_vector).
 * @param[in]  length_vector The number of elements in the array.
 * @param[out] destination Pointer to the destination array (length length_vector).
 */
void copy_array(const double* source, const int length_vector, double* destination);

/**
 * @brief Adds two vectors and stores the result in a third vector.
 *
 * Adds two vectors of length length_vector and stores the result in the provided result array.
 *
 * @param[in]  vector_1 Pointer to the first input vector (length length_vector).
 * @param[in]  vector_2 Pointer to the second input vector (length length_vector).
 * @param[in]  length_vector The number of elements in each vector.
 * @param[out] result Pointer to the output vector (length length_vector) where the result is
 * stored.
 */
void add_vectors(const double* vector_1, const double* vector_2, const int length_vector,
                 double* result);

/**
 * @brief Multiplies a vector by a scalar and stores the result in a third vector.
 *
 * Multiplies each component of the input vector by the scalar and stores the result in the provided
 * result array.
 *
 * @param[in]  vector Pointer to the input vector (length length_vector).
 * @param[in]  scalar The scalar value to multiply with.
 * @param[in]  length_vector The number of elements in the vector.
 * @param[out] result Pointer to the output vector (length length_vector) where the result is
 * stored.
 */
void vector_scalar_multiply(const double* vector, const double scalar, const int length_vector,
                            double* result);

/**
 * @brief Computes the outer product of two vectors and stores the result in a matrix.
 *
 * Computes the outer product of two vectors of length length_vector and stores the result in a
 * matrix. The result is stored in a 1D array in row-major order.
 *
 * @param[in]  vector_1 Pointer to the first input vector (length length_vector).
 * @param[in]  vector_2 Pointer to the second input vector (length length_vector).
 * @param[in]  length_vector The number of elements in each vector.
 * @param[out] result Pointer to the output matrix (length length_vector×length length_vector) where
 * the result is stored.
 */
void vector_outer_product(const double* vector_1, const double* vector_2, const int length_vector,
                          double* result);

#endif  // UTILS_H
