#pragma once

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

/**
 * @brief Inverts a 3x3 matrix.
 *
 * @param[in]  matrix  Input matrix (3x3) in row-major order.
 * @param[out] inverse Inverse of the matrix (3x3) in row-major order.
 * @return 1 on success, 0 if the matrix is singular (|det| < SMALL_VALUE).
 */
int invert_matrix_3x3(const double matrix[9], double inverse[9]);

/**
 * @brief Solves the linear system matrix * solution = rhs by Gaussian elimination with partial
 * pivoting.
 *
 * The system is the leading n x n block of a row-major matrix with row_stride columns, so small
 * systems of varying size can be stored in a fixed size array. matrix and rhs are overwritten.
 *
 * @param[in]     n          Number of equations.
 * @param[in]     row_stride Number of columns of the storage of matrix (row_stride >= n).
 * @param[in,out] matrix     System matrix, overwritten by the elimination.
 * @param[in,out] rhs        Right hand side (length n), overwritten by the elimination.
 * @param[out]    solution   Solution vector (length n).
 * @return 1 on success, 0 if the matrix is singular (pivot < SMALL_VALUE).
 */
int solve_linear_system(int n, int row_stride, double* matrix, double* rhs, double* solution);
