#include "globals.h"
#include "utils.h"

double calculate_determinant_voigt_vector_3d(const double vector[VOIGTSIZE_3D])
{
    // Calculate the determinant of a 3x3 matrix represented as a Voigt vector
    // vector = [sxx, syy, szz, sxy, syz, sxz]
    double det = 0.0;
    det = vector[XX] * (vector[YY] * vector[ZZ] - vector[XY] * vector[XY]) -
          vector[XY] * (vector[XY] * vector[ZZ] - vector[YZ] * vector[XZ]) +
          vector[XZ] * (vector[XY] * vector[YZ] - vector[YY] * vector[XZ]);

    return det;
}

double vector_dot_product(const double* vector_1, const double* vector_2, int length_vector)
{
    // Simple dot product for vectors
    double dot = 0.0;
    for (int i = 0; i < length_vector; ++i)
    {
        dot += vector_1[i] * vector_2[i];
    }
    return dot;
}

void matrix_vector_multiply(const double* matrix, const double* vector, const int length_vector,
                            double* result)
{
    // Assumes matrix is row-major 1D array of size length_vector*length_vector
    for (int i = 0; i < length_vector; ++i)
    {
        result[i] = 0.0;
        for (int j = 0; j < length_vector; ++j)
        {
            result[i] += matrix[i * length_vector + j] * vector[j];
        }
    }
}

void copy_array(const double* source, const int length_vector, double* destination)
{
    // Copy array of length length_vector
    for (int i = 0; i < length_vector; ++i)
    {
        destination[i] = source[i];
    }
}

void add_vectors(const double* vector_1, const double* vector_2, const int length_vector,
                 double* result)
{
    // Add two vectors of length length_vector
    for (int i = 0; i < length_vector; ++i)
    {
        result[i] = vector_1[i] + vector_2[i];
    }
}

void vector_scalar_multiply(const double* vector, const double scalar, const int length_vector,
                            double* result)
{
    // Multiply vector by scalar
    for (int i = 0; i < length_vector; ++i)
    {
        result[i] = vector[i] * scalar;
    }
}

void vector_outer_product(const double* vector_1, const double* vector_2, const int length_vector,
                          double* result)
{
    for (int i = 0; i < length_vector; ++i)
    {
        for (int j = 0; j < length_vector; ++j)
        {
            result[i * length_vector + j] = vector_1[i] * vector_2[j];
        }
    }
}
