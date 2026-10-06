#include <math.h>
#include <stdlib.h>

#include "globals.h"
#include "utils.h"

double calculate_determinant_voigt_vector_3d(const double vector[VOIGTSIZE_3D])
{
    // Calculate the determinant of a 3x3 matrix represented as a Voigt vector
    // vector = [sxx, syy, szz, sxy, syz, sxz]
    double det = 0.0;
    det = vector[XX] * (vector[YY] * vector[ZZ] - vector[YZ] * vector[YZ]) -
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

int invert_matrix(const double* matrix, const int size, double* inverse)
{
    // Work copy of the matrix, reduced to the identity while the inverse is built up
    double* a = (double*)malloc((size_t)size * (size_t)size * sizeof(double));
    if (a == NULL) return 0;

    double max_entry = 0.0;
    for (int i = 0; i < size * size; ++i)
    {
        a[i] = matrix[i];
        inverse[i] = (i % (size + 1) == 0) ? 1.0 : 0.0;
        if (fabs(a[i]) > max_entry) max_entry = fabs(a[i]);
    }

    int success = (max_entry > 0.0);
    for (int col = 0; success && col < size; ++col)
    {
        // Partial pivoting: row with the largest entry in this column on or below the diagonal
        int piv = col;
        for (int r = col + 1; r < size; ++r)
        {
            if (fabs(a[r * size + col]) > fabs(a[piv * size + col])) piv = r;
        }
        if (fabs(a[piv * size + col]) < 1.0e-14 * max_entry)
        {
            success = 0;
            break;
        }

        if (piv != col)
        {
            for (int k = 0; k < size; ++k)
            {
                double tmp = a[col * size + k];
                a[col * size + k] = a[piv * size + k];
                a[piv * size + k] = tmp;
                tmp = inverse[col * size + k];
                inverse[col * size + k] = inverse[piv * size + k];
                inverse[piv * size + k] = tmp;
            }
        }

        // Normalise the pivot row and eliminate the column from all other rows
        const double pivot = a[col * size + col];
        for (int k = 0; k < size; ++k)
        {
            a[col * size + k] /= pivot;
            inverse[col * size + k] /= pivot;
        }
        for (int r = 0; r < size; ++r)
        {
            if (r == col) continue;
            const double factor = a[r * size + col];
            if (factor == 0.0) continue;
            for (int k = 0; k < size; ++k)
            {
                a[r * size + k] -= factor * a[col * size + k];
                inverse[r * size + k] -= factor * inverse[col * size + k];
            }
        }
    }

    free(a);
    return success;
}
