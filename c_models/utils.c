#include <math.h>

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

int invert_matrix_3x3(const double matrix[9], double inverse[9])
{
    const double* A = matrix;
    double c00 = A[4] * A[8] - A[5] * A[7];
    double c01 = A[5] * A[6] - A[3] * A[8];
    double c02 = A[3] * A[7] - A[4] * A[6];
    double det = A[0] * c00 + A[1] * c01 + A[2] * c02;
    if (fabs(det) < SMALL_VALUE) return 0;

    double inv = 1.0 / det;
    inverse[0] = c00 * inv;
    inverse[1] = (A[2] * A[7] - A[1] * A[8]) * inv;
    inverse[2] = (A[1] * A[5] - A[2] * A[4]) * inv;
    inverse[3] = c01 * inv;
    inverse[4] = (A[0] * A[8] - A[2] * A[6]) * inv;
    inverse[5] = (A[2] * A[3] - A[0] * A[5]) * inv;
    inverse[6] = c02 * inv;
    inverse[7] = (A[1] * A[6] - A[0] * A[7]) * inv;
    inverse[8] = (A[0] * A[4] - A[1] * A[3]) * inv;
    return 1;
}

int solve_linear_system(int n, int row_stride, double* matrix, double* rhs, double* solution)
{
    double* J = matrix;
    for (int col = 0; col < n; ++col)
    {
        // partial pivoting
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (fabs(J[r * row_stride + col]) > fabs(J[piv * row_stride + col])) piv = r;
        if (fabs(J[piv * row_stride + col]) < SMALL_VALUE) return 0;

        if (piv != col)
        {
            for (int k = 0; k < n; ++k)
            {
                double tmp = J[col * row_stride + k];
                J[col * row_stride + k] = J[piv * row_stride + k];
                J[piv * row_stride + k] = tmp;
            }
            double tmp = rhs[col];
            rhs[col] = rhs[piv];
            rhs[piv] = tmp;
        }

        // forward elimination
        for (int r = col + 1; r < n; ++r)
        {
            double factor = J[r * row_stride + col] / J[col * row_stride + col];
            for (int k = col; k < n; ++k) J[r * row_stride + k] -= factor * J[col * row_stride + k];
            rhs[r] -= factor * rhs[col];
        }
    }

    // back substitution
    for (int r = n - 1; r >= 0; --r)
    {
        double sum = rhs[r];
        for (int k = r + 1; k < n; ++k) sum -= J[r * row_stride + k] * solution[k];
        solution[r] = sum / J[r * row_stride + r];
    }
    return 1;
}
