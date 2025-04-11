#include "utils.h"
#define VOIGTSIZE_3D 6


void calculate_determinant_voigt_vector_3d(const double vector[VOIGTSIZE_3D], double* determinant)
{
    // Calculate the determinant of a 3x3 matrix represented as a Voigt vector
    // vector = [sxx, syy, szz, sxy, syz, sxz]
    *determinant = vector[0] * (vector[1] * vector[2] - vector[3] * vector[3])
                - vector[3] * (vector[3] * vector[2] - vector[4] * vector[5])
                + vector[5] * (vector[3] * vector[4] - vector[1] * vector[5]);

//       *determinant = vector[XX] * (vector[YY] * vector[ZZ] - vector[XY] * vector[XY])
//    - vector[XY] * (vector[XY] * vector[ZZ] - vector[YZ] * vector[XZ])
//    + vector[XZ] * (vector[XY] * vector[YZ] - vector[YY] * vector[XZ]);

}

double vector_dot_product(const double* vec1, const double* vec2, int NTENS) {
    // Simple dot product for Voigt vectors (no factors of 2 needed here)
    double dot = 0.0;
    for (int i = 0; i < NTENS; ++i) {
        dot += vec1[i] * vec2[i];
    }
    return dot;
    }

void matrix_vector_multiply(const double* matrix, const double* vector, int NTENS, double* result) {
    // Assumes matrix is row-major 1D array of size NTENS*NTENS
    for (int i = 0; i < NTENS; ++i) {
        result[i] = 0.0;
        for (int j = 0; j < NTENS; ++j) {
            result[i] += matrix[i * NTENS + j] * vector[j];
        }
    }
}
