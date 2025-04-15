#ifndef UTILS_H
#define UTILS_H

#ifdef __cplusplus
extern "C"
{
#endif

    void calculate_determinant_voigt_vector_3d(const double* vector, double* determinant);

    double vector_dot_product(const double* vec1, const double* vec2, int NTENS);

    void matrix_vector_multiply(const double* matrix, const double* vector, int NTENS,
                                double* result);

#ifdef __cplusplus
}
#endif

#endif  // UTILS_H
