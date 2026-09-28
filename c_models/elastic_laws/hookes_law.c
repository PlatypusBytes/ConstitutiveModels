#include "../globals.h"
#include "hookes_law.h"

void calculate_elastic_stiffness_matrix_3d(double E, double nu,
                                           double elastic_matrix[VOIGTSIZE_3D * VOIGTSIZE_3D])
{
    double G = E / (2.0 * (1.0 + nu));                         // Shear modulus
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));  // Lame's first param

    double factor = lambda + 2.0 * G;

    // Initialize to zero
    for (int i = 0; i < VOIGTSIZE_3D * VOIGTSIZE_3D; ++i)
    {
        elastic_matrix[i] = 0.0;
    }

    // Populate using row-major indexing: elastic_matrix[row * VOIGTSIZE_3D + col]
    // Normal Stresses
    elastic_matrix[0 * VOIGTSIZE_3D + 0] = factor;  // C_1111 (row 0, col 0)
    elastic_matrix[0 * VOIGTSIZE_3D + 1] = lambda;  // C_1122 (row 0, col 1)
    elastic_matrix[0 * VOIGTSIZE_3D + 2] = lambda;  // C_1133 (row 0, col 2)

    elastic_matrix[1 * VOIGTSIZE_3D + 0] = lambda;  // C_2211 (row 1, col 0)
    elastic_matrix[1 * VOIGTSIZE_3D + 1] = factor;  // C_2222 (row 1, col 1)
    elastic_matrix[1 * VOIGTSIZE_3D + 2] = lambda;  // C_2233 (row 1, col 2)

    elastic_matrix[2 * VOIGTSIZE_3D + 0] = lambda;  // C_3311 (row 2, col 0)
    elastic_matrix[2 * VOIGTSIZE_3D + 1] = lambda;  // C_3322 (row 2, col 1)
    elastic_matrix[2 * VOIGTSIZE_3D + 2] = factor;  // C_3333 (row 2, col 2)

    // Shear Stresses (using engineering shear strain convention gamma = 2*epsilon_shear)
    // The stiffness term C_ij = G for i != j in engineering strain notation
    elastic_matrix[3 * VOIGTSIZE_3D + 3] = G;  // C_1212 (row 3, col 3)
    elastic_matrix[4 * VOIGTSIZE_3D + 4] = G;  // C_2323 (row 4, col 4)
    elastic_matrix[5 * VOIGTSIZE_3D + 5] = G;  // C_1313 (row 5, col 5)
}

void calculate_elastic_stiffness_matrix_2d_interface(double E, double nu,
                                           double elastic_matrix[VOIGTSIZE_2D_INTERFACE * VOIGTSIZE_2D_INTERFACE])
{
    double M =  E * (1 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)); // P wave modulus
    double G = E / (2.0 * (1.0 + nu));                         // Shear modulus

    // Initialize to zero
    for (int i = 0; i < VOIGTSIZE_2D_INTERFACE * VOIGTSIZE_2D_INTERFACE; ++i)
    {
        elastic_matrix[i] = 0.0;
    }

    // Populate using row-major indexing: elastic_matrix[row * VOIGTSIZE_2D_INTERFACE + col]
    // Normal Stresses
    elastic_matrix[0 * VOIGTSIZE_2D_INTERFACE + 0] = M;  // (row 0, col 0)
    elastic_matrix[1 * VOIGTSIZE_2D_INTERFACE + 1] = G;  // (row 1, col 1)

}

void calculate_elastic_stiffness_matrix_3d_interface(double E, double nu,
                                           double elastic_matrix[VOIGTSIZE_3D_INTERFACE * VOIGTSIZE_3D_INTERFACE])
{
    double M =  E * (1 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)); // P wave modulus
    double G = E / (2.0 * (1.0 + nu));                         // Shear modulus

    // Initialize to zero
    for (int i = 0; i < VOIGTSIZE_3D_INTERFACE * VOIGTSIZE_3D_INTERFACE; ++i)
    {
        elastic_matrix[i] = 0.0;
    }

    // Populate using row-major indexing: elastic_matrix[row * VOIGTSIZE_3D_INTERFACE + col]
    // Normal Stresses
    elastic_matrix[0 * VOIGTSIZE_3D_INTERFACE + 0] = M;  // (row 0, col 0)
    elastic_matrix[1 * VOIGTSIZE_3D_INTERFACE + 1] = G;  // (row 1, col 1)
    elastic_matrix[2 * VOIGTSIZE_3D_INTERFACE + 2] = G;  // (row 2, col 2)

}