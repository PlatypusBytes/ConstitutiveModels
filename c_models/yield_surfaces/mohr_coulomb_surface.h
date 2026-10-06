#pragma once

/**
 * @file mohr_coulomb_surface.h
 * @brief Mohr-Coulomb yield function and plastic potential in principal stress space.
 *
 * Compression positive, with principal stresses ordered s1 >= s2 >= s3. For a pair of principal
 * stresses (i, j):
 *
 *   f_ij = (s_i - s_j) / 2 - (s_i + s_j) / 2 sin(phi) - c cos(phi)
 *
 * The same function with the dilation angle psi in place of phi is the plastic potential g_ij;
 * the cohesion term does not affect its gradient.
 */

/**
 * @brief Mohr-Coulomb function f_ij of the principal stress pair (i, j).
 *
 * @param[in]  s           Principal stresses (compression positive).
 * @param[in]  i           Index of the larger principal stress of the pair.
 * @param[in]  j           Index of the smaller principal stress of the pair.
 * @param[in]  sin_angle   sin(phi) for the yield function, sin(psi) for the plastic potential.
 * @param[in]  c_cos_angle c cos(phi); irrelevant for the plastic potential gradient.
 * @param[out] gradient    df_ij/ds (3 components), or NULL.
 * @return value of f_ij
 */
double mohr_coulomb_principal_function(const double s[3], int i, int j, double sin_angle,
                                       double c_cos_angle, double gradient[3]);

/**
 * @brief Gradient of f_ij (or g_ij) with respect to the principal stresses, which is constant.
 * Its plastic shear strain increment is d(eps_i^p) - d(eps_j^p) = dLambda_ij.
 *
 * @param[in]  i         Index of the larger principal stress of the pair.
 * @param[in]  j         Index of the smaller principal stress of the pair.
 * @param[in]  sin_angle sin(phi) for the yield function, sin(psi) for the plastic potential.
 * @param[out] gradient  df_ij/ds (3 components).
 */
void mohr_coulomb_principal_gradient(int i, int j, double sin_angle, double gradient[3]);

/**
 * @brief Mobilised friction sin(phi_m) = (s1 - s3) / (s1 + s3 + 2a) of ordered principal stresses
 * (compression positive).
 *
 * @param[in]  s        Principal stresses, s1 >= s2 >= s3.
 * @param[in]  a        Shift of the stress origin, c cot(phi).
 * @param[in]  fallback Value returned when s1 + s3 + 2a is not positive.
 * @return sin(phi_m)
 */
double mohr_coulomb_mobilised_sin_phi(const double s[3], double a, double fallback);

/**
 * @brief Factor k_f of the deviator stress at failure, q_f = s1 - s3 = k_f (s3 + a), with
 * k_f = 2 sin(phi) / (1 - sin(phi)).
 */
double mohr_coulomb_failure_deviator_factor(double sin_phi);
