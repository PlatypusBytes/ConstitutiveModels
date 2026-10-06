#pragma once

/**
 * @file hyperbolic_shear_surface.h
 * @brief Shear hardening yield function of the Hardening Soil model (Schanz et al., 1999, Eq. 8)
 * in principal stress space, compression positive, s1 >= s2 >= s3:
 *
 *   f = 2/Ei * q / (1 - q/qa) - 2 q / Eur - gamma_p,   q = s1 - s3,   qa = k_a (s3 + a)
 *
 * The function is evaluated multiplied by (qa - q), which removes the pole at the asymptote and
 * gives the quadratic form of Sec. 3 (Eq. 25). For q < qa the sign is unchanged, and for q >= qa
 * the scaled function is strictly positive, so a stress beyond the asymptote is always detected
 * as yielding. The hardening law gamma_p(q) is in hardening_rules/hyperbolic_shear_hardening.h.
 */

/**
 * @brief Scaled hyperbolic shear yield function (qa - q) f.
 *
 * @param[in]  s         Principal stresses (compression positive), s1 >= s2 >= s3.
 * @param[in]  gamma_p   Plastic shear strain (hardening parameter).
 * @param[in]  Ei        Initial stiffness of the hyperbola.
 * @param[in]  Eur       Unloading / reloading stiffness.
 * @param[in]  k_a       Asymptote factor, qa = k_a (s3 + a), i.e. k_f / Rf.
 * @param[in]  a         Shift of the stress origin, c cot(phi).
 * @param[out] gradient  d((qa - q) f)/ds (3 components), or NULL.
 * @param[out] df_dgamma d((qa - q) f)/d gamma_p, or NULL.
 * @return value of the scaled yield function
 */
double hyperbolic_shear_yield_function(const double s[3], double gamma_p, double Ei, double Eur,
                                       double k_a, double a, double gradient[3], double* df_dgamma);
