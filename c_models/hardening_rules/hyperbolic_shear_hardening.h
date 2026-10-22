#pragma once

/**
 * @file hyperbolic_shear_hardening.h
 * @brief Hyperbolic shear hardening law of the Hardening Soil model (Schanz et al., 1999,
 * Eqs. 1, 7-9): the plastic shear strain gamma_p = eps1_p - eps2_p - eps3_p follows from the
 * hyperbolic stress-strain relation of a drained triaxial test,
 *
 *   gamma_p = 2/Ei * q / (1 - q/qa) - 2 q / Eur,   q < qa
 *
 * The corresponding yield function is in yield_surfaces/hyperbolic_shear_surface.h.
 */

/**
 * @brief Plastic shear strain on the hyperbola for a deviator stress q < qa (Eq. 8 with f = 0).
 *
 * @param[in]  q   Deviator stress s1 - s3.
 * @param[in]  qa  Asymptotic deviator stress qf / Rf.
 * @param[in]  Ei  Initial stiffness of the hyperbola.
 * @param[in]  Eur Unloading / reloading stiffness.
 * @return plastic shear strain gamma_p
 */
double hyperbolic_plastic_shear_strain(double q, double qa, double Ei, double Eur);

/**
 * @brief Initial stiffness Ei = 2 E50 / (2 - Rf) of the hyperbola (Eq. 1), such that E50 is the
 * secant stiffness at q = qf / 2 (Sec. 2.1).
 *
 * @param[in]  E50 Secant stiffness at 50% of the strength.
 * @param[in]  Rf  Failure ratio qf / qa.
 * @return initial stiffness Ei
 */
double hyperbolic_initial_stiffness(double E50, double Rf);
