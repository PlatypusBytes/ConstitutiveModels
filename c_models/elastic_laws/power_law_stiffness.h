#pragma once

/**
 * @file power_law_stiffness.h
 * @brief Stress dependent stiffness of the power law type (Ohde / Janbu), as used in the
 * Hardening Soil model (Schanz et al., 1999, Eqs. 3-5):
 *
 *   E = E_ref * ((sigma + a) / (p_ref + a))^m,   a = c cot(phi)
 *
 * with sigma a stress measure (e.g. the minor principal stress, compression positive).
 */

/**
 * @brief Stress dependency factor ((sigma + a) / (p_ref + a))^m.
 *
 * The stress ratio is limited to min_ratio from below, so the stiffness stays positive at zero
 * confinement and in tension.
 *
 * @param[in]  sigma     Stress measure (compression positive).
 * @param[in]  a         Shift of the stress origin, c cot(phi).
 * @param[in]  p_ref     Reference stress.
 * @param[in]  m         Power of the stress dependency.
 * @param[in]  min_ratio Lower limit of (sigma + a) / (p_ref + a).
 * @return stiffness factor, multiplies the reference stiffness.
 */
double calculate_power_law_stiffness_factor(double sigma, double a, double p_ref, double m,
                                            double min_ratio);
