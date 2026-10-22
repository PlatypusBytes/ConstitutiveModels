#pragma once

/**
 * @file logarithmic_elasticity.h
 * @brief Pressure dependent (hypo-)elasticity with a linear relation between the elastic
 * volumetric strain and the logarithm of the mean stress, as in Cam-clay type and soft soil
 * models (e.g. Stolle, Vermeer & Bonnier, 1999, Eq. 6):
 *
 *   d eps_v^e = kappa d p / p,   i.e.   K = p / kappa
 *
 * with kappa the modified swelling index (kappa*, called A in Stolle et al.). Compression positive.
 * The shear modulus usually follows from K with a constant Poisson's ratio, see
 * calculate_shear_modulus_from_bulk_modulus in hookes_law.h.
 *
 * The mean stress is limited to p_min from below, so the stiffness stays positive at zero stress.
 */

/**
 * @brief Tangent bulk modulus K = max(p, p_min) / kappa.
 *
 * @param[in]  p     Mean stress (compression positive).
 * @param[in]  kappa Modified swelling index.
 * @param[in]  p_min Lower limit of the mean stress.
 * @return bulk modulus
 */
double logarithmic_bulk_modulus(double p, double kappa, double p_min);

/**
 * @brief Mean stress after an elastic volumetric strain increment, from the exact integration of
 * d p = (p / kappa) d eps_v^e:
 *
 *   p = max(p0, p_min) exp(deps_v^e / kappa)
 *
 * The result is always positive; its derivative with respect to deps_v^e is p / kappa.
 *
 * @param[in]  p0       Mean stress at the start of the increment (compression positive).
 * @param[in]  deps_v_e Elastic volumetric strain increment (compression positive).
 * @param[in]  kappa    Modified swelling index.
 * @param[in]  p_min    Lower limit of the mean stress at the start of the increment.
 * @return mean stress at the end of the increment
 */
double logarithmic_mean_stress(double p0, double deps_v_e, double kappa, double p_min);
