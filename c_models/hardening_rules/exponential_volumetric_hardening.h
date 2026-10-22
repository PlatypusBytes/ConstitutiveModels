#pragma once

/**
 * @file exponential_volumetric_hardening.h
 * @brief Exponential volumetric hardening of an isotropic pre-consolidation pressure, as in
 * Cam-clay type and soft soil models (Stolle, Vermeer & Bonnier, 1999, Eq. 3), compression
 * positive:
 *
 *   p_p = p_p0 exp(eps_v / B)
 *
 * with eps_v the irreversible (plastic or creep) volumetric strain and B = lambda* - kappa* the
 * difference of the modified compression and swelling indices.
 */

/**
 * @brief Pre-consolidation pressure after an irreversible volumetric strain increment,
 * p_p = p_p0 exp(deps_v / B) (Eq. 3). Its derivative with respect to deps_v is p_p / B.
 *
 * @param[in]  p_p0   Pre-consolidation pressure at the start of the increment.
 * @param[in]  deps_v Irreversible volumetric strain increment (compression positive).
 * @param[in]  B      Hardening index lambda* - kappa*.
 * @return pre-consolidation pressure at the end of the increment
 */
double exponential_volumetric_hardening(double p_p0, double deps_v, double B);
