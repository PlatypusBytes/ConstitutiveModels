#pragma once

/**
 * @file isotache_creep.h
 * @brief Volumetric creep law of the soft soil creep model (Stolle, Vermeer & Bonnier, 1999,
 * Eqs. 2, 3 and 5), compression positive. The volumetric creep strain rate depends on the ratio
 * of the equivalent pressure p_eq to the pre-consolidation pressure p_p (isotaches):
 *
 *   d eps_c / dt = C / tau (p_eq / p_p)^(B / C)                                        (Eq. 2)
 *
 * and the accumulated creep strain hardens the pre-consolidation pressure,
 * p_p = p_p0 exp(eps_c / B) (Eq. 3, hardening_rules/exponential_volumetric_hardening.h).
 *
 * Parameters: B = lambda* - kappa*, C = mu* (modified creep index) and the reference time tau,
 * which fixes the time scale together with the initial over-consolidation ratio.
 */

/**
 * @brief Volumetric creep strain rate C / tau (p_eq / p_p)^(B / C) (Eq. 2).
 *
 * @param[in]  p_eq Equivalent pressure.
 * @param[in]  p_p  Pre-consolidation pressure.
 * @param[in]  B    Hardening index lambda* - kappa*.
 * @param[in]  C    Modified creep index mu*.
 * @param[in]  tau  Reference time.
 * @return creep strain rate (0 for a non-positive p_eq or p_p)
 */
double isotache_creep_rate(double p_eq, double p_p, double B, double C, double tau);

/**
 * @brief Volumetric creep strain increment over a time increment dt at constant p_eq, from the
 * integration of Eq. 2 with the hardening of Eq. 3 (Eq. 5):
 *
 *   deps_c = C ln(1 + dt / tau (p_eq / p_p)^(B / C))
 *
 * with p_p the pre-consolidation pressure at the start of the increment. Since the hardening is
 * integrated exactly, consecutive increments at constant p_eq add up to one increment over the
 * total time. The power is evaluated in logarithmic form, so large ratios p_eq / p_p do not
 * overflow.
 *
 * @param[in]  p_eq              Equivalent pressure.
 * @param[in]  p_p               Pre-consolidation pressure at the start of the increment.
 * @param[in]  B                 Hardening index lambda* - kappa*.
 * @param[in]  C                 Modified creep index mu*.
 * @param[in]  tau               Reference time.
 * @param[in]  dt                Time increment.
 * @param[out] d_increment_d_peq d deps_c / d p_eq, or NULL.
 * @return creep strain increment (0 for dt <= 0 or a non-positive p_eq or p_p)
 */
double isotache_creep_increment(double p_eq, double p_p, double B, double C, double tau,
                                double dt, double* d_increment_d_peq);
