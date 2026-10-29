#pragma once

/**
 * @file cap_hardening.h
 * @brief Volumetric hardening of the pre-consolidation stress of a cap (Schanz et al., 1999,
 * Eqs. 31-35):
 *
 *   dp_c = H deps_v^pc,   H = Ks Kc / (Ks - Kc) = Ks / (Ks/Kc - 1)
 *
 * with Ks the swelling (elastic) and Kc the elasto-plastic bulk modulus in isotropic compression.
 * A stress dependency of H can be added with elastic_laws/power_law_stiffness.h (Eq. 35).
 */

/**
 * @brief Cap hardening modulus H = Ks / (K_ratio - 1) (Eq. 32).
 *
 * @param[in]  Ks      Swelling bulk modulus.
 * @param[in]  K_ratio Ks / Kc, larger than 1.
 * @return hardening modulus H
 */
double cap_hardening_modulus(double Ks, double K_ratio);
