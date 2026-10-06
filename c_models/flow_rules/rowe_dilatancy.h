#pragma once

/**
 * @file rowe_dilatancy.h
 * @brief Rowe's stress-dilatancy theory as used in the Hardening Soil model (Schanz et al., 1999,
 * Eqs. 11-13):
 *
 *   sin(psi_m) = (sin(phi_m) - sin(phi_cv)) / (1 - sin(phi_m) sin(phi_cv))
 *
 * The mobilised dilatancy is contractant (negative) for phi_m < phi_cv. The critical state
 * friction angle follows from the same relation at failure (phi_m = phi, psi_m = psi).
 */

/**
 * @brief Mobilised dilatancy sin(psi_m) (Eq. 11).
 *
 * @param[in]  sin_phi_m  Mobilised friction sin(phi_m).
 * @param[in]  sin_phi_cv Critical state friction sin(phi_cv).
 * @return sin(psi_m)
 */
double rowe_mobilised_sin_psi(double sin_phi_m, double sin_phi_cv);

/**
 * @brief Critical state friction sin(phi_cv) = (sin(phi) - sin(psi)) / (1 - sin(phi) sin(psi))
 * (Eq. 13).
 *
 * @param[in]  sin_phi Friction sin(phi) at failure.
 * @param[in]  sin_psi Dilatancy sin(psi) at failure.
 * @return sin(phi_cv)
 */
double rowe_critical_state_sin_phi(double sin_phi, double sin_psi);
