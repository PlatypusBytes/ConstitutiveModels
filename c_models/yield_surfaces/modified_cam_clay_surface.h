#pragma once

#include "../globals.h"

/**
 * @file modified_cam_clay_surface.h
 * @brief Elliptic surface of the Modified Cam-Clay model in p-q space, compression positive:
 *
 *   p_eq = p + q^2 / (M^2 p)
 *
 * p_eq is the isotropic stress on the ellipse through (p, q); the ellipse has its top on the line
 * q = M p. It serves as the yield function f = p_eq - p_c of Modified Cam-Clay and as the stress
 * measure and plastic potential of soft soil (creep) models (Stolle, Vermeer & Bonnier, 1999,
 * Eq. 1). The corresponding flow rule is in flow_rules/modified_cam_clay_flow.h.
 *
 * All functions require p > 0.
 */

/**
 * @brief Equivalent pressure p_eq = p + q^2 / (M^2 p) and its partial derivatives.
 *
 * @param[in]  p       Mean stress (compression positive, > 0).
 * @param[in]  q       Deviatoric stress sqrt(3 J2).
 * @param[in]  M       Slope of the line through the tops of the ellipses.
 * @param[out] dpeq_dp d p_eq / d p = 1 - q^2 / (M^2 p^2), or NULL.
 * @param[out] dpeq_dq d p_eq / d q = 2 q / (M^2 p), or NULL.
 * @return equivalent pressure p_eq
 */
double modified_cam_clay_equivalent_pressure(double p, double q, double M, double* dpeq_dp,
                                             double* dpeq_dq);

/**
 * @brief Equivalent pressure of a 3D stress tensor and its gradient with respect to the stress.
 *
 * @param[in]  stress   3D stress tensor (Voigt notation, compression positive, p > 0).
 * @param[in]  M        Slope of the line through the tops of the ellipses.
 * @param[out] gradient d p_eq / d sigma (Voigt notation, shear components doubled so that it
 * contracts with stresses like a strain with engineering shear components), or NULL.
 * @param[out] dpeq_dp  d p_eq / d p, or NULL.
 * @return equivalent pressure p_eq
 */
double modified_cam_clay_equivalent_pressure_3d(const double stress[VOIGTSIZE_3D], double M,
                                                double gradient[VOIGTSIZE_3D], double* dpeq_dp);
