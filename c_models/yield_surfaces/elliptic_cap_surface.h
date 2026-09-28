#pragma once

#include "../stress_utils.h"

/**
 * @file elliptic_cap_surface.h
 * @brief Elliptic cap yield function of the Hardening Soil model (Schanz et al., 1999,
 * Eqs. 27-29) in principal stress space, compression positive, s1 >= s2 >= s3:
 *
 *   fc = q~^2 / M^2 + (p + a)^2 - (p_c + a)^2,   q~ = w . s,   w = [1, alpha - 1, -alpha]
 *
 * With associated flow, fc is also the plastic potential (Eq. 30).
 */

/**
 * @brief Cap shape factor alpha = (3 + sin(phi)) / (3 - sin(phi)) of q~ (Eq. 29), such that q~
 * equals s1 - s3 in triaxial compression and alpha (s1 - s3) in triaxial extension.
 */
double elliptic_cap_shape_factor(double sin_phi);

/**
 * @brief Weights w of q~ = w . s (Eq. 28).
 *
 * q~ is not smooth at the triaxial corners; there the average over both orderings of the equal
 * principal stresses is used, which gives the values of the paper at the corners:
 * q~ = s1 - s3 (s2 = s3) and q~ = alpha (s1 - s3) (s1 = s2).
 *
 * @param[in]  alpha  Cap shape factor.
 * @param[in]  corner Triaxial corner at which q~ is averaged, or PRINCIPAL_CORNER_NONE.
 * @param[out] w      Weights (3 components).
 */
void elliptic_cap_weights(double alpha, PrincipalCorner corner, double w[3]);

/**
 * @brief Equivalent deviator stress q~ = w . s of the cap.
 */
double elliptic_cap_equivalent_deviator(const double w[3], const double s[3]);

/**
 * @brief Cap yield function fc.
 *
 * @param[in]  s        Principal stresses (compression positive).
 * @param[in]  w        Weights of q~ (see elliptic_cap_weights).
 * @param[in]  p_c      Pre-consolidation stress (cap size).
 * @param[in]  M        Cap aspect ratio.
 * @param[in]  a        Shift of the stress origin, c cot(phi).
 * @param[out] gradient dfc/ds (3 components), or NULL. Also the flow direction (associated).
 * @param[out] df_dpc   dfc/dp_c, or NULL.
 * @return value of fc
 */
double elliptic_cap_yield_function(const double s[3], const double w[3], double p_c, double M,
                                   double a, double gradient[3], double* df_dpc);
