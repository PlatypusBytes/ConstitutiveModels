
#pragma once

#include "../globals.h"

/**
 * @file matsuoka_nakai_surface.h
 * @brief Header file for Matsuoka-Nakai yield surface calculations.
 *
 * This file contains function declarations for calculating the Matsuoka-Nakai yield function,
 * its first and second derivatives, and the constants used in the calculations. The
 * implementation follows the formulation defined in \cite Lagioia_2016, written tension
 * positive:
 *
 *     f = M p - K + J alpha cos(acos(beta sin(3 theta)) / 3 - gamma pi / 6)
 *
 * with p the mean stress (tension positive), J = sqrt(J2) and sin(3 theta) = 3 sqrt(3) / 2 J3 /
 * J2^(3/2) (theta = -pi/6 in triaxial compression). For the Matsuoka-Nakai criterion
 * I1 I2 / I3 = 9 + 8 tan^2(phi) of the stresses shifted by c cot(phi) \cite Matsuoka_1974:
 *
 *     M     = 2 sqrt(3) sin(phi) / (3 - sin(phi))
 *     K     = 2 sqrt(3) c cos(phi) / (3 - sin(phi))   (= M c cot(phi))
 *     alpha = 2 sqrt(3 + sin^2(phi)) / (3 - sin(phi))
 *     beta  = sin(phi) (9 - sin^2(phi)) / (3 + sin^2(phi))^(3/2)
 *     gamma = 0
 *
 * These expressions are continuous in phi; for phi = 0 the surface is the von Mises cylinder
 * J = 2 c / sqrt(3), which passes through the corners of the Tresca hexagon, as the
 * Matsuoka-Nakai surface passes through the corners of the Mohr-Coulomb pyramid.
 */

/**
 * @brief Structure to hold Matsuoka-Nakai constants.
 *
 * This structure contains the constants used in the Matsuoka-Nakai yield function calculations.
 */
typedef struct
{
    double alpha;  ///< scaling of the deviatoric term
    double beta;   ///< Lode angle dependency, 0 <= |beta| < 1 (0: circular deviatoric section)
    double gamma;  ///< Lode angle shift, 0 for Matsuoka-Nakai
    double K;      ///< intercept, related to the cohesion and the angle (f = M p - K at J = 0)
    double M;      ///< slope with the mean stress, related to the angle
} MatsuokaNakaiConstants;

/**
 * @brief Function to calculate the Matsuoka-Nakai yield function.
 *
 * @param[in]  p Mean stress (pressure), tension positive.
 * @param[in]  theta Lode angle.
 * @param[in]  J Square root of the second deviatoric invariant of the stress tensor.
 * @param[in]  constants struct containing Matsuoka-Nakai constants (alpha, beta, gamma, K, M).
 * @param[out] f Pointer to the output yield function value.
 */
void calculate_yield_function(const double p, const double theta, const double J,
                              const MatsuokaNakaiConstants constants, double* f);

/**
 * @brief Function to calculate the first and second derivatives of the Matsuoka-Nakai yield
 * function (or plastic potential) with respect to the stress.
 *
 * The derivatives are evaluated through the invariants p, J2 and J3, in which the function is
 * smooth for J2 > 0 and |beta| < 1, also in triaxial states where the derivative of the Lode angle
 * is singular. The derivatives are taken with respect to the Voigt stress components, so the
 * shear components of the gradient are engineering (strain-like) components.
 *
 * @param[in]  s_dev Deviatoric stress (Voigt notation).
 * @param[in]  j2 Second invariant of the deviatoric stress.
 * @param[in]  j3 Third invariant of the deviatoric stress.
 * @param[in]  constants struct containing Matsuoka-Nakai constants (alpha, beta, gamma, K, M).
 * @param[out] grad Gradient df/dsigma (Voigt notation).
 * @param[out] hessian Second derivative d2f/dsigma2 (6x6, row-major), may be NULL.
 */
void calculate_yield_derivatives(const double s_dev[VOIGTSIZE_3D], const double j2,
                                 const double j3, const MatsuokaNakaiConstants constants,
                                 double grad[VOIGTSIZE_3D],
                                 double hessian[VOIGTSIZE_3D * VOIGTSIZE_3D]);

/**
 * @brief Support function of the deviatoric section of the surface: the largest value of s : n
 * over the deviatoric stresses s with J alpha cos(acos(beta sin(3 theta)) / 3 - gamma pi / 6) = 1.
 *
 * This is the dual norm of the deviatoric part of the function, which bounds the deviatoric part
 * of its subgradients at the apex (h*(n) <= 1). The section is isotropic, so the maximising s is
 * coaxial with n and the maximisation is over the Lode angle only.
 *
 * @param[in]  j2 Second invariant of the deviatoric tensor n.
 * @param[in]  j3 Third invariant of the deviatoric tensor n.
 * @param[in]  constants struct containing Matsuoka-Nakai constants (alpha, beta, gamma, K, M).
 * @return The support function h*(n) (tensor contraction s : n).
 */
double calculate_deviatoric_support_function(const double j2, const double j3,
                                             const MatsuokaNakaiConstants constants);

/**
 * @brief Function to calculate the Matsuoka-Nakai constants.
 *
 * For a negative angle (contractant plastic potential) the expressions are continued
 * analytically, which gives M < 0 and a mirrored deviatoric section.
 *
 * @param[in]  angle_rad Angle in radians, |angle_rad| < pi / 2.
 * @param[in]  c Cohesion.
 * @return struct containing the calculated constants (alpha, beta, gamma, K, M).
 */
MatsuokaNakaiConstants calculate_matsuoka_nakai_constants(const double angle_rad, const double c);
