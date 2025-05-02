
#pragma once

#include "../globals.h"

/**
 * @file matsuoka_nakai_surface.h
 * @brief Header file for Matsuoka-Nakai yield surface calculations.
 *
 * This file contains function declarations for calculating the Matsuoka-Nakai yield function,
 * its gradient, and the constants used in the calculations. The implementation follows the
 * formulation defined in \cite Lagioia_2016.
 */

/**
 * @brief Structure to hold Matsuoka-Nakai constants.
 *
 * This structure contains the constants used in the Matsuoka-Nakai yield function calculations.
 */
typedef struct
{
    double alpha;      ///< Matsuoka-Nakai constant
    double beta;       ///< Matsuoka-Nakai constant
    double gamma;      ///< Matsuoka-Nakai constant
    double K;          ///< Matsuoka-Nakai constant, related to the angle and cohesion
    double M;          ///< Matsuoka-Nakai constant, related to the angle
} MatsuokaNakaiConstants;

/**
 * @brief Function to calculate the Matsuoka-Nakai yield function.
 *
 * @param[in]  p Mean stress (pressure).
 * @param[in]  theta Lode angle.
 * @param[in]  J Square root of the second deviatoric invariant of the stress tensor.
 * @param[in]  constants struct containing Matsuoka-Nakai constants (angle, alpha, beta, gamma K,
 * M).
 * @param[out] f Pointer to the output yield function value.
 */
void calculate_yield_function(const double p, const double theta, const double J,
                              const MatsuokaNakaiConstants constants, double* f);

/**
 * @brief Function to calculate the gradient of the Matsuoka-Nakai yield function.
 *
 * @param[in]  theta Lode angle.
 * @param[in]  J Square root of the second deviatoric invariant of the stress tensor.
 * @param[in]  constants struct containing Matsuoka-Nakai constants (angle, alpha, beta, gamma K,
 * M).
 * @param[in]  dp_dsig Derivative of mean stress with respect to stress tensor.
 * @param[in]  dJ_dsig Derivative of J with respect to stress tensor .
 * @param[in]  dtheta_dsig Derivative of theta with respect to stress tensor.
 * @param[out] grad Pointer to the output gradient vector.
 */
void calculate_yield_gradient(const double theta, const double J,
                              const MatsuokaNakaiConstants constants,
                              const double dp_dsig[VOIGTSIZE_3D],
                              const double dJ_dsig[VOIGTSIZE_3D],
                              const double dtheta_dsig[VOIGTSIZE_3D], double grad[VOIGTSIZE_3D]);

/**
 * @brief Function to calculate the Matsuoka-Nakai constants.
 *
 * @param[in]  angle_rad Angle in radians.
 * @param[in]  c Cohesion.
 * @return struct containing the calculated constants (angle_rad, alpha, beta, gamma, K, M).
 */
MatsuokaNakaiConstants calculate_matsuoka_nakai_constants(const double angle_rad, const double c);
