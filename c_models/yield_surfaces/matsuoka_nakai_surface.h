
#ifndef MATSUOKA_NAKAI_SURFACE_H
#define MATSUOKA_NAKAI_SURFACE_H

/**
 * @file matsuoka_nakai_surface.h
 * @brief Header file for Matsuoka-Nakai yield surface calculations.
 *
 * This file contains function declarations for calculating the Matsuoka-Nakai yield function,
 * its gradient, and the constants used in the calculations. The implementation follows the
 * formulation defined in \cite Lagioia_2016.
 */

/**
 * @brief Function to calculate the Matsuoka-Nakai yield function.
 *
 * @param p Mean stress (pressure).
 * @param theta Lode angle.
 * @param J Square root of the second deviatoric invariant of the stress tensor.
 * @param alpha Matsuoka-Nakai constant.
 * @param beta Matsuoka-Nakai constant.
 * @param gamma Matsuoka-Nakai constant.
 * @param K Matsuoka-Nakai constant, related to the angle and cohesion.
 * @param M Matsuoka-Nakai constant, related to the angle.
 * @param f Pointer to the output yield function value.
 */
void calculate_yield_function(const double p, const double theta, const double J,
                              const double alpha, const double beta, const double gamma,
                              const double K, const double M, double* f);

/**
 * @brief Function to calculate the gradient of the Matsuoka-Nakai yield function.
 *
 * @param theta Lode angle.
 * @param J Square root of the second deviatoric invariant of the stress tensor.
 * @param constants Array containing Matsuoka-Nakai constants (angle, M, alpha, beta, gamma).
 * @param dp_dsig Derivative of mean stress with respect to stress tensor.
 * @param dJ_dsig Derivative of J with respect to stress tensor .
 * @param dtheta_dsig Derivative of theta with respect to stress tensor.
 * @param grad Pointer to the output gradient vector.
 */
void calculate_yield_gradient(const double theta, const double J, const double* constants,
                              const double* dp_dsig, const double* dJ_dsig,
                              const double* dtheta_dsig, double* grad);

/**
 * @brief Function to calculate the Matsuoka-Nakai constants.
 *
 * @param angle_rad Angle in radians.
 * @param c Cohesion.
 * @param alpha Pointer to output alpha constant (output).
 * @param beta Pointer to output beta constant (output).
 * @param gamma Pointer to output gamma constant (output).
 * @param K Pointer to output K constant, related to the angle and cohesion (output).
 * @param M Pointer to output M constant, related to the angle (output).
 */
void calculate_matsuoka_nakai_constants(const double angle_rad, const double c, double* alpha,
                                        double* beta, double* gamma, double* K, double* M);

#endif  // MATSUOKA_NAKAI_SURFACE_H
