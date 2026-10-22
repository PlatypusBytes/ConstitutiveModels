#pragma once

/**
 * @file modified_cam_clay_flow.h
 * @brief Flow rule associated with the Modified Cam-Clay ellipse
 * (yield_surfaces/modified_cam_clay_surface.h), compression positive. The inelastic strain
 * increment is normal to the ellipse through the current stress, scaled by its volumetric part
 * (Stolle, Vermeer & Bonnier, 1999, Eq. 7):
 *
 *   deps = deps_v / (d p_eq/d p) d p_eq / d sigma
 *
 * so that the conjugate deviatoric strain increment is
 *
 *   deps_q = deps_v (d p_eq/d q) / (d p_eq/d p) = deps_v 2 p q / (M^2 p^2 - q^2)
 *
 * In an implicit (backward Euler) update with elastic shear modulus G this gives a radial return
 * in the deviatoric plane, q = q_trial - 3 G deps_q (Stolle et al., Eq. 12; Borja & Lee, 1990).
 */

/**
 * @brief Deviatoric stress of the implicit update
 *
 *   q = q_trial - 3 G deps_v 2 p q / (M^2 p^2 - q^2)                                (Eq. 12)
 *
 * for a given mean stress p and a compactive volumetric inelastic strain increment deps_v >= 0
 * (wet side, e.g. a creep strain). Multiplied by (M^2 p^2 - q^2) this is a cubic in q, which has
 * exactly one root in [0, min(q_trial, M p)]; that root is returned. For deps_v <= 0 no
 * inelastic strain develops and q = q_trial.
 *
 * @param[in]  q_trial    Elastic trial deviatoric stress (>= 0).
 * @param[in]  p          Mean stress (compression positive, > 0).
 * @param[in]  G          Elastic shear modulus.
 * @param[in]  deps_v     Volumetric inelastic strain increment.
 * @param[in]  M          Slope of the line through the tops of the ellipses.
 * @param[out] dq_dp      d q / d p at fixed q_trial and deps_v, or NULL.
 * @param[out] dq_ddeps_v d q / d deps_v at fixed q_trial and p, or NULL.
 * @return deviatoric stress q
 */
double modified_cam_clay_deviatoric_return(double q_trial, double p, double G, double deps_v,
                                           double M, double* dq_dp, double* dq_ddeps_v);
