#pragma once

#include "../globals.h"

/**
 * @file mohr_coulomb_return_mapping.h
 * @brief Return mapping for perfectly plastic Mohr-Coulomb plasticity with a non-associated
 * Mohr-Coulomb plastic potential, in principal stress space (compression positive, principal
 * stresses ordered s1 >= s2 >= s3). The yield functions and potentials are in
 * yield_surfaces/mohr_coulomb_surface.h.
 *
 * The yield functions are linear in the principal stresses, so for a constant elastic stiffness
 * the return to a plane is a single correction (e.g. Stolle, Vermeer & Bonnier, 1999, Eq. 15):
 *
 *   s = s_trial - dLambda D dg/ds,   dLambda = f(s_trial) / (df/ds . D dg/ds)
 *
 * When the principal stress order is lost, the return is made to the edge of the pyramid in
 * triaxial compression (s2 = s3, planes f13 and f12) or extension (s1 = s2, planes f13 and f23),
 * and in tension beyond these to the apex s1 = s2 = s3 = -c cot(phi).
 */

/**
 * @brief Region of the Mohr-Coulomb pyramid to which the stress is returned.
 */
typedef enum
{
    MC_RETURN_ELASTIC,        ///< inside the yield surface, no return
    MC_RETURN_PLANE,          ///< main plane f13
    MC_RETURN_EDGE_S2_EQ_S3,  ///< edge of f13 and f12, triaxial compression
    MC_RETURN_EDGE_S1_EQ_S2,  ///< edge of f13 and f23, triaxial extension
    MC_RETURN_APEX            ///< apex s1 = s2 = s3 = -c cot(phi)
} MohrCoulombReturnType;

/**
 * @brief Result of the Mohr-Coulomb return mapping.
 */
typedef struct
{
    MohrCoulombReturnType type;  ///< region of the returned stress
    int n_active;                ///< number of active planes (1 on a plane, 2 on an edge, else 0)
    int planes[2][2];            ///< principal stress pairs (i, j) of the active planes f_ij
    double dlambda[2];           ///< plastic multipliers of the active planes
} MohrCoulombReturn;

/**
 * @brief Returns the principal trial stress to the Mohr-Coulomb yield surface.
 *
 * @param[in]  s_trial Principal trial stresses (compression positive), s1 >= s2 >= s3.
 * @param[in]  D       Elastic stiffness matrix acting on the principal stresses (3x3, row-major).
 * @param[in]  sin_phi sin of the friction angle.
 * @param[in]  cos_phi cos of the friction angle.
 * @param[in]  sin_psi sin of the dilatancy angle (plastic potential).
 * @param[in]  c       Cohesion.
 * @param[out] s       Returned principal stresses (equal to s_trial for an elastic state).
 * @param[out] ret     Region and active planes of the returned stress.
 * @return 1 on success, 0 if no admissible return is found (only for sin_phi = 0 in tension).
 */
int mohr_coulomb_return_mapping(const double s_trial[3], const double D[9], double sin_phi,
                                double cos_phi, double sin_psi, double c, double s[3],
                                MohrCoulombReturn* ret);

/**
 * @brief Applies the plastic part of the Mohr-Coulomb tangent to a tangent stiffness matrix:
 *
 *   T <- (I - De N (A^T De N)^-1 A^T) T
 *
 * with A and N the yield function gradients and flow directions of the active planes, rotated to
 * the global frame with the principal directions Q. For T = De this gives the elasto-plastic
 * (continuum) tangent of perfect plasticity; for another tangent T of a preceding stress update it
 * gives the tangent of that update followed by the plastic return. At the apex the stress does
 * not change and T is set to zero.
 *
 * @param[in]     ret     Result of mohr_coulomb_return_mapping.
 * @param[in]     Q       Principal directions stored as columns (see calculate_principal_system).
 * @param[in]     De      Elastic stiffness matrix used in the return (6x6, row-major).
 * @param[in]     sin_phi sin of the friction angle.
 * @param[in]     sin_psi sin of the dilatancy angle.
 * @param[in,out] tangent Tangent stiffness matrix (6x6, row-major).
 */
void mohr_coulomb_project_tangent(const MohrCoulombReturn* ret, double Q[3][3],
                                  const double De[VOIGTSIZE_3D * VOIGTSIZE_3D], double sin_phi,
                                  double sin_psi, double tangent[VOIGTSIZE_3D * VOIGTSIZE_3D]);
