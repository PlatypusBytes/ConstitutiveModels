"""
Python counterpart of c_models/yield_surfaces/hyperbolic_shear_surface.c, mohr_coulomb_surface.c and
elliptic_cap_surface.c.

All functions work on principal stresses s = [s1, s2, s3] (s1 >= s2 >= s3, compression positive) and
return the yield function value together with its derivatives.
"""

import numpy as np

from python_prototypes.stress_utils import PrincipalCorner
from python_prototypes.utils import ZERO_TOL


class YieldFunctions:

    # ------------------------------------------------------------------
    # Hyperbolic shear hardening surface (Eq. 8)
    # ------------------------------------------------------------------
    @staticmethod
    def hyperbolic_shear_yield_function(s, gamma_p, Ei, Eur, k_a, a):
        """
        Shear hardening yield function f13 (Eq. 8), multiplied by (qa - q) to remove the asymptote:
            f = 2/Ei q qa - (2 q / Eur + gamma_p) (qa - q),   q = s1 - s3,  qa = k_a (s3 + a)

        :return: f, df/ds (3), df/dgamma_p
        """
        q = s[0] - s[2]
        qa = k_a * (s[2] + a)
        strain = 2.0 * q / Eur + gamma_p

        df_dq = 2.0 / Ei * qa - 2.0 / Eur * (qa - q) + strain
        df_dqa = 2.0 / Ei * q - strain
        gradient = np.array([df_dq, 0.0, -df_dq + df_dqa * k_a])
        df_dgamma = -(qa - q)

        return 2.0 / Ei * q * qa - strain * (qa - q), gradient, df_dgamma

    # ------------------------------------------------------------------
    # Mohr-Coulomb failure surface and plastic potential (Eqs. 2, 14)
    # ------------------------------------------------------------------
    @staticmethod
    def mohr_coulomb_principal_gradient(i, j, sin_angle):
        """Gradient of f_ij = (s_i - s_j)/2 - (s_i + s_j)/2 sin(angle) - c cos(angle)."""
        gradient = np.zeros(3)
        gradient[i] = 0.5 - 0.5 * sin_angle
        gradient[j] = -0.5 - 0.5 * sin_angle
        return gradient

    @staticmethod
    def mohr_coulomb_principal_function(s, i, j, sin_angle, c_cos_angle):
        """
        :return: f_ij, df_ij/ds (3)
        """
        f = 0.5 * (s[i] - s[j]) - 0.5 * (s[i] + s[j]) * sin_angle - c_cos_angle
        return f, YieldFunctions.mohr_coulomb_principal_gradient(i, j, sin_angle)

    @staticmethod
    def mohr_coulomb_mobilised_sin_phi(s, a, fallback):
        """sin(phi_m) = (s1 - s3) / (s1 + s3 + 2a) (Eq. 12), fallback if the denominator vanishes."""
        denominator = s[0] + s[2] + 2.0 * a
        if denominator > ZERO_TOL:
            return (s[0] - s[2]) / denominator
        return fallback

    @staticmethod
    def mohr_coulomb_failure_deviator_factor(sin_phi):
        """qf = k_f (s3 + a) with k_f = 2 sin(phi) / (1 - sin(phi)) (Eq. 2)."""
        return 2.0 * sin_phi / (1.0 - sin_phi)

    # ------------------------------------------------------------------
    # Elliptic cap (Eqs. 27-29)
    # ------------------------------------------------------------------
    @staticmethod
    def elliptic_cap_shape_factor(sin_phi):
        """alpha = (3 + sin(phi)) / (3 - sin(phi)) (Eq. 29)."""
        return (3.0 + sin_phi) / (3.0 - sin_phi)

    @staticmethod
    def elliptic_cap_weights(alpha, corner):
        """
        Weights w of q~ = w . s (Eq. 28). At a triaxial corner the average of both orderings is used.
        """
        if corner == PrincipalCorner.S2_EQ_S3:
            return np.array([1.0, -0.5, -0.5])
        if corner == PrincipalCorner.S1_EQ_S2:
            return np.array([0.5 * alpha, 0.5 * alpha, -alpha])
        return np.array([1.0, alpha - 1.0, -alpha])

    @staticmethod
    def elliptic_cap_equivalent_deviator(w, s):
        return w @ s

    @staticmethod
    def elliptic_cap_yield_function(s, w, p_c, M, a):
        """
        f = q~^2 / M^2 + (p + a)^2 - (p_c + a)^2 (Eq. 27, shifted by a = c cot(phi)).

        :return: f, df/ds (3), df/dp_c
        """
        M2 = M * M
        q_tilde = w @ s
        p = s.sum() / 3.0

        gradient = 2.0 * q_tilde / M2 * w + 2.0 * (p + a) / 3.0
        df_dpc = -2.0 * (p_c + a)

        return q_tilde * q_tilde / M2 + (p + a) * (p + a) - (p_c + a) * (p_c + a), gradient, df_dpc
