"""
Python counterpart of c_models/elastic_laws/hookes_law.c and power_law_stiffness.c.
"""

import numpy as np


class ElasticLaws:

    @staticmethod
    def shear_modulus(E, nu):
        return E / (2.0 * (1.0 + nu))

    @staticmethod
    def bulk_modulus(E, nu):
        return E / (3.0 * (1.0 - 2.0 * nu))

    @staticmethod
    def elastic_stiffness_matrix_3d(E, nu):
        """6x6 isotropic stiffness matrix for engineering shear strains."""
        G = ElasticLaws.shear_modulus(E, nu)
        lame_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        D = np.zeros((6, 6))
        D[:3, :3] = lame_lambda
        D[0, 0] = D[1, 1] = D[2, 2] = lame_lambda + 2.0 * G
        D[3, 3] = D[4, 4] = D[5, 5] = G
        return D

    @staticmethod
    def elastic_stiffness_matrix_principal(E, nu):
        """3x3 isotropic stiffness matrix in principal stress space."""
        G = ElasticLaws.shear_modulus(E, nu)
        lame_lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        return np.full((3, 3), lame_lambda) + 2.0 * G * np.eye(3)

    @staticmethod
    def power_law_stiffness_factor(sigma, a, p_ref, m, min_ratio):
        """
        ((sigma + a) / (p_ref + a))^m, with the stress ratio limited to min_ratio (Eqs. 3-5).
        """
        ratio = (sigma + a) / (p_ref + a)
        if ratio < min_ratio:
            ratio = min_ratio
        return ratio ** m
