"""
Python counterpart of c_models/stress_utils.c.

Voigt ordering: [xx, yy, zz, xy, yz, xz] (see globals.h).
"""

from enum import IntEnum

import numpy as np

from python_prototypes.utils import XX, YY, ZZ, XY, YZ, XZ


class PrincipalCorner(IntEnum):
    """
    Triaxial corners of isotropic yield surfaces in principal stress space, for principal stresses
    ordered s1 >= s2 >= s3.
    """
    NONE = 0
    S2_EQ_S3 = 1  # s2 = s3 (triaxial compression when compression is positive)
    S1_EQ_S2 = 2  # s1 = s2 (triaxial extension when compression is positive)


class StressUtils:

    @staticmethod
    def voigt_to_matrix(stress):
        return np.array([[stress[XX], stress[XY], stress[XZ]],
                         [stress[XY], stress[YY], stress[YZ]],
                         [stress[XZ], stress[YZ], stress[ZZ]]])

    @staticmethod
    def matrix_to_voigt(m):
        return np.array([m[0, 0], m[1, 1], m[2, 2], m[0, 1], m[1, 2], m[0, 2]])

    @staticmethod
    def p(stress):
        """Mean stress p = trace(sigma) / 3."""
        return (stress[XX] + stress[YY] + stress[ZZ]) / 3.0

    @staticmethod
    def q(stress):
        """Von Mises equivalent stress q = sqrt(3 J2), including the shear stresses."""
        s_dev = StressUtils.voigt_to_matrix(stress) - StressUtils.p(stress) * np.eye(3)
        return np.sqrt(1.5 * np.sum(s_dev * s_dev))

    @staticmethod
    def principal_system(stress):
        """
        Principal stresses and directions.

        :return: principal stresses sorted descending (s1 >= s2 >= s3), Q with the corresponding
                 directions as columns
        """
        eig, Q = np.linalg.eigh(StressUtils.voigt_to_matrix(stress))
        return eig[::-1], Q[:, ::-1]

    @staticmethod
    def stress_from_principal_system(principal_stress, Q):
        """stress_ij = sum_k principal_stress_k Q_ik Q_jk, in Voigt notation."""
        return StressUtils.matrix_to_voigt((Q * principal_stress) @ Q.T)

    @staticmethod
    def min_principal_stress(stress):
        """Minor principal stress s3."""
        return np.linalg.eigvalsh(StressUtils.voigt_to_matrix(stress))[0]
