"""
Python counterpart of c_models/strain_utils.c.
"""

import numpy as np

from python_prototypes.utils import XX, YY, ZZ


class StrainUtils:

    @staticmethod
    def volumetric_strain(strain):
        return strain[XX] + strain[YY] + strain[ZZ]

    @staticmethod
    def strain_from_principal_system(principal_strain, Q):
        """Voigt strain from principal strains and directions, with engineering shear strains."""
        m = (Q * principal_strain) @ Q.T
        return np.array([m[0, 0], m[1, 1], m[2, 2], 2.0 * m[0, 1], 2.0 * m[1, 2], 2.0 * m[0, 2]])

    @staticmethod
    def void_ratio(e0, eps_v):
        """Void ratio for a volumetric strain eps_v (compression positive), Eq. 39."""
        return (1.0 + e0) * np.exp(-eps_v) - 1.0
