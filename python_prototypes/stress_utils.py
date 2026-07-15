
import numpy as np

class StressUtils:

    @staticmethod
    def p(sigma_voigt):
        """
        Calculates the mean stress (pressure) from the voigt notation stress vector.
        """

        return (sigma_voigt[0] + sigma_voigt[1] + sigma_voigt[2]) / 3.0


    @staticmethod
    def q(sigma_voigt):
        """
        Von Mises equivalent (deviatoric) stress q = sqrt(3 * J2).

        Uses the full stress state, including shear components. The stress is
        stored in engineering Voigt notation
        [s_xx, s_yy, s_zz, s_yz, s_xz, s_xy], so the shear stresses appear with
        a unit factor in J2 (they are true tensor shear stresses).
        """
        p = StressUtils.p(sigma_voigt)
        s = sigma_voigt[:3] - p
        J2 = 0.5 * (s @ s) + (
            sigma_voigt[3] ** 2 + sigma_voigt[4] ** 2 + sigma_voigt[5] ** 2
        )
        return np.sqrt(3.0 * J2)