
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


        s = sigma_voigt[:3] - StressUtils.p(sigma_voigt)
        J2 = 0.5 * s @ s
        return np.sqrt(3.0 * J2)