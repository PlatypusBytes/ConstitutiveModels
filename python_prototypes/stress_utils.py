
import numpy as np

class StressUtils:

    @staticmethod
    def p(sigma):
        return np.trace(sigma) / 3.0

    @staticmethod
    def q(sigma):
        s = sigma - StressUtils.p(sigma) * np.eye(3)
        J2 = 0.5 * np.einsum('ij,ij', s, s)

        # return sigma[0,0] - sigma[2,2]
        return np.sqrt(3.0 * J2)
