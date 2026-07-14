
import numpy as np

from python_prototypes.stress_utils import StressUtils

class HardeningLaws:


    @staticmethod
    def dh_dlambda_s(dgdsigma):
        """

        hardening law according to Benz

        d_gamma_p/ dlambda_s = sqrt(1/2 (dgdsigma1 - dgdsigma2)**2 + (dgdsigma2 - dgdsigma3)**2 + (dgdsigma3 - dgdsigma1)**2)) = 3/2

        :param dh_dlambda_s:
        :return:
        """
        # return 1
        return 2
        # return 3
        # return np.sqrt(2.0 / 3.0 * np.einsum('i,i', dgdsigma, dgdsigma))
        # return 3/2


    @staticmethod
    def dgamma(eps_plastic):
        """Calculate the equivalent plastic strain increment dγ from the plastic strain increment tensor."""
        # For J2 plasticity, we can use the following formula:
        # dγ = sqrt(2/3 * eps_plastic:eps_plastic)
        return np.sqrt(2.0 / 3.0 * np.einsum('ij,ij', eps_plastic, eps_plastic))

    @staticmethod
    def dpreconsolidation_stress(eps_vol_plastic,  nu, Eu_ref, K_ratio, p0):

        Ks = Eu_ref / (3.0 * (1.0 - 2.0 * nu))
        H = (1/(K_ratio-1)) * Ks

        return p0 * np.exp(-H * eps_vol_plastic)

    @staticmethod
    def j2_dh_dlambda_s(dgdsigma):
        """dγ_p/dλ_s = sqrt(2/3 dgdsigma:dgdsigma)
        strain hardening as written in Sloan et al 2001 (written in terms of invariants)
        (J2 hardening law)
        """
        return np.sqrt(2.0 /3.0 * np.einsum('ij,ij', dgdsigma, dgdsigma))

    @staticmethod
    def dh_dlambda_c(sigma, Eu_ref, nu, m, pref, p_t, K_ratio):
        """
        hardening cap law according to Benz
        :return:
        """

        p = StressUtils.p(sigma)

        Ks = Eu_ref / (3.0 * (1.0 - 2.0 * nu))
        H = (1/(K_ratio-1)) * Ks

        return 2*H *((sigma[2,2] + p_t)/(pref + p_t))**m * p

    @staticmethod
    def maximum_gamma_p(q, qa, Ei, Eur):


        gamma_p = 2*q/Ei * qa/(qa-q) - 2*q/Eur

        return gamma_p








