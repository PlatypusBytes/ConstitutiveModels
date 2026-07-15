
import numpy as np

class YieldFunctions:

    @staticmethod
    def HS_cone_F_s( q, gamma_p, Ei, Eur, qa):
        """Cone (shear) yield function."""
        if q <= 0.0:
            return -gamma_p

        if q >= 0.99 * qa:
            # stay away from asymptote
            return 1e10
        term = 2*q / (1.0 - q/ qa)
        return (term / Ei - 2 * q / Eur - gamma_p)



    @staticmethod
    def F_mohr_coulomb(sigma, phi, c):
        """Mohr-Coulomb yield function."""
        f = (2 * sigma[0] - sigma[1] - sigma[2])/4  - (2* sigma[0] + sigma[1] + sigma[2])/4 * np.sin(phi) - c * np.cos(phi)
        return f

    @staticmethod
    def dfds_mohr_coulomb(sigma, phi):
        """Derivative of the Mohr-Coulomb yield function with respect to the stress tensor."""
        dfdsigma = np.zeros(6)
        dfdsigma[0] = 0.5 - 0.5 * np.sin(phi)
        dfdsigma[1] = -0.25 - 0.25 * np.sin(phi)
        dfdsigma[2] = -0.25 - 0.25 * np.sin(phi)
        return dfdsigma



    @staticmethod
    def d_HS_cone_dh():
        """Derivative of the cone yield function with respect to the hardening variable gamma_p."""
        return -1

    @staticmethod
    def dq_dsigma(sigma, p, q):
        """
        Derivative of the von Mises equivalent stress q with respect to the
        stress tensor, expressed in engineering Voigt notation.

        For a scalar function of stress, the energy-conjugate Voigt gradient
        carries a factor 2 on the shear components so that plain Voigt dot
        products (e.g. dF/dsigma . De . dg/dsigma) reproduce the correct tensor
        contraction.
        """
        if q < 1e-12:
            return np.zeros(6)
        s = sigma - p * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        dq = (3.0 / (2.0 * q)) * s
        dq[3:] *= 2.0  # engineering shear factor
        return dq

    @staticmethod
    def HS_cap_F_c(sigma, p, p_p, phi, M):
        """
        Cap (volumetric) yield function.

        Elliptic cap in the p-q~ plane:
            F_c = (q~ / M)^2 + p^2 - p_p^2
        where q~ is the special deviatoric measure used in the Hardening Soil
        model and p_p is the pre-consolidation pressure (cap size).
        """
        alpha = (3 + np.sin(phi)) / (3 - np.sin(phi))  # cap shape parameter
        q_special = sigma[0] + (alpha - 1) * sigma[1] - alpha * sigma[2]

        return (q_special / M) ** 2 + p ** 2 - p_p ** 2


    @staticmethod
    def d_HS_cap_dh(p_p):
        """Derivative of the cap yield function with respect to the hardening variable p_p."""
        return -2 * p_p

    @staticmethod
    def cap_df_dsigma(sigma, p, q, M, phi):
        """
        Gradient of the cap yield function (associated flow) with respect to the
        stress, consistent with ``HS_cap_F_c``:

            F_c = (q~ / M)^2 + p^2 - p_p^2

        with q~ = s0 + (alpha - 1) s1 - alpha s2 and p = tr(sigma)/3.
        """
        alpha = (3 + np.sin(phi)) / (3 - np.sin(phi))
        q_special = sigma[0] + (alpha - 1) * sigma[1] - alpha * sigma[2]

        dq_special = np.array([1.0, alpha - 1.0, -alpha, 0.0, 0.0, 0.0])
        dp_dsigma = (1.0 / 3.0) * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

        return (2.0 * q_special / M ** 2) * dq_special + 2.0 * p * dp_dsigma


    @staticmethod
    def drucker_prager_df_dsigma(sigma, p, q, phi, psi, p_t):
        """
        Compute the derivative of the Drucker-Prager plastic potential with respect to the stress tensor.



        :param sigma: Stress tensor (3x3 numpy array)
        :param phi: Friction angle in radians
        :param p_t: c cotangent of the friction angle (tension cutoff)
        :return:
        """

        eta = q / (p + p_t)
        sin_phi_mob = 3.0 * eta / (6.0 + eta)
        # Guard against round-off / large-eta states pushing the argument of
        # arcsin outside its valid [-1, 1] domain.
        sin_phi_mob = np.clip(sin_phi_mob, -1.0, 1.0)


        sin_phi = np.sin(phi)
        sin_psi = np.sin(psi)

        sin_phi_cs = (sin_phi - sin_psi) / (1 - sin_phi * sin_psi)

        if np.arcsin(sin_phi_mob) < phi:
            sin_psi_mob = 0
        else:
            sin_psi_mob = (sin_phi_mob - sin_phi_cs) / \
                          (1 - sin_phi_mob * sin_phi_cs)


        if q < 1e-12:
            return np.zeros(6)

        M = YieldFunctions.mohr_coulomb_to_drucker_prager(sin_psi_mob, 0)[0]

        # Deviatoric part (engineering Voigt gradient of q) minus the
        # volumetric part applied to the identity vector only.
        dq_dsigma = YieldFunctions.dq_dsigma(sigma, p, q)
        identity = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
        dgdsigma = dq_dsigma - (M / 3.0) * identity

        return dgdsigma

    @staticmethod
    def dgd_dsigma(sigma, p, q):

        s = sigma - p * np.array([1.0,1.0,1.0,0,0,0])

        # Identity in Voigt
        I6 = np.eye(6)

        # J matrix
        J = np.zeros((6, 6))
        J[:3, :3] = 1.0

        Pdev_voigt = I6 - (1 / 3) * J

        term1_voigt = (3.0 / (2.0 * q)) * Pdev_voigt
        term2_voigt = (9.0 / (4.0 * q ** 3)) * np.outer(s, s)

        return term1_voigt - term2_voigt


    @staticmethod
    def mohr_coulomb_to_drucker_prager( phi, c):
        """
        Convert Mohr-Coulomb parameters to Drucker-Prager parameters.
        :param phi: Friction angle in rad
        :param c: Cohesion
        :return: M (slope of the DP yield surface) and p_t (tension cutoff)
        """
        M = 6 * np.sin(phi) / (3 - np.sin(phi))
        p_t = 6 * c * np.cos(phi) / (3 - np.sin(phi))
        return M, p_t

    # @staticmethod
    # def get_fs_gradient(sigma,q,p, p_t, Ei, Eur, Rf, phi):
    #
    #     dp_dsigma = 1/3 * np.eye(3)
    #     s = sigma - p * np.eye(3)
    #     dq_dsigma = (3.0 / 2.0) * s / q
    #     # Constant term
    #     fact_phi = (1.0 - np.sin(phi)) / np.sin(phi)
    #     C = Rf * fact_phi
    #
    #     # Simplified fact_mob based on eta
    #     eta = q / (p + p_t)
    #     fact_mob = (2.0 / eta) - (2.0 / 3.0)
    #
    #     # 1. Partial derivative wrt p
    #     denom = (fact_mob - C) ** 2
    #     dfs_dp = -3.0 * C / (Ei * denom)
    #
    #     # 2. Partial derivative wrt q
    #     dfs_dq = (3.0 / (2.0 * Ei)) * (fact_mob ** 2 + 2.0 * C / 3.0) / denom - (3.0 / (2.0 * Eur))
    #
    #     # 3. Apply chain rule for full tensor derivative
    #     # Note: If dp_dsigma and dq_dsigma are arrays (e.g., 3x3 or 6x1),
    #     # this will correctly broadcast to a matching tensor shape.
    #     dfs_dsigma = (dfs_dp * dp_dsigma) + (dfs_dq * dq_dsigma)
    #     return dfs_dsigma

    @staticmethod
    def df_dsigma_mohr_coulomb(sigma,p, q, Ei, Eur, qa):
        """
        Gradient of the Hardening Soil cone (shear) yield function with respect
        to the stress tensor, in engineering Voigt notation.

        :param sigma:  Stress vector (engineering Voigt)
        :param Ei:  Primary loading stiffness
        :param Eur:  Unloading/reloading stiffness
        :param qa:  Asymptotic deviator stress
        :return:  Gradient of the yield function w.r.t. stress (6-vector)
        """
        if q < 1e-12:
            return np.zeros(6)

        dq_dsigma = YieldFunctions.dq_dsigma(sigma, p, q)

        denom = 1.0 - q / qa
        if denom <= 1e-12:
            denom = 1e-12

        dF_dq = 2.0 / (Ei * denom ** 2) - 2.0 / Eur

        return dF_dq * dq_dsigma
