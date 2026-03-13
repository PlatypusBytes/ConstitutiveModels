import numpy as np


class ConstitutiveModelUtils:


    @staticmethod
    def calculate_implicit_dlambda_two_yields(dFdsigma1, dFdsigma2, De_dg_dsigma1, De_dg_dsigma2, dFdh1,dFdh2, dhdlambda1,dhdlambda2, F1, F2, active_1, active_2, dlambda1, dlambda2):
        """
        Note that no full coupling between the two yield functions is assumed. no cross derivative for hardening variables is included

        :param dFdsigma1:
        :param dFdsigma2:
        :param De_dg_dsigma1:
        :param De_dg_dsigma2:
        :param dFdh1:
        :param dFdh2:
        :param dhdlambda1:
        :param dhdlambda2:
        :param F1:
        :param F2:
        :return:
        """


        # Build system of equations for plastic multipliers dλ_s, dλ_c
        # We solve: A * [dλ_s; dλ_c] = [F_s; F_c]
        A = np.zeros((2, 2))
        rhs = np.zeros(2)


        if active_1:
            # ∂F_s/∂σ : De : m_s
            a_11 = np.einsum('ij,ij', dFdsigma1, De_dg_dsigma1)
            # ∂F_s/∂γ_p * h_s
            a_11 -= dFdh1 * dhdlambda1

            A[0, 0] = a_11
            rhs[0] = F1

            if active_2:
                a_12 = np.einsum('ij,ij', dFdsigma1, De_dg_dsigma2)
                A[0, 1] = a_12


        if active_2:
            if active_1:
                a_21 = np.einsum('ij,ij', dFdsigma2, De_dg_dsigma1)
                A[1, 0] = a_21

            a_22 = np.einsum('ij,ij', dFdsigma2, De_dg_dsigma2)
            # ∂F_c/∂p_p * h_c
            a_22 -= dFdh2 * dhdlambda2
            A[1, 1] = a_22
            rhs[1] = F2

        # Solve for dλ_s, dλ_c
        try:
            dlambda = np.linalg.solve(A, rhs)
        except np.linalg.LinAlgError:
            # Singular matrix – fallback to diagonal
            dlambda = np.zeros(2)
            if active_1:
                dlambda[0] = rhs[0] / (A[0, 0] + 1e-12)
            if active_2:
                dlambda[1] = rhs[1] / (A[1, 1] + 1e-12)

        # dlambda_1_tmp = dlambda[0]
        # dlambda_2_tmp = dlambda[1]
        dlambda_1_tmp = dlambda[0] if active_1 else 0.0
        dlambda_2_tmp = dlambda[1] if active_2 else 0.0
        #
        # dlambda_1_tmp = np.max([dlambda_1_tmp, 0.0])
        # dlambda_2_tmp = np.max([dlambda_2_tmp, 0.0])

        dlambda1 += dlambda_1_tmp
        dlambda2 += dlambda_2_tmp

        return dlambda1, dlambda2, dlambda_1_tmp, dlambda_2_tmp
