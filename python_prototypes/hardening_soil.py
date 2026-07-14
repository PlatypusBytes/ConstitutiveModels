import numpy as np

from python_prototypes.utils import ConstitutiveModelUtils
from python_prototypes.stress_utils import StressUtils
from python_prototypes.hardening_laws import HardeningLaws
from python_prototypes.yield_functions import YieldFunctions

def matrix_to_voigt(m):
    return np.array([m[0,0], m[1,1], m[2,2], m[1,2], m[0,2], m[0,1]])

def voigt_to_matrix(v):
    m = np.zeros((3,3))
    m[0,0] = v[0]
    m[1,1] = v[1]
    m[2,2] = v[2]
    m[1,2] = m[2,1] = v[3]
    m[0,2] = m[2,0] = v[4]
    m[0,1] = m[1,0] = v[5]
    return m

class HardeningSoilWithCap:
    """
    Hardening Soil model with cone hardening (shear) and cap hardening (volumetric).
    No tension cutoff for simplicity.
    """

    def __init__(self, params):
        """
        Parameters
        ----------
        params : dict
            Material parameters:
                E50_ref  : reference secant stiffness at 50% strength [kPa]
                Eur_ref  : reference unloading/reloading stiffness [kPa]
                m        : stress exponent for stiffness
                c        : cohesion [kPa]
                phi      : friction angle [degrees]
                psi      : dilation angle [degrees]
                pref     : reference stress (usually 100 kPa)
                Rf       : failure ratio (0.9 typical)
                nu       : Poisson's ratio (default 0.2)
                M        : cap shape parameter (slope of cap ellipse, default 1.0)
                H_cap    : cap hardening modulus (positive) [kPa]
        """
        self.E50_ref = params['E50_ref']
        self.Eur_ref = params['Eur_ref']

        self.E50_current = self.E50_ref
        self.Eur_current = self.Eur_ref
        self.m = params['m']
        self.c = params['c']
        self.phi = np.radians(params['phi'])
        self.psi = np.radians(params['psi'])
        self.pref = params['pref']
        self.Rf = params['Rf']
        self.nu = params.get('nu', 0.2)
        self.M = params.get('M', 1.0)          # cap aspect ratio
        self.H_cap = params.get('H_cap', 1e4)  # cap hardening modulus
        self.e_init = params.get('e_init', 0.5)  # initial void ratio
        self.e_max = params.get('e_max', 1.0)    # maximum void ratio (loosest state)
        self.e = self.e_init
        self.k0nc = params.get('k0nc', 1 - np.sin(self.phi))  # at-rest earth pressure coefficient for normally consolidated state
        self.K_ratio = params.get('K_ratio', 0.7)  # ratio of bulk to shear stiffness for elastic part

        self.use_mohr_coulomb =False

        # Derived constants
        if self.phi > 0:
            self.cot_phi = 1.0 / np.tan(self.phi)
            self.p_t = self.c * self.cot_phi   # tensile strength (positive)
        else:
            self.cot_phi = 1e12
            self.p_t = 0.0
        # self.Ei_ref = self.E50_ref
        self.Ei_ref = 65488
        # self.Ei_ref = self.E50_ref * (3 - self.Rf)/(2-self.Rf)  # reference initial stiffness (for elastic predictor)

        # State variables
        self.sigma = None          # current stress (3x3 tensor)
        self.gamma_p = 0.0         # plastic shear strain (cone hardening)
        self.p_p = None            # preconsolidation pressure (cap hardening)
        self.initialized = False

        self.recalculate_E50 = True
        self.recalculate_Eur = True
        self.recalculate_Ei = True

    # ------------------------------------------------------------------
    # Stress invariants and derived quantities
    # ------------------------------------------------------------------


    def _Ei(self, sigma_voigt):
        """Stress‑dependent initial stiffness for elastic predictor."""
        # p = StressUtils.p(sigma)
        denom = self.p_t + self.pref
        if denom == 0.0:
            denom = 1e-12
        if self.recalculate_Ei:
            res_1= self.Ei_ref * ((sigma_voigt[2]  + self.p_t)/(self.pref + self.p_t) )** self.m
            self.Ei_current = res_1
        else:
            res_1 = self.Ei_current

        self.recalculate_Ei = False

        return res_1


    def _E50(self, sigma_voigt):
        """Stress‑dependent primary loading stiffness (used in cone yield)."""

        denom = self.p_t + self.pref
        if denom == 0.0:
            denom = 1e-12

        if self.recalculate_E50:

            res_1 = self.E50_ref * ((sigma_voigt[2]  + self.p_t)/(self.pref + self.p_t) )** self.m

            self.E50_current = res_1
        else:
            res_1 =self.E50_current

        self.recalculate_E50 = False

        return res_1

    def _Eur(self, sigma):
        """Stress‑dependent unloading/reloading stiffness."""
        denom = self.p_t + self.pref
        if denom == 0.0:
            denom = 1e-12

        ratio =  (sigma[2]  + self.p_t)/(self.pref + self.p_t)
        if ratio <= 0.0:
            ratio = 1e-12

        if self.recalculate_Eur:

            res_1 = self.Eur_ref * ((sigma[2] + self.p_t) / (self.pref + self.p_t) )** self.m
            self.Eur_current = res_1
        else:
            res_1 = self.Eur_current

        self.recalculate_Eur = False

        return res_1

    def _qf(self, sigma_voigt):
        """Failure deviator stress (triaxial compression)."""
        sin_phi = np.sin(self.phi)

        qf  = (2 * sin_phi / (1 - sin_phi)) * (sigma_voigt[2] + self.p_t)
        return qf

    def _sin_phi_mob(self, sigma_voigt):
        sin_phi_mob = (sigma_voigt[0] - sigma_voigt[2]) / (sigma_voigt[0] + sigma_voigt[2] + 2* self.p_t)
        return sin_phi_mob

    def _qa(self, sigma):
        """Asymptotic deviator stress for the hyperbola."""

        qf =  self._qf(sigma)
        qa = qf / self.Rf

        return qa

    # ------------------------------------------------------------------
    # Elastic stiffness (using Eur)
    # ------------------------------------------------------------------
    def elastic_moduli(self, Eur, nu):
        # Eur = self._E50(sigma)
        G = Eur / (2.0 * (1.0 + nu))
        K = Eur / (3.0 * (1.0 - 2.0 * nu))
        return G, K

    def elastic_product(self, dqdsigma, Eur, nu):
        """Compute De : m (tensor)."""
        G, K = self.elastic_moduli(Eur, nu)
        trace_dqdsigma = dqdsigma[0] + dqdsigma[1] + dqdsigma[2]
        # trace_dqdsigma = np.trace(dqdsigma)
        lame_lambda = K - 2.0 * G / 3.0
        return 2.0 * G * dqdsigma + lame_lambda * trace_dqdsigma * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])



    # ------------------------------------------------------------------
    # Initialisation and integration driver
    # ------------------------------------------------------------------
    def set_initial_state(self, sigma0, p_p0=None, ocr=1.0):
        """
        Set initial stress and preconsolidation pressure.
        If p_p0 is not given, it is set so that the initial stress lies on or inside the cap.
        """
        self.sigma = np.array(sigma0, dtype=float)
        self.gamma_p = 0.0

        p = StressUtils.p(self.sigma)
        q = StressUtils.q(self.sigma)

        qa = self._qa(self.sigma)
        Ei = self._Ei(self.sigma)
        # E50 = self._E50(self.sigma)

        term = 2*q / (1.0 - q / qa)
        self.gamma_p =  term / Ei - 2 * q / self._Eur(sigma0)

        if p_p0 is None:
            # Place cap just outside current stress to start elastically
            p0 = StressUtils.p(sigma0)
            q0 = StressUtils.q(sigma0)
            # Ensure F_c <= 0
            # rhs = np.sqrt((q0/self.M)**2 + (p0 + self.p_t)**2) - self.p_t
            fact = (q0 / self.M)**2  - (p0) ** 2
            if fact <= 0:
                self.p_p = 1.0  # small positive value to ensure we start inside the cap
            else:
                rhs = np.sqrt((q0 / self.M)**2  - (p0) ** 2)
                self.p_p = max(rhs, 1.0)  # at least a small positive value
        else:
            self.p_p = p_p0
        self.initialized = True

    def integrate(self, deps, tol=1e-5, max_iter=20, max_subdivisions=20):
        """
        Integrate stress for a given strain increment with automatic substepping.
        Handles combined cone‑cap activation using multi‑surface cutting plane.
        """
        if not self.initialized:
            raise ValueError("Initial state not set. Call set_initial_state first.")

        sigma = self.sigma.copy()

        self.update_variables(sigma)

        gamma_p = self.gamma_p
        p_p = self.p_p



        q= StressUtils.q(sigma)

        qf = self._qf(sigma)

        # min_d_deviatoric = 0.01
        # q_step = min_d_deviatoric*qf
        # max_it1 = np.ceil(q / q_step)

        delta_sigma_el = self._calculate_dsigma_elastic(self._Eur(sigma), self.nu, deps)
        # qa_trial = self._qa(sigma + delta_sigma_el)
        # q_trial = StressUtils.q(sigma + delta_sigma_el)
        # determine substep size based on proximity to yield surface to prevent overshooting
        alpha = 0.9
        cutoff_ratio = 0.99
        ratio = q / qf
        q_step = alpha * qf * (1 - min(ratio, cutoff_ratio))
        max_it2 = np.ceil(q / q_step)

        min_it = 5
        max_it =max(max_it2, min_it)

        if self.use_mohr_coulomb:
            max_it=1
        # it = 0
        # print(max_it)
        max_trials = 5
        n_trial = 0
        converged = False
        while not converged and n_trial < max_trials:
            # max_it = max_it * 2** n_trial
            max_it = max_it * (n_trial+1)
            it = 0
            sigma_0 = sigma.copy()
            gamma_p_0 = gamma_p
            p_p_0 = p_p
            while  it < max_it:
                it += 1

                deps_step = deps /max_it
                self.update_variables(sigma)
                # Perform plastic correction on this substep
                sigma, gamma_p, p_p, converged = self._single_step_euler_backward(sigma, gamma_p, p_p, deps_step, tol, max_iter)

                if not converged:
                    q = StressUtils.q(sigma)
                    Ei = self._Ei(sigma)
                    Eur = self._Eur(sigma)
                    f = YieldFunctions.HS_cone_F_s(q, gamma_p, Ei, Eur, self._qa(sigma))


                    sigma = sigma_0.copy()
                    gamma_p = gamma_p_0
                    p_p = p_p_0
                    break

            n_trial+=1

        if not converged:
            raise RuntimeError(f"Integration did not converge after {n_trial} trials with max_it up to {max_it}.")


        self.sigma = sigma
        self.gamma_p = gamma_p
        self.p_p = p_p
        return self.sigma



    def _calculate_dsigma_elastic(self, Eur, nu, deps):
        G, K = self.elastic_moduli(Eur, nu)

        deps_v = np.sum(deps[:3])

        lame_lambda = K - 2.0 * G / 3.0

        dsigma_elastic = 2.0 * G * deps
        dsigma_elastic[:3] += lame_lambda * deps_v
        return dsigma_elastic


    def update_variables(self, sigma):

        self.recalculate_Ei = True
        self.recalculate_E50 = True
        self.recalculate_Eur = True

        return self._Ei(sigma), self._E50(sigma), self._Eur(sigma)



    def update_stress_and_state_2_yields_implicit(self, sigma_trial, active_1, active_2, dlambda_1, dlambda_2,dg_dsigma1, dg_dsigma2, De_dg_dsigma1,
                                                  De_dg_dsigma2, gamma_p_0,
                                                  p_p_0):
        # Update stress and hardening variables
        sigma_curr_tmp = np.copy(sigma_trial)
        d_eps_p = np.zeros(6)

        gamma_p_curr = gamma_p_0
        p_p_curr = p_p_0

        if active_1:
            sigma_curr_tmp = sigma_curr_tmp - dlambda_1 * De_dg_dsigma1
            dhdlambda_1 = HardeningLaws.dh_dlambda_s(dg_dsigma1)

            d_eps_p = dlambda_1 * dg_dsigma1
            gamma_p_curr = gamma_p_0 + dlambda_1 * dhdlambda_1



        if active_2:
            sigma_curr_tmp = sigma_curr_tmp - dlambda_2 * De_dg_dsigma2

            dhdlambda_2 = HardeningLaws.dh_dlambda_c(sigma_curr_tmp,self.Eur_ref, self.nu, self.m, self.pref, self.p_t, self.K_ratio)
            # d_eps_v = dlambda_2 * dhdlambda_2


            p_p_curr += dlambda_2 * dhdlambda_2
            d_eps_p = d_eps_p + dlambda_2 * dg_dsigma2

        # d_eps_p_vol  = np.trace(d_eps_p)
        d_eps_p_vol = d_eps_p[0] + d_eps_p[1] + d_eps_p[2]
        # p_p_curr = HardeningLaws.dpreconsolidation_stress(d_eps_p_vol,  self.nu, self.Eur_ref, self.K_ratio, p_p_0)
        # gamma_p_curr = gamma_p_0 + HardeningLaws.dgamma(d_eps_p)

        sigma_curr = sigma_curr_tmp
        gamma_p_curr = max(gamma_p_curr, gamma_p_0)
        p_p_curr = max(p_p_curr, p_p_0)  # ensure preconsolidation pressure does not decrease

        return sigma_curr, d_eps_p, gamma_p_curr, p_p_curr



    def _single_step_euler_backward(self, sigma, gamma_p, p_p, deps, tol, max_iter):
        # Elastic trial

        sigma_0 = sigma.copy()
        gamma_p_0 = gamma_p
        p_p_0 = p_p

        # Ei, E50, Eur = self.update_variables(sigma)
        Ei = self.Ei_current
        E50 = self.E50_current
        Eur = self.Eur_current

        delta_sigma_el = self._calculate_dsigma_elastic(Eur, self.nu, deps)
        gamma_p_trial, p_p_trial = gamma_p, p_p
        gamma_p0, p_p0 = gamma_p, p_p
        sigma_trial = sigma + delta_sigma_el




        # sigma_trial, gamma_p_trial, p_p_trial = self._elastic_predictor(sigma, gamma_p, p_p, deps)
        p_trial = StressUtils.p(sigma_trial)
        q_trial = StressUtils.q(sigma_trial)

        if not self.use_mohr_coulomb:
            f_s_trial = YieldFunctions.HS_cone_F_s(q_trial, gamma_p_trial, Ei, Eur, self._qa(sigma_trial))
        else:
            f_s_trial = YieldFunctions.F_mohr_coulomb(sigma_trial, self.phi, self.c)

        f_c_trial = YieldFunctions.HS_cap_F_c(sigma_trial,p_trial, p_p_trial, self.phi, self.M)
        # Check if any yielding
        # tol1= tol
        tol1 = abs(f_s_trial) * tol
        tol2 = abs(f_c_trial) * tol
        if f_s_trial <= tol1 and f_c_trial <= tol2:
            return sigma_trial, gamma_p_trial, p_p_trial, True

        # Initialize
        sigma_curr = sigma_trial.copy()

        gamma_p_curr = gamma_p0
        p_p_curr = p_p0

        d_eps_p = np.zeros(6)

        # initial guess for yield function values (using trial state)
        f_s_curr = f_s_trial
        f_c_curr = f_c_trial

        active_s = f_s_curr > tol1
        active_c = f_c_curr > tol2
        active_c = False

        converged = False

        for i in range(2):
            # first only check 1 surface, then both.
            if i == 0 and active_s:
                active_c = False
            if not active_s and not active_c:
                return sigma_curr, gamma_p_curr, p_p_curr, True

            dlambda_s = 0
            dlambda_c = 0

            for it in range(max_iter):

                deps = deps - d_eps_p
                d_eps_p = np.zeros(6)

                # Compute necessary tensors and scalars
                p_curr = StressUtils.p(sigma_curr)
                q_curr = StressUtils.q(sigma_curr)

                qa = self._qa(sigma_curr)

                dg_dsigma_s, dg_dsigma_c = None, None
                De_dg_dsigma_s, De_dg_dsigma_c = None, None
                dF_s_dsigma, dF_c_dsigma = None, None
                dF_s_dh, dF_c_dh = None, None
                dh_dlambda_s, dh_dlambda_c = None, None

                if active_s:
                    dg_dsigma_s = YieldFunctions.drucker_prager_df_dsigma(sigma_curr, p_curr, q_curr, self.phi,self.psi, self.p_t)
                    dg_ddsigma_s = YieldFunctions.dgd_dsigma(sigma_curr,p_curr,q_curr )
                    # fourth_order_identity = np.einsum('ik,jl->ijkl', np.eye(3), np.eye(3))
                    #
                    G, K= self.elastic_moduli(Eur, self.nu)
                    # De =fourth_order_identity * K + 2 * G * (np.einsum('ik,jl->ijkl', np.eye(3), np.eye(3)) - (1/3) * np.einsum('ij,kl->ijkl', np.eye(3), np.eye(3)))
                    #
                    # A = fourth_order_identity + dlambda_s * np.einsum(
                    #     "ijkl,klmn->ijmn",
                    #     De,
                    #     dg_ddsigma_s
                    # )
                    #
                    #
                    # term1 = np.linalg.inv(A.reshape(9, 9)).reshape(3, 3, 3, 3)
                    # term2 = dg_dsigma_s
                    De_dg_dsigma_s = self.elastic_product(dg_dsigma_s, Eur, self.nu)
                    #
                    # res = sigma_curr - (sigma_trial - dlambda_s * De_dg_dsigma_s)

                    # --- Elastic matrix ---
                    lame_lambda = K - 2.0 * G / 3.0

                    De_voigt = np.zeros((6, 6))
                    De_voigt[:3, :3] = lame_lambda
                    np.fill_diagonal(De_voigt[:3, :3], lame_lambda + 2 * G)
                    De_voigt[3:, 3:] = 2* G * np.eye(3)

                    # --- A matrix ---
                    A_voigt = np.eye(6) + dlambda_s * (De_voigt @ dg_ddsigma_s)
                    # A_voigt = np.eye(6) + dlambda_s * De_dg_ddsigma_s

                    # term1 = np.linalg.inv(A_voigt)



                    # --- Elastic product ---
                    # De_dg_dsigma_s = De_voigt @ dg_dsigma_s

                    # --- Residual ---
                    res = sigma_curr - (sigma_trial - dlambda_s * De_dg_dsigma_s)

                    dF_s_dh = YieldFunctions.d_HS_cone_dh()
                    dh_dlambda_s = HardeningLaws.dh_dlambda_s(dg_dsigma_s)
                    if not self.use_mohr_coulomb:
                        dF_s_dsigma = YieldFunctions.df_dsigma_mohr_coulomb(sigma_curr, p_curr, q_curr, Ei, Eur, qa)
                        Ainv_dF_s_dsigma = np.linalg.solve(A_voigt, dF_s_dsigma)
                        if it > 0:
                            # dd_lambda_s = f_s_curr / (np.einsum('ij,ij', dF_s_dsigma, De_dg_dsigma_s) - dF_s_dh * dh_dlambda_s)
                            # dd_lambda_s = f_s_curr / (dF_s_dsigma  @ De_dg_dsigma_s - dF_s_dh * dh_dlambda_s)

                            dd_lambda_s = (f_s_curr- res @ (Ainv_dF_s_dsigma)) / (De_dg_dsigma_s @ (Ainv_dF_s_dsigma) - dF_s_dh * dh_dlambda_s ) * 0.5

                            # dd_lambda_s = ((f_s_curr - np.einsum('ij,ijmn,mn->', res, term1, dF_s_dsigma))
                            #                / (np.einsum('ijkl,kl,ijmn,mn->', De, term2, term1,
                            #                             dF_s_dsigma) - dF_s_dh * dh_dlambda_s))

                            dlambda_s += dd_lambda_s
                    else:
                        dg_dsigma_s = YieldFunctions.dfds_mohr_coulomb(sigma,self.psi)
                        dF_s_dsigma = YieldFunctions.dfds_mohr_coulomb(sigma, self.phi)
                        dF_s_dh = 0

                        De_dg_dsigma_s = self.elastic_product(dg_dsigma_s, Eur, self.nu)

                        if it> 0:
                            # dd_lambda_s = ((f_s_curr - np.einsum('ij,ijmn,mn->', res, term1, dF_s_dsigma))
                            #                / (np.einsum('ijkl,kl,ijmn,mn->', De, term2, term1,
                            #                             dF_s_dsigma) - dF_s_dh * dh_dlambda_s))

                            A = dF_s_dsigma @ De_dg_dsigma_s - dF_s_dh * dh_dlambda_s
                            dd_lambda_s = f_s_curr / A
                            dlambda_s += dd_lambda_s

                if active_c:
                    dg_dsigma_c = YieldFunctions.cap_df_dsigma(sigma_curr,p_curr, q_curr, self.M)
                    De_dg_dsigma_c = self.elastic_product(dg_dsigma_c, Eur, self.nu)
                    dF_c_dsigma = YieldFunctions.cap_df_dsigma(sigma_curr,p_curr, q_curr, self.M)

                    dF_c_dh = YieldFunctions.d_HS_cap_dh(p_p_curr)
                    dh_dlambda_c = HardeningLaws.dh_dlambda_c(sigma_curr, self.Eur_ref, self.nu, self.m, self.pref,
                                                              self.p_t, self.K_ratio)


                if it == 0:
                    dlambda_s,dlambda_c, res_1, res_2 =  ConstitutiveModelUtils.calculate_implicit_dlambda_two_yields(dF_s_dsigma, dF_c_dsigma,
                                                                                                        De_dg_dsigma_s, De_dg_dsigma_c,
                                                                                                        dF_s_dh,dF_c_dh, dh_dlambda_s,dh_dlambda_c, f_s_curr, f_c_curr, active_s, active_c, dlambda_s, dlambda_c)

                sigma_curr, d_eps_p, gamma_p_curr, p_p_curr = self.update_stress_and_state_2_yields_implicit(sigma_trial, active_s, active_c, dlambda_s, dlambda_c, dg_dsigma_s,
                                                                                                              dg_dsigma_c, De_dg_dsigma_s,
                                                                                                              De_dg_dsigma_c, gamma_p0,
                                                                                                              p_p0)

                # check if we are within tolerance
                p_curr = StressUtils.p(sigma_curr)
                q_curr = StressUtils.q(sigma_curr)
                qa = self._qa(sigma_curr)
                qf = self._qf(sigma_curr)
                Ei, E50, Eur = self.update_variables(sigma_curr)

                if active_s:

                    q = q_curr

                    if not self.use_mohr_coulomb:
                        f_s_curr = YieldFunctions.HS_cone_F_s(q_curr, gamma_p_curr, Ei, Eur, qa)
                    else:
                        f_s_curr = YieldFunctions.F_mohr_coulomb(sigma_curr, self.phi, self.c)

                        # no hardening for mohr coulomb, so if we are beyond qf, we set gamma_p to the value at qf
                        if q > qf:
                            gamma_p_curr = HardeningLaws.maximum_gamma_p(qf, qa, Ei, Eur)

                    # f_s_curr = YieldFunctions.HS_cone_F_s(q_curr, gamma_p_curr, Ei, Eur, qa)

                # f_s_curr =YieldFunctions.HS_cone_F_s(sigma_curr, p_curr, q_curr, gamma_p_trial, Ei, Eur, self.p_t, self.Rf,
                #                            self.phi)
                if active_c:
                    f_c_curr = YieldFunctions.HS_cap_F_c(sigma_curr, p_curr, p_p_curr, self.phi, self.M)


                if active_s and not active_c:
                    if abs(f_s_curr) < tol1: #or  abs(res_1) < 1e-10: # or abs(f_s_curr) < tol:
                        converged = True
                        active_s = False

                        if q / qf < 0.999:
                            if self.use_mohr_coulomb:
                                print(
                                    f"Switching back to HS cone yield function at q={q:.3e} (qf={qf:.3e}). Final f_s = {f_s_curr:.3e}, f_c = {f_c_curr:.3e}")

                            self.use_mohr_coulomb = False
                        else:
                            if not self.use_mohr_coulomb:
                                print(
                                    f"Switching to Mohr-Coulomb yield function at q={q:.3e} (qf={qf:.3e}). Final f_s = {f_s_curr:.3e}, f_c = {f_c_curr:.3e}")

                            self.use_mohr_coulomb = True

                        break

                if active_c and not active_s:
                    if abs(f_c_curr) < tol2:
                        converged = True
                        break

                if active_s and active_c:
                    if abs(f_s_curr) <= tol1 or abs(f_c_curr) <= tol2:
                        converged = True
                        break

            q = StressUtils.q(sigma_curr)

            f_c_curr = YieldFunctions.HS_cap_F_c(sigma_curr, p_curr, p_p_curr, self.phi, self.M)
            active_s = f_s_curr > tol1

            # if f_c_curr > tol2:
            #     print(f"Cap still active after cutting plane iterations: f_c = {f_c_curr:.3e}, f_s = {f_s_curr:.3e}, q/qf = {q/qf:.3e}")

            # active_c = f_c_curr > tol2


        # not converged after max iterations, return last state and flag
        if it == max_iter-1:
            return sigma_0, gamma_p_0, p_p_0, False

        return sigma_curr, gamma_p_curr, p_p_curr, converged


    def ddsde(self, sigma):
        """Compute consistent tangent operator dσ/dε for given state."""
        # This is a complex derivation involving the active yield surfaces and their derivatives.
        # For simplicity, we return the elastic stiffness here. A full implementation would require
        # differentiating the cutting plane algorithm and accounting for the active surfaces.

        Ei, E50, Eur = self.update_variables(sigma)

        G, K = self.elastic_moduli(Eur, self.nu)
        lame_lambda = K - 2.0 * G / 3.0
        De = np.zeros((6,6))
        De[:3,:3] = lame_lambda + 2.0 * G * np.eye(3)
        De[3:,3:] = G * np.eye(3)
        return De

