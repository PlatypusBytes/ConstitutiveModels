"""
Python prototype of the Hardening Soil model of Schanz, Vermeer & Bonnier (1999), "The hardening soil
model: Formulation and verification". Equation numbers refer to that paper.

This is a one-to-one port of c_models/hardening_soil/hardening_soil.c; see that file for the
description of the model, the integration scheme and the additions for boundary value problems.
The functions below carry the names of their C counterparts (without the hs_ prefix).

Differences with the C implementation are limited to the interface:
 - The model is a Python class instead of a UMAT. Stresses and strains are compression positive
   (the internal convention of the C model), Voigt ordering [xx, yy, zz, xy, yz, xz] with
   engineering shear strains.
 - The total volumetric strain (for the void ratio in the dilatancy cut-off) is tracked by the
   model instead of being passed in as STRAN.
 - When the integration fails, integrate() returns False and keeps the stress and the state at the
   start of the increment, with the elastic tangent (the C UMAT then requests PNEWDT = 0.5).

Material parameters (dict keys, or the C PROPS order with HardeningSoil.from_props)
------------------------------------------------------------------------------------
  E50_ref  - reference secant stiffness at 50% strength   [stress]
  Eur_ref  - reference unloading/reloading stiffness       [stress]
  m        - stress exponent for stiffness                 [-]
  c        - cohesion                                      [stress]
  phi      - friction angle                                [degrees]
  psi      - dilation angle                                [degrees]
  p_ref    - reference pressure for stiffness              [stress]
  Rf       - failure ratio qf/qa (typ. 0.9)                [-]
  nu       - Poisson's ratio (unloading/reloading)         [-]
  M_cap    - cap aspect ratio M (Eq. 27), 0 disables the cap [-]
  K_ratio  - Ks/Kc (Eq. 32)                                [-]
  e0       - (optional) initial void ratio (Eq. 39)        [-]
  e_cv     - (optional) critical void ratio of the dilatancy cut-off (Eq. 38), 0 disables it [-]

State variables
---------------
  gamma_p    - plastic shear strain (Eq. 9)
  p_c        - pre-consolidation stress (cap size)
  at_failure - True when the stress is on the Mohr-Coulomb failure surface
"""

from dataclasses import dataclass
from enum import IntEnum
from typing import Optional

import numpy as np

from python_prototypes.elastic_laws import ElasticLaws
from python_prototypes.flow_rules import FlowRules
from python_prototypes.hardening_laws import HardeningLaws
from python_prototypes.strain_utils import StrainUtils
from python_prototypes.stress_utils import PrincipalCorner, StressUtils
from python_prototypes.utils import SMALL_VALUE
from python_prototypes.yield_functions import YieldFunctions

HS_MAX_LOCAL_ITER = 50
HS_MAX_LINE_SEARCH = 30
HS_MAX_ACTIVE_SET_ITER = 10
HS_MAX_SUBSTEPS = 1000
HS_MAX_SUBSTEP_TRIALS = 5

# Maximum elastic stress increment per sub-step, relative to sigma_3 + a.
HS_SUBSTEP_STRESS_FRACTION = 0.01

# Relative tolerance on the yield functions, and the tolerance at which a stagnating local
# iteration (round-off level) is still accepted.
HS_YIELD_TOL = 1.0e-12
HS_YIELD_TOL_ACCEPT = 1.0e-8

# Lower limit of (sigma_3 + a) / (p_ref + a) in the stress dependency of the stiffness.
HS_MIN_STRESS_RATIO = 1.0e-2

# Fraction of the elastic stiffness kept in the tangent for the corner mode and at the apex.
HS_CORNER_STIFFNESS_FRACTION = 1.0e-2

# The C implementation currently returns the elastic stiffness as tangent (the call to
# hs_elastoplastic_tangent is commented out). Set to True to use the elasto-plastic tangent.
HS_USE_ELASTOPLASTIC_TANGENT = False


class SurfaceType(IntEnum):
    CONE = 0    # shear hardening surface f_13, Eqs. 7-8
    MC = 1      # Mohr-Coulomb failure surface
    CAP = 2     # cap, Eq. 27
    CORNER = 3  # principal stress equality at a triaxial corner


class Shear(IntEnum):
    NONE = 0
    CONE = 1
    MC = 2


@dataclass
class HSSurface:
    """Active surface of the return mapping."""
    type: SurfaceType
    corner: PrincipalCorner = PrincipalCorner.NONE  # corner at which the flow is averaged
    w: Optional[np.ndarray] = None                  # cap: q~ = w . sigma (Eq. 28)
    # constant during a step, set by _prepare_surfaces:
    n: Optional[np.ndarray] = None                  # shear and corner: flow direction
    Dn: Optional[np.ndarray] = None                 # shear and corner: D n
    DN: Optional[np.ndarray] = None                 # cap: D N, with n = N sigma + 2a/3 [1 1 1]


@dataclass
class HSStep:
    """Quantities that are frozen during one (sub-)step (Euler explicit, Sec. 3)."""
    Eur: float
    Ei: float
    D: np.ndarray       # elasticity in principal stress space (3x3)
    sin_psi_m: float    # flow of the shear surfaces, Eqs. 11 and 38
    H: float            # cap hardening modulus, Eqs. 32 and 35
    gamma_p0: float
    p_c0: float


@dataclass
class HSIterate:
    """Stress and hardening state of the return mapping for given plastic multipliers."""
    s: np.ndarray
    gamma_p: float
    p_c: float
    A: np.ndarray       # system matrix I + sum_{cap k} dLambda_k D N_k
    N: np.ndarray       # flow directions dg/dsigma as columns (3 x n_active)
    f: np.ndarray       # yield function values


class HardeningSoil:
    """
    Hardening Soil model with shear hardening, Mohr-Coulomb failure and an elliptic cap.
    """

    PROPS_ORDER = ["E50_ref", "Eur_ref", "m", "c", "phi", "psi", "p_ref", "Rf", "nu", "M_cap",
                   "K_ratio", "e0", "e_cv"]

    def __init__(self, params):
        """
        :param params: dict with the material parameters, see the module docstring
        """
        self._check_properties(params)
        self._set_parameters(params)

        # state
        self.sigma = None
        self.gamma_p = 0.0
        self.p_c = 0.0
        self.at_failure = False
        self.eps_v = 0.0
        self.ddsdde = None
        self.initialized = False

    @classmethod
    def from_props(cls, props):
        """Creates the model from a property list in the order of the C PROPS array."""
        return cls(dict(zip(cls.PROPS_ORDER, props)))

    # ------------------------------------------------------------------
    # Parameters and initial state
    # ------------------------------------------------------------------
    @staticmethod
    def _check_properties(params):
        errors = []
        missing = [key for key in HardeningSoil.PROPS_ORDER[:11] if key not in params]
        if missing:
            raise ValueError(f"Hardening Soil requires the properties {missing}.")

        if params["E50_ref"] <= 0.0:
            errors.append("E50_ref must be positive.")
        if params["Eur_ref"] <= 0.0:
            errors.append("Eur_ref must be positive.")
        if params["c"] < 0.0:
            errors.append("cohesion must be non-negative.")
        if params["phi"] <= 0.0 or params["phi"] >= 90.0:
            errors.append("phi in (0, 90).")
        if params["psi"] > params["phi"]:
            errors.append("psi must not exceed phi.")
        if params["p_ref"] <= 0.0:
            errors.append("p_ref must be positive.")
        if params["Rf"] <= 0.0 or params["Rf"] >= 1.0:
            errors.append("Rf in (0, 1).")
        elif params["Eur_ref"] <= 2.0 * params["E50_ref"] / (2.0 - params["Rf"]):
            errors.append("Eur_ref must exceed Ei_ref = 2 E50_ref / (2 - Rf).")
        if params["nu"] < 0.0 or params["nu"] >= 0.5:
            errors.append("nu in [0, 0.5).")
        if params["M_cap"] < 0.0:
            errors.append("M_cap must be non-negative (0 disables the cap).")
        if params["M_cap"] > 0.0 and params["K_ratio"] <= 1.0:
            errors.append("K_ratio must be > 1.")
        if params.get("e_cv", 0.0) > 0.0 and params.get("e0", 0.0) <= 0.0:
            errors.append("e0 must be positive when the dilatancy cut-off is used.")
        if errors:
            raise ValueError("Hardening Soil: " + " ".join(errors))

    def _set_parameters(self, params):
        self.E50_ref = params["E50_ref"]
        self.Eur_ref = params["Eur_ref"]
        self.m = params["m"]
        self.c = params["c"]
        self.phi = np.radians(params["phi"])
        self.psi = np.radians(params["psi"])
        self.p_ref = params["p_ref"]
        self.Rf = params["Rf"]
        self.nu = params["nu"]
        self.M_cap = params["M_cap"]
        self.K_ratio = params["K_ratio"]
        self.e0 = params.get("e0", 0.0)
        self.e_cv = params.get("e_cv", 0.0)

        self.sin_phi = np.sin(self.phi)
        self.cos_phi = np.cos(self.phi)
        self.sin_psi = np.sin(self.psi)
        self.sin_phi_cv = FlowRules.rowe_critical_state_sin_phi(self.sin_phi, self.sin_psi)
        self.a = self.c * self.cos_phi / self.sin_phi
        self.k_f = YieldFunctions.mohr_coulomb_failure_deviator_factor(self.sin_phi)
        self.alpha = YieldFunctions.elliptic_cap_shape_factor(self.sin_phi)

        # Hyperbola (Eq. 1) with E50 the secant stiffness at q = qf / 2 (Sec. 2.1).
        self.Ei_ref = HardeningLaws.hyperbolic_initial_stiffness(self.E50_ref, self.Rf)

        self.use_cap = self.M_cap > 0.0
        self.use_cutoff = self.e_cv > 0.0

        # H = Ks Kc / (Ks - Kc) = Ks / (K_ratio - 1), Eq. 32
        self.H_ref = 0.0
        if self.use_cap:
            Ks_ref = ElasticLaws.bulk_modulus(self.Eur_ref, self.nu)
            self.H_ref = HardeningLaws.cap_hardening_modulus(Ks_ref, self.K_ratio)

    def _stiffness_factor(self, sigma):
        """((sigma + a) / (p_ref + a))^m, Eqs. 3-4 and 35."""
        return ElasticLaws.power_law_stiffness_factor(sigma, self.a, self.p_ref, self.m,
                                                      HS_MIN_STRESS_RATIO)

    def _initialise_state(self, stress, gamma_p):
        """
        Hardening parameters such that the stress lies on the shear hardening surface (Eq. 8) and on
        the cap (Eq. 27), i.e. a normally consolidated state.

        :return: gamma_p, p_c
        """
        s, _ = StressUtils.principal_system(stress)

        p = (s[0] + s[1] + s[2]) / 3.0
        w = YieldFunctions.elliptic_cap_weights(self.alpha, PrincipalCorner.NONE)
        q_tilde = YieldFunctions.elliptic_cap_equivalent_deviator(w, s)
        q_term = q_tilde / self.M_cap if self.use_cap else 0.0
        p_c0 = np.sqrt(q_term * q_term + (p + self.a) * (p + self.a)) - self.a
        p_c_min = 1.0e-2 * self.p_ref
        p_c = p_c0 if p_c0 > p_c_min else p_c_min

        qf = self.k_f * (s[2] + self.a)
        if qf > 0.0:
            factor = self._stiffness_factor(s[2])
            q = min(s[0] - s[2], qf)
            gamma_0 = HardeningLaws.hyperbolic_plastic_shear_strain(q, qf / self.Rf, self.Ei_ref * factor,
                                                                    self.Eur_ref * factor)
            if gamma_0 > gamma_p:
                gamma_p = gamma_0
        return gamma_p, p_c

    def set_initial_state(self, sigma0, gamma_p=0.0, p_c=0.0, eps_v=0.0):
        """
        Sets the initial stress and state. When p_c <= 0, gamma_p and p_c are initialised such that
        the initial stress lies on the shear hardening surface and on the cap (normally consolidated
        state), as on the first call of the C UMAT.

        :param sigma0: initial stress (compression positive, Voigt [xx, yy, zz, xy, yz, xz])
        :param gamma_p: initial plastic shear strain
        :param p_c: initial pre-consolidation stress
        :param eps_v: initial volumetric strain (compression positive), for the void ratio
        """
        self.sigma = np.array(sigma0, dtype=float)
        self.gamma_p = gamma_p
        self.p_c = p_c
        self.at_failure = False
        self.eps_v = eps_v
        if self.p_c <= 0.0:
            self.gamma_p, self.p_c = self._initialise_state(self.sigma, self.gamma_p)

        factor = self._stiffness_factor(StressUtils.min_principal_stress(self.sigma))
        self.ddsdde = ElasticLaws.elastic_stiffness_matrix_3d(self.Eur_ref * factor, self.nu)
        self.initialized = True

    # ------------------------------------------------------------------
    # Stress-dependent dilatancy
    # ------------------------------------------------------------------
    def _sin_psi_mobilised(self, s, void_ratio):
        """
        Mobilised dilatancy angle from Rowe's stress-dilatancy theory (Eqs. 11-12), including the
        dilatancy cut-off of Eq. 38. Contractant (negative) for phi_m < phi_cv.
        """
        if self.use_cutoff and void_ratio >= self.e_cv:
            return 0.0

        sin_phi_m = YieldFunctions.mohr_coulomb_mobilised_sin_phi(s, self.a, self.sin_phi)
        sin_phi_m = min(max(sin_phi_m, 0.0), self.sin_phi)
        return FlowRules.rowe_mobilised_sin_psi(sin_phi_m, self.sin_phi_cv)

    # ------------------------------------------------------------------
    # Yield surfaces and plastic potentials (principal stresses)
    # ------------------------------------------------------------------
    @staticmethod
    def _corner_function(corner, s):
        """Principal stress equality at a triaxial corner, s2 - s3 = 0 or s1 - s2 = 0."""
        if corner == PrincipalCorner.S2_EQ_S3:
            return s[1] - s[2], np.array([0.0, 1.0, -1.0])
        return s[0] - s[1], np.array([1.0, -1.0, 0.0])

    def _surface_function(self, st, sf, s, gamma_p, p_c):
        """
        Yield function of an active surface (cone: Eq. 8 scaled by (qa - q); MC: Eq. 2; cap: Eqs. 27-29).

        :return: f, df/ds (3), df/dgamma_p, df/dp_c
        """
        if sf.type == SurfaceType.CONE:
            f, df_ds, df_dgamma = YieldFunctions.hyperbolic_shear_yield_function(
                s, gamma_p, st.Ei, st.Eur, self.k_f / self.Rf, self.a)
            return f, df_ds, df_dgamma, 0.0
        if sf.type == SurfaceType.MC:
            f, df_ds = YieldFunctions.mohr_coulomb_principal_function(s, 0, 2, self.sin_phi,
                                                                     self.c * self.cos_phi)
            return f, df_ds, 0.0, 0.0
        if sf.type == SurfaceType.CAP:
            f, df_ds, df_dpc = YieldFunctions.elliptic_cap_yield_function(s, sf.w, p_c, self.M_cap, self.a)
            return f, df_ds, 0.0, df_dpc
        f, df_ds = self._corner_function(sf.corner, s)
        return f, df_ds, 0.0, 0.0

    @staticmethod
    def _shear_flow(st, sf):
        """
        Flow direction of the shear surfaces (Mohr-Coulomb potential g13 with the mobilised dilatancy,
        Eqs. 14-15, averaged with g12 or g23 at a triaxial corner) and of the corner equality.
        """
        if sf.type == SurfaceType.CORNER:
            if sf.corner == PrincipalCorner.S2_EQ_S3:
                return np.array([0.0, 0.5, -0.5])
            return np.array([0.5, -0.5, 0.0])

        sin_psi = st.sin_psi_m
        n = YieldFunctions.mohr_coulomb_principal_gradient(0, 2, sin_psi)
        if sf.corner != PrincipalCorner.NONE:
            if sf.corner == PrincipalCorner.S2_EQ_S3:
                n_corner = YieldFunctions.mohr_coulomb_principal_gradient(0, 1, sin_psi)
            else:
                n_corner = YieldFunctions.mohr_coulomb_principal_gradient(1, 2, sin_psi)
            n = 0.5 * (n + n_corner)
        return n

    def _prepare_surfaces(self, st, surf):
        """
        Quantities of the active surfaces that are constant during a step: the flow direction n and D n
        of the shear and corner surfaces, and D N of the cap, whose associated flow
        n = N sigma + 2a/3 [1 1 1] with N = 2/M^2 w w^T + 2/9 [1 1 1][1 1 1]^T is affine in sigma.
        """
        for sf in surf:
            if sf.type == SurfaceType.CAP:
                N = 2.0 / self.M_cap ** 2 * np.outer(sf.w, sf.w) + 2.0 / 9.0 * np.ones((3, 3))
                sf.DN = st.D @ N
            else:
                sf.n = self._shear_flow(st, sf)
                sf.Dn = st.D @ sf.n

    # ------------------------------------------------------------------
    # Return mapping in principal stress space
    # ------------------------------------------------------------------
    def _evaluate_iterate(self, st, surf, s_tr, dlambda):
        """
        Stress and hardening variables for given plastic multipliers of the active surfaces
        (Eqs. 19, 21, 26 and 33):
          sigma   = sigma_tr - sum_k dLambda_k D n_k
          gamma_p = gamma_p0 + sum_{shear k} dLambda_k
          p_c     = p_c0 + 2 H sum_{cap k} dLambda_k (p + a)
        The associated cap flow is affine in sigma, so the stress follows from the 3x3 linear system
          (I + sum_{cap k} dLambda_k D N_k) sigma
              = sigma_tr - sum_{other k} dLambda_k D n_k - sum_{cap k} dLambda_k 2a/3 D [1 1 1]

        :return: HSIterate, or None if the system matrix is singular
        """
        A = np.eye(3)
        rhs = np.array(s_tr, dtype=float)
        gamma_p = st.gamma_p0
        dlambda_cap = 0.0

        for sf, dl in zip(surf, dlambda):
            if sf.type == SurfaceType.CAP:
                A += dl * sf.DN
                rhs -= dl * 2.0 * self.a / 3.0 * st.D.sum(axis=1)
                dlambda_cap += dl
                continue
            rhs -= dl * sf.Dn
            if sf.type != SurfaceType.CORNER:
                gamma_p += dl

        try:
            s = np.linalg.solve(A, rhs)
        except np.linalg.LinAlgError:
            return None

        p_c = st.p_c0 + 2.0 * st.H * dlambda_cap * (s.sum() / 3.0 + self.a)

        # yield functions, and the flow directions (the gradient for the associated cap flow)
        f = np.zeros(len(surf))
        N = np.zeros((3, len(surf)))
        for k, sf in enumerate(surf):
            f[k], df_ds, _, _ = self._surface_function(st, sf, s, gamma_p, p_c)
            N[:, k] = df_ds if sf.type == SurfaceType.CAP else sf.n
        return HSIterate(s=s, gamma_p=gamma_p, p_c=p_c, A=A, N=N, f=f)

    @staticmethod
    def _merit(f, scale):
        return np.sum((f / scale) ** 2)

    @staticmethod
    def _is_converged(f, scale, tol):
        return np.all(np.abs(f) <= tol * scale)

    def _solve_active_set(self, st, surf, s_tr, scale):
        """
        Newton iteration with a backtracking line search on the plastic multipliers of a fixed set of
        active surfaces, such that all active yield functions vanish (Eqs. 25-26 and 34).

        :return: (dlambda, HSIterate) on convergence, otherwise None
        """
        self._prepare_surfaces(st, surf)
        is_cap = np.array([sf.type == SurfaceType.CAP for sf in surf])
        is_shear = np.array([sf.type in (SurfaceType.CONE, SurfaceType.MC) for sf in surf])

        dlambda = np.zeros(len(surf))
        it = self._evaluate_iterate(st, surf, s_tr, dlambda)
        if it is None:
            return None
        merit = self._merit(it.f, scale)

        for _ in range(HS_MAX_LOCAL_ITER):
            if self._is_converged(it.f, scale, HS_YIELD_TOL):
                return dlambda, it

            # derivatives of the stress and of p_c w.r.t. the multipliers (columns l):
            #   d sigma / d dLambda_l = -(I + sum_{cap k} dLambda_k D N_k)^-1 D n_l
            #   d p_c / d dLambda_l   = 2 H ((p + a) [l is cap] + sum_{cap k} dLambda_k dp / d dLambda_l)
            ds_dl = -np.linalg.solve(it.A, st.D @ it.N)
            dpc_dl = 2.0 * st.H * (is_cap * (it.s.sum() / 3.0 + self.a)
                                   + dlambda[is_cap].sum() * ds_dl.sum(axis=0) / 3.0)

            # Jacobian d f_k / d dLambda_l, with d gamma_p / d dLambda_l = 1 for the shear surfaces
            derivatives = [self._surface_function(st, sf, it.s, it.gamma_p, it.p_c)[1:] for sf in surf]
            df_ds = np.array([d[0] for d in derivatives])
            df_dgamma = np.array([d[1] for d in derivatives])
            df_dpc = np.array([d[2] for d in derivatives])
            J = df_ds @ ds_dl + np.outer(df_dpc, dpc_dl) + np.outer(df_dgamma, is_shear)
            try:
                step = np.linalg.solve(J, -it.f)
            except np.linalg.LinAlgError:
                return None

            # backtracking line search on the scaled residual
            accepted = False
            t = 1.0
            for _ls in range(HS_MAX_LINE_SEARCH):
                lambda_try = dlambda + t * step
                trial = self._evaluate_iterate(st, surf, s_tr, lambda_try)
                if trial is not None:
                    merit_try = self._merit(trial.f, scale)
                    if merit_try < merit:
                        accepted = True
                        merit = merit_try
                        it = trial
                        dlambda = lambda_try
                        break
                t *= 0.5
            if not accepted:
                break  # no further decrease: stagnation at round-off level

        if self._is_converged(it.f, scale, HS_YIELD_TOL_ACCEPT):
            return dlambda, it
        return None

    def _return_mapping(self, st, s_tr):
        """
        Return mapping of the principal trial stress s_tr with an active set strategy over the shear
        hardening surface / Mohr-Coulomb surface, the cap and the triaxial corners (see the C file).

        :return: None on failure, otherwise (s, gamma_p, p_c, plastic, on_failure, active surfaces)
        """
        w_ordered = YieldFunctions.elliptic_cap_weights(self.alpha, PrincipalCorner.NONE)
        qa_factor = self.k_f / self.Rf  # qa = qa_factor (sigma_3 + a), Eqs. 2, 23
        mc_cohesion_term = self.c * self.cos_phi

        def cone_function(s, gamma_p):
            return YieldFunctions.hyperbolic_shear_yield_function(s, gamma_p, st.Ei, st.Eur, qa_factor,
                                                                  self.a)[0]

        def mc_function(s):
            return YieldFunctions.mohr_coulomb_principal_function(s, 0, 2, self.sin_phi, mc_cohesion_term)[0]

        def cap_function(s, p_c):
            return YieldFunctions.elliptic_cap_yield_function(s, w_ordered, p_c, self.M_cap, self.a)[0]

        # fixed residual scales for the convergence checks
        stress_scale = max(abs(s_tr[0]), abs(s_tr[2])) + self.a + 1.0e-3 * self.p_ref
        scale_cone = qa_factor * stress_scale
        scale_mc = stress_scale
        scale_cap = (stress_scale + abs(st.p_c0) + self.a) * (stress_scale + abs(st.p_c0) + self.a)
        order_tol = 1.0e-10 * stress_scale

        shear = Shear.NONE
        corner = PrincipalCorner.NONE
        cap = False

        if cone_function(s_tr, st.gamma_p0) > HS_YIELD_TOL * scale_cone:
            shear = Shear.CONE
        elif mc_function(s_tr) > HS_YIELD_TOL * scale_mc:
            shear = Shear.MC
        if self.use_cap and cap_function(s_tr, st.p_c0) > HS_YIELD_TOL * scale_cap:
            cap = True

        for _pass in range(HS_MAX_ACTIVE_SET_ITER):
            surf = []
            scale = []
            if shear != Shear.NONE:
                surf.append(HSSurface(SurfaceType.CONE if shear == Shear.CONE else SurfaceType.MC, corner))
                scale.append(scale_cone if shear == Shear.CONE else scale_mc)
            if cap:
                surf.append(HSSurface(SurfaceType.CAP, corner,
                                      YieldFunctions.elliptic_cap_weights(self.alpha, corner)))
                scale.append(scale_cap)
            if corner != PrincipalCorner.NONE and surf:
                surf.append(HSSurface(SurfaceType.CORNER, corner))
                scale.append(stress_scale)
            n_act = len(surf)
            scale = np.array(scale)

            result = self._solve_active_set(st, surf, s_tr, scale)
            if result is None:
                return None
            dlambda, it = result

            # 1. surfaces with a negative multiplier are not active
            changed = False
            for k in range(n_act):
                if surf[k].type == SurfaceType.CORNER or dlambda[k] >= 0.0:
                    continue
                changed = True
                if surf[k].type == SurfaceType.CAP:
                    cap = False
                else:
                    shear = Shear.NONE
            if changed:
                continue

            # 2. at a corner, the antisymmetric plastic strain mu must be carried by the pairs
            if corner != PrincipalCorner.NONE and n_act > 0:
                compression = corner == PrincipalCorner.S2_EQ_S3
                capacity = 0.0
                mu = 0.0
                for k in range(n_act):
                    if surf[k].type == SurfaceType.CORNER:
                        mu = dlambda[k]
                    elif surf[k].type == SurfaceType.CAP:
                        q_tilde = YieldFunctions.elliptic_cap_equivalent_deviator(surf[k].w, it.s)
                        capacity += (2.0 * abs(q_tilde) / (self.M_cap * self.M_cap)
                                     * ((2.0 * self.alpha - 1.0) if compression else (2.0 - self.alpha))
                                     * dlambda[k])
                    else:
                        sin_psi = st.sin_psi_m
                        capacity += 0.5 * ((1.0 + sin_psi) if compression else (1.0 - sin_psi)) * dlambda[k]
                if abs(mu) > capacity * (1.0 + 1.0e-8) + SMALL_VALUE:
                    return None

            # 3. principal stress order lost: return to the triaxial corner
            if n_act > 0 and corner == PrincipalCorner.NONE:
                if it.s[1] < it.s[2] - order_tol:
                    corner = PrincipalCorner.S2_EQ_S3
                    changed = True
                elif it.s[0] < it.s[1] - order_tol:
                    corner = PrincipalCorner.S1_EQ_S2
                    changed = True

            # 4. failure criterion q <= qf, otherwise return to the Mohr-Coulomb surface (Sec. 3)
            if shear == Shear.CONE and mc_function(it.s) > HS_YIELD_TOL * scale_mc:
                shear = Shear.MC
                changed = True

            # 5. surfaces that are violated by the returned stress
            if shear == Shear.NONE:
                if cone_function(it.s, it.gamma_p) > HS_YIELD_TOL * scale_cone:
                    shear = Shear.CONE
                    changed = True
                elif mc_function(it.s) > HS_YIELD_TOL * scale_mc:
                    shear = Shear.MC
                    changed = True
            if self.use_cap and not cap and cap_function(it.s, it.p_c) > HS_YIELD_TOL * scale_cap:
                cap = True
                changed = True

            if not changed:
                return it.s, it.gamma_p, it.p_c, n_act > 0, shear == Shear.MC, surf
        return None

    # ------------------------------------------------------------------
    # Single step and integration with automatic sub-stepping
    # ------------------------------------------------------------------
    def _apex_return(self, st, s_tr):
        """
        Return to the apex of the Mohr-Coulomb pyramid, sigma_1 = sigma_2 = sigma_3 = -a, for trial
        stresses in tension that cannot be returned to its faces or edges.

        :return: s, gamma_p, p_c
        """
        G = ElasticLaws.shear_modulus(st.Eur, self.nu)
        s = np.full(3, -self.a)
        gamma_p = st.gamma_p0 + (s_tr[0] - s_tr[2]) / (2.0 * G)
        return s, gamma_p, st.p_c0

    def _koiter_tangent(self, st, active, s, gamma_p, p_c, Q):
        """
        Elasto-plastic (continuum) tangent for the given active surfaces (Koiter):
          D_ep = De - De N (G^T De N + H)^-1 G^T De

        :return: (ok, ddsdde); ddsdde = De if the Koiter matrix is singular
        """
        De = ElasticLaws.elastic_stiffness_matrix_3d(st.Eur, self.nu)
        ddsdde = De.copy()
        n_active = len(active)
        if n_active == 0:
            return True, ddsdde

        # flow directions (rows of N6) and yield function gradients (rows of G6) in Voigt notation
        N6, G6 = np.zeros((n_active, 6)), np.zeros((n_active, 6))
        df_dgamma, df_dpc = np.zeros(n_active), np.zeros(n_active)
        for k, sf in enumerate(active):
            _, g, df_dgamma[k], df_dpc[k] = self._surface_function(st, sf, s, gamma_p, p_c)
            n = g if sf.type == SurfaceType.CAP else self._shear_flow(st, sf)  # associated cap flow, Eq. 30

            g = g.copy()
            if sf.type != SurfaceType.CORNER and sf.corner == PrincipalCorner.S2_EQ_S3:
                g[1] = g[2] = 0.5 * (g[1] + g[2])
            elif sf.type != SurfaceType.CORNER and sf.corner == PrincipalCorner.S1_EQ_S2:
                g[0] = g[1] = 0.5 * (g[0] + g[1])

            N6[k] = StrainUtils.strain_from_principal_system(n, Q)
            G6[k] = StrainUtils.strain_from_principal_system(g, Q)
        De_N = N6 @ De  # rows De n_k (De is symmetric)
        De_G = G6 @ De

        # M = G^T De N + H
        is_cap = np.array([sf.type == SurfaceType.CAP for sf in active])
        is_shear = np.array([sf.type in (SurfaceType.CONE, SurfaceType.MC) for sf in active])
        p = s.sum() / 3.0
        M =(G6 @ De_N.T - np.outer(df_dpc, is_cap * 2.0 * st.H * (p + self.a))
             - np.outer(df_dgamma, is_shear))
        try:
            M_inv = np.linalg.inv(M)
        except np.linalg.LinAlgError:
            return False, ddsdde

        ddsdde -= De_N.T @ M_inv @ De_G

        # shear in the plane of the two equal principal stresses at a corner: no stiffness
        for sf in active:
            if sf.type != SurfaceType.CORNER:
                continue
            a = 1 if sf.corner == PrincipalCorner.S2_EQ_S3 else 0
            dyad = np.outer(Q[:, a], Q[:, a + 1])
            S = StressUtils.matrix_to_voigt(0.5 * (dyad + dyad.T))
            G = ElasticLaws.shear_modulus(st.Eur, self.nu)
            ddsdde -= 4.0 * G * np.outer(S, S)
        return True, ddsdde

    def _elastoplastic_tangent(self, st, active, s, gamma_p, p_c, Q):
        """
        Koiter tangent, with a fraction HS_CORNER_STIFFNESS_FRACTION of the stiffness of the corner
        mode kept at a triaxial corner.
        """
        smooth = [sf for sf in active if sf.type != SurfaceType.CORNER]

        ok, ddsdde = self._koiter_tangent(st, active, s, gamma_p, p_c, Q)
        if not ok or len(smooth) == len(active):
            if len(smooth) != len(active):
                ddsdde = self._koiter_tangent(st, smooth, s, gamma_p, p_c, Q)[1]
            return ddsdde

        D_smooth = self._koiter_tangent(st, smooth, s, gamma_p, p_c, Q)[1]
        return (1.0 - HS_CORNER_STIFFNESS_FRACTION) * ddsdde + HS_CORNER_STIFFNESS_FRACTION * D_smooth

    def _single_step(self, stress, gamma_p, p_c, deps, void_ratio, compute_tangent):
        """
        One (sub-)step: elastic predictor (Eq. 20) and return mapping on the principal trial stresses.

        :return: None on failure, otherwise (stress, gamma_p, p_c, at_failure, ddsdde or None)
        """
        # stiffness and dilatancy at the start of the step (Euler explicit, Sec. 3)
        s0, _ = StressUtils.principal_system(stress)
        factor = self._stiffness_factor(s0[2])
        Eur = self.Eur_ref * factor
        st = HSStep(
            Eur=Eur,
            Ei=self.Ei_ref * factor,
            D=ElasticLaws.elastic_stiffness_matrix_principal(Eur, self.nu),
            sin_psi_m=self._sin_psi_mobilised(s0, void_ratio),
            # cap hardening modulus (Eq. 32), stress dependent through p_c as in Eq. 35
            H=self.H_ref * self._stiffness_factor(p_c) if self.use_cap else 0.0,
            gamma_p0=gamma_p,
            p_c0=p_c,
        )

        # elastic predictor in the global frame, Eq. 20
        Ce = ElasticLaws.elastic_stiffness_matrix_3d(st.Eur, self.nu)
        stress_trial = stress + Ce @ deps

        # return mapping in principal stress space
        s_tr, Q = StressUtils.principal_system(stress_trial)
        result = self._return_mapping(st, s_tr)
        ok = result is not None
        plastic, on_failure, active = False, False, []
        s = None
        if ok:
            s, gamma_p, p_c, plastic, on_failure, active = result

        # a trial stress in tension beyond the apex has no valid return to the faces or edges
        at_apex = False
        apex_tol = 1.0e-10 * (abs(s_tr[0]) + self.a + self.p_ref)
        if not ok or (plastic and s[2] + self.a < -apex_tol):
            if s_tr[2] + self.a >= 0.0:
                return None
            s, gamma_p, p_c = self._apex_return(st, s_tr)
            plastic = True
            on_failure = True
            active = []
            at_apex = True

        new_stress = StressUtils.stress_from_principal_system(s, Q) if plastic else stress_trial

        ddsdde = None
        if compute_tangent:
            if HS_USE_ELASTOPLASTIC_TANGENT:
                ddsdde = self._elastoplastic_tangent(st, active if plastic else [], s, gamma_p, p_c, Q)
            else:
                ddsdde = Ce.copy()
            # at the apex the stress does not change under further loading
            if at_apex:
                ddsdde *= HS_CORNER_STIFFNESS_FRACTION
        return new_stress, gamma_p, p_c, on_failure, ddsdde

    def _integrate(self, stress, gamma_p, p_c, at_failure, dstrain, eps_v0):
        """
        Integrates the strain increment with sub-steps of a fixed size (HS_SUBSTEP_STRESS_FRACTION of
        sigma_3 + a in terms of the elastic stress increment), the first sub-step taking the remaining
        fraction. The number of sub-steps is doubled when a sub-step fails.

        :return: None on failure, otherwise (stress, gamma_p, p_c, at_failure, ddsdde)
        """
        # (real) number of sub-steps from the size of the elastic stress increment
        sigma_3 = StressUtils.min_principal_stress(stress)
        Eur = self.Eur_ref * self._stiffness_factor(sigma_3)
        delta_sigma = ElasticLaws.elastic_stiffness_matrix_3d(Eur, self.nu) @ dstrain

        increment = max(StressUtils.q(delta_sigma), abs(StressUtils.p(delta_sigma)))
        stress_scale = max(sigma_3 + self.a, HS_MIN_STRESS_RATIO * (self.p_ref + self.a))
        n_real = increment / (HS_SUBSTEP_STRESS_FRACTION * stress_scale)
        if not n_real >= 1.0:
            n_real = 1.0
        if n_real > HS_MAX_SUBSTEPS:
            n_real = float(HS_MAX_SUBSTEPS)

        deps_v = StrainUtils.volumetric_strain(dstrain)

        for _trial in range(HS_MAX_SUBSTEP_TRIALS):
            sigma = np.array(stress, dtype=float)
            gp, pc, fail_flag = gamma_p, p_c, at_failure
            ddsdde = None

            # n_full sub-steps of size 1 / n_real, preceded by the remainder
            n_full = int(np.floor(n_real))
            remainder = n_real - n_full
            n_sub = n_full + (1 if remainder > 0.0 else 0)
            fraction_done = 0.0

            ok = True
            for step in range(n_sub):
                fraction = (remainder if (step == 0 and remainder > 0.0) else 1.0) / n_real
                deps_step = dstrain * fraction

                # void ratio at the start of the sub-step, Eq. 39 (eps_v = 0 at e = e0)
                eps_v = eps_v0 + deps_v * fraction_done
                void_ratio = StrainUtils.void_ratio(self.e0, eps_v)
                result = self._single_step(sigma, gp, pc, deps_step, void_ratio, step == n_sub - 1)
                if result is None:
                    ok = False
                    break
                sigma, gp, pc, fail_flag, ddsdde = result
                fraction_done += fraction

            if ok:
                return sigma, gp, pc, fail_flag, ddsdde

            if n_real >= HS_MAX_SUBSTEPS:
                break
            n_real = min(2.0 * n_real, float(HS_MAX_SUBSTEPS))  # refine, retry
        return None

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    def integrate(self, deps):
        """
        Integrates the stress for a strain increment (compression positive, engineering shear strains).

        On success the stress, the state variables and the tangent (ddsdde) are updated. On failure the
        stress and state are kept and ddsdde is the elastic stiffness at the current stress.

        :param deps: strain increment, Voigt [xx, yy, zz, xy, yz, xz]
        :return: True if the integration converged
        """
        if not self.initialized:
            raise ValueError("Initial state not set. Call set_initial_state first.")

        deps = np.array(deps, dtype=float)
        result = self._integrate(self.sigma, self.gamma_p, self.p_c, self.at_failure, deps, self.eps_v)
        if result is None:
            factor = self._stiffness_factor(StressUtils.min_principal_stress(self.sigma))
            self.ddsdde = ElasticLaws.elastic_stiffness_matrix_3d(self.Eur_ref * factor, self.nu)
            return False

        self.sigma, self.gamma_p, self.p_c, self.at_failure, self.ddsdde = result
        self.eps_v += StrainUtils.volumetric_strain(deps)
        return True

    @property
    def state_variables(self):
        """State variables in the order of the C STATEV array: [gamma_p, p_c, at_failure]."""
        return np.array([self.gamma_p, self.p_c, float(self.at_failure)])

    @state_variables.setter
    def state_variables(self, statev):
        self.gamma_p = statev[0]
        self.p_c = statev[1]
        self.at_failure = statev[2] > 0.5
