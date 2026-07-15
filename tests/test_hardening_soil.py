"""
Tests for the Hardening Soil UMAT (c_models/hardening_soil/hardening_soil.c).

The model is a cone (shear) hardening + Mohr-Coulomb failure formulation
(the volumetric cap is present in the source but disabled by default).

Conventions (compression positive, engineering shear):
    PROPS  = [E50_ref, Eur_ref, m, c, phi_deg, psi_deg, p_ref, Rf, nu, M_cap, K_ratio]
    STATEV = [gamma_p, p_p, at_failure]
    Voigt  = [xx, yy, zz, xy, yz, xz]

The tests below are driven directly through the UMAT (Utils.run_c_umat) using
robust, deterministic strain paths so they do not depend on the fragile
fixed-point stress-control driver.
"""

import os
import sys
import shutil
import subprocess

import numpy as np
import pytest

from tests.utils import Utils

# --- material parameters used throughout -------------------------------------
E50_REF, EUR_REF, M, C, PHI_DEG, PSI_DEG = 30000.0, 90000.0, 0.55, 0.0, 42.0, 16.0
P_REF, RF, NU, M_CAP, K_RATIO = 100.0, 0.9, 0.25, 1.0, 1.84
PROPS = [E50_REF, EUR_REF, M, C, PHI_DEG, PSI_DEG, P_REF, RF, NU, M_CAP, K_RATIO]

PHI = np.radians(PHI_DEG)
EI_REF = 2.0 * E50_REF / (2.0 - RF)  # matches the C UMAT
P_T = 0.0                             # c = 0 -> no tensile intercept


# --------------------------------------------------------------------------- #
# Build / locate the shared library                                            #
# --------------------------------------------------------------------------- #
def _dll_path():
    ext = "dll" if sys.platform == "win32" else "so"
    return os.path.join(os.getcwd(), "build_C", "lib", f"hardening_soil.{ext}")


@pytest.fixture(scope="module")
def dll():
    path = _dll_path()
    if os.path.exists(path):
        return path

    # Build with gcc (self-contained: only needs globals, utils, hookes_law).
    if shutil.which("gcc") is None:
        pytest.skip("hardening_soil shared library not built and gcc not available")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    root = os.getcwd()
    sources = [
        os.path.join(root, "c_models", "hardening_soil", "hardening_soil.c"),
        os.path.join(root, "c_models", "globals.c"),
        os.path.join(root, "c_models", "utils.c"),
        os.path.join(root, "c_models", "elastic_laws", "hookes_law.c"),
    ]
    cmd = ["gcc", "-O2", "-shared", "-o", path, *sources, "-lm"]
    subprocess.run(cmd, check=True)
    assert os.path.exists(path)
    return path


# --------------------------------------------------------------------------- #
# Helpers                                                                       #
# --------------------------------------------------------------------------- #
def mean_stress(s):
    return (s[0] + s[1] + s[2]) / 3.0


def dev_q(s):
    p = mean_stress(s)
    d = np.array([s[0] - p, s[1] - p, s[2] - p, s[3], s[4], s[5]])
    j2 = 0.5 * (d[0] ** 2 + d[1] ** 2 + d[2] ** 2) + d[3] ** 2 + d[4] ** 2 + d[5] ** 2
    return np.sqrt(3.0 * j2)


def qf_of(sigma3):
    return (2.0 * np.sin(PHI) / (1.0 - np.sin(PHI))) * (sigma3 + P_T)


def cone_F(s, gamma_p):
    """Cone yield function value (should stay <= tol on the returned state)."""
    sigma3 = min(s[0], s[1], s[2])
    factor = ((sigma3 + P_T) / (P_REF + P_T)) ** M
    Ei = EI_REF * factor
    Eur = EUR_REF * factor
    qa = qf_of(sigma3) / RF
    q = dev_q(s)
    if q <= 0.0:
        return -gamma_p
    if q >= 0.99 * qa:
        return 1e10
    term = 2.0 * q / (1.0 - q / qa)
    return term / Ei - 2.0 * q / Eur - gamma_p


def step(dll, stress, statev, strain, dstrain):
    return Utils.run_c_umat(dll, stress, statev.copy(), strain, dstrain, PROPS, 1)


# --------------------------------------------------------------------------- #
# Tests                                                                         #
# --------------------------------------------------------------------------- #
def test_elastic_isotropic_response(dll):
    """A tiny isotropic strain must give the elastic (Eur-based) response."""
    stress = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    statev = np.array([0.0, 0.0, 0.0])
    e = 1e-6
    dstrain = np.array([e, e, e, 0.0, 0.0, 0.0])

    s_new, ddsdde, sv_new = step(dll, stress, statev, np.zeros(6), dstrain)

    # At sigma3 = p_ref the stiffness factor is 1 -> Eur = Eur_ref.
    K = EUR_REF / (3.0 * (1.0 - 2.0 * NU))
    expected = 100.0 + 3.0 * K * e  # each normal stress increases by 3*K*e

    assert np.allclose(s_new[:3], expected, rtol=1e-6)
    assert np.allclose(s_new[3:], 0.0, atol=1e-9)
    # purely elastic -> no hardening
    assert sv_new[0] == pytest.approx(0.0, abs=1e-12)
    assert sv_new[2] == pytest.approx(0.0, abs=1e-12)
    # DDSDDE is the elastic Eur matrix (evaluated at the slightly updated
    # sigma3, hence a loose tolerance on the stress-dependency factor).
    G = EUR_REF / (2.0 * (1.0 + NU))
    lam = K - 2.0 * G / 3.0
    assert ddsdde[0, 0] == pytest.approx(lam + 2.0 * G, rel=1e-2)
    assert ddsdde[3, 3] == pytest.approx(G, rel=1e-2)


def test_shear_yields_immediately(dll):
    """
    In the Hardening Soil model the cone surface starts at q = 0 (gamma_p = 0),
    so deviatoric loading is elasto-plastic: plastic shear strain accumulates
    and the shear stress is softer than the purely elastic G*gamma response.
    """
    stress = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    statev = np.array([0.0, 0.0, 0.0])
    gamma = 2e-3
    dstrain = np.array([0.0, 0.0, 0.0, gamma, 0.0, 0.0])

    s_new, _, sv_new = step(dll, stress, statev, np.zeros(6), dstrain)

    G = EUR_REF / (2.0 * (1.0 + NU))
    # softer than the purely elastic response because of plastic yielding
    assert 0.0 < s_new[3] < G * gamma
    # plastic shear strain accumulated
    assert sv_new[0] > 0.0
    # the returned state does not violate the cone yield surface
    assert cone_F(s_new, sv_new[0]) <= 1e-3 * (dev_q(s_new) + 1.0)


def test_hardening_monotonic_and_yield_consistent(dll):
    """
    Strain-controlled deviatoric compression: the deviator stress must increase
    monotonically (hardening), the return-mapped state must never violate the
    cone yield surface, and q must stay below the Mohr-Coulomb failure value.
    """
    stress = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    statev = np.array([0.0, 0.0, 0.0])
    strain = np.zeros(6)

    # net-compressive, mostly deviatoric increment
    dstrain = np.array([0.0015, -0.0003, -0.0003, 0.0, 0.0, 0.0])

    q_prev = 0.0
    gamma_prev = 0.0
    for _ in range(150):
        stress, _, statev = step(dll, stress, statev, strain, dstrain)
        strain = strain + dstrain

        assert not np.any(np.isnan(stress)), "NaN in stress"

        q = dev_q(stress)
        sigma3 = min(stress[0], stress[1], stress[2])

        # (1) hardening: q never decreases (allow tiny numerical slack)
        assert q >= q_prev - 1e-3 * (abs(q_prev) + 1.0)
        # (2) hardening variable never decreases
        assert statev[0] >= gamma_prev - 1e-12
        # (3) never above the Mohr-Coulomb failure envelope
        assert q <= qf_of(sigma3) * (1.0 + 1e-3)
        # (4) the returned state satisfies the cone yield condition (F_s <= tol)
        F = cone_F(stress, statev[0])
        tol = 1e-3 * (abs(dev_q(stress)) + 1.0)
        assert F <= tol, f"cone yield violated: F_s = {F:.3e}"

        q_prev = q
        gamma_prev = statev[0]

    # meaningful plastic mobilisation occurred
    assert statev[0] > 0.02
    assert dev_q(stress) > 0.0


def test_triaxial_reaches_and_respects_failure(dll):
    """
    Constant-confinement (100 kPa) drained triaxial via a robust Newton control
    on the lateral strain. The deviator stress must mobilise significantly and
    never exceed the Mohr-Coulomb failure deviator.
    """
    cell = 100.0
    stress = np.array([cell, cell, cell, 0.0, 0.0, 0.0])
    statev = np.array([0.0, 0.0, 0.0])
    strain = np.zeros(6)

    n = 120
    d_axial = 0.12 / n
    q_max = 0.0
    for _ in range(n):
        # Newton on the (symmetric) lateral strain so lateral stress == cell.
        e = 0.0
        for _ in range(50):
            ds = np.array([d_axial, e, e, 0.0, 0.0, 0.0])
            s_try, _, _ = step(dll, stress, statev, strain, ds)
            r = s_try[1] - cell
            if abs(r) < 1e-7:
                break
            de = 1e-6
            ds2 = np.array([d_axial, e + de, e + de, 0.0, 0.0, 0.0])
            s2, _, _ = step(dll, stress, statev, strain, ds2)
            deriv = (s2[1] - s_try[1]) / de
            if abs(deriv) < 1e-9:
                break
            e -= r / deriv

        ds = np.array([d_axial, e, e, 0.0, 0.0, 0.0])
        stress, _, statev = step(dll, stress, statev, strain, ds)
        strain = strain + ds

        assert not np.any(np.isnan(stress))
        sigma3 = min(stress[0], stress[1], stress[2])
        q = stress[0] - stress[1]
        # never above failure
        assert q <= qf_of(sigma3) * (1.0 + 5e-3)
        q_max = max(q_max, q)

    # confinement held reasonably and significant strength mobilised
    assert abs(stress[1] - cell) < 0.15 * cell
    assert q_max > 0.4 * qf_of(cell)
    # plastic shear strain accumulated
    assert statev[0] > 0.05




