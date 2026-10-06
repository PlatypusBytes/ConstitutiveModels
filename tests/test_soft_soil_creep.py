"""
Tests for the Soft Soil Creep UMAT (c_models/soft_soil/soft_soil_creep.c).

The model follows Stolle, Vermeer & Bonnier (1999), "Time integration of a constitutive law for
soft clays": volumetric creep on the Modified Cam-Clay equivalent pressure p_eq, a pressure
dependent bulk modulus K = p / A and Mohr-Coulomb failure with zero dilatancy. Equation numbers
refer to that paper; the material properties are those of its Table I.

Conventions (compression positive, engineering shear):
    PROPS  = [A, B, C, nu, tau, M, phi_deg, OCR0(, p_min)]
    STATEV = [eps_c, p_p, p_eq, OCR, at_failure]
    Voigt  = [xx, yy, zz, xy, yz, xz]
    DTIME  = time increment, in the unit of tau (minutes here)

The UMAT interface itself is tension positive (Abaqus convention, like the other models in this
library); the helper `step` converts, so the tests below are written compression positive as the
paper. All tests drive the UMAT directly along strain paths.
"""

import os
import sys
import shutil
import subprocess

import numpy as np
import pytest

from tests.utils import Utils

# --- material parameters of Table I --------------------------------------------
A, B, C, NU, TAU, M = 0.016, 0.090, 0.004, 0.25, 0.70, 1.29
PHI_DEG = 30.0  # M_mc = 6 sin(phi) / (3 - sin(phi)) = 1.20, the failure envelope of the paper
M_MC = 6.0 * np.sin(np.radians(PHI_DEG)) / (3.0 - np.sin(np.radians(PHI_DEG)))


def props(ocr0=1.0, phi_deg=PHI_DEG):
    return [A, B, C, NU, TAU, M, phi_deg, ocr0]


# --------------------------------------------------------------------------- #
# Build / locate the shared library                                            #
# --------------------------------------------------------------------------- #
def _dll_path():
    ext = "dll" if sys.platform == "win32" else "so"
    return os.path.join(os.getcwd(), "build_C", "lib", f"soft_soil_creep.{ext}")


@pytest.fixture(scope="module")
def dll():
    path = _dll_path()
    if os.path.exists(path):
        return path

    # Build with gcc from the model and the shared modules it uses.
    if shutil.which("gcc") is None:
        pytest.skip("soft_soil_creep shared library not built and gcc not available")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    root = os.path.join(os.getcwd(), "c_models")
    sources = [
        os.path.join(root, "soft_soil", "soft_soil_creep.c"),
        os.path.join(root, "globals.c"),
        os.path.join(root, "utils.c"),
        os.path.join(root, "stress_utils.c"),
        os.path.join(root, "strain_utils.c"),
        os.path.join(root, "elastic_laws", "hookes_law.c"),
        os.path.join(root, "elastic_laws", "logarithmic_elasticity.c"),
        os.path.join(root, "yield_surfaces", "modified_cam_clay_surface.c"),
        os.path.join(root, "yield_surfaces", "mohr_coulomb_surface.c"),
        os.path.join(root, "flow_rules", "modified_cam_clay_flow.c"),
        os.path.join(root, "hardening_rules", "exponential_volumetric_hardening.c"),
        os.path.join(root, "creep_laws", "isotache_creep.c"),
        os.path.join(root, "return_mappings", "mohr_coulomb_return_mapping.c"),
    ]
    subprocess.run(["gcc", "-O2", "-shared", "-o", path, *sources, "-lm"], check=True)
    assert os.path.exists(path)
    return path


# --------------------------------------------------------------------------- #
# Helpers                                                                       #
# --------------------------------------------------------------------------- #
def step(dll, stress, statev, dstrain, prm, dt):
    """UMAT call with stresses and strains in the compression-positive convention of the paper."""
    stress_new, ddsdde, statev_new = Utils.run_c_umat(
        dll, -np.asarray(stress, dtype=float), np.array(statev, dtype=float), np.zeros(6),
        -np.asarray(dstrain, dtype=float), prm, dt)
    return -stress_new, ddsdde, statev_new


def p_q(s):
    p = np.mean(s[:3])
    d = np.asarray(s[:3]) - p
    return p, np.sqrt(1.5 * (d @ d) + 3.0 * (np.asarray(s[3:]) @ np.asarray(s[3:])))


def equivalent_pressure(s):
    p, q = p_q(s)
    return p + q * q / (M * M * p)


def elastic_matrix(K, G):
    D = np.zeros((6, 6))
    D[:3, :3] = K - 2.0 * G / 3.0
    D[range(3), range(3)] += 2.0 * G
    D[range(3, 6), range(3, 6)] = G
    return D


def shear_modulus(p):
    return 3.0 * (p / A) * (1.0 - 2.0 * NU) / (2.0 * (1.0 + NU))


def creep_increment(p_eq, p_p, dt):
    """Eq. 5."""
    return C * np.log1p(dt / TAU * (p_eq / p_p) ** (B / C))


def bisect(f, lo, hi, n=200):
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def fd_tangent(dll, stress, statev, dstrain, prm, dt, h=1e-7):
    J = np.zeros((6, 6))
    for j in range(6):
        d_plus, d_minus = dstrain.copy(), dstrain.copy()
        d_plus[j] += h
        d_minus[j] -= h
        J[:, j] = (step(dll, stress, statev, d_plus, prm, dt)[0] -
                   step(dll, stress, statev, d_minus, prm, dt)[0]) / (2.0 * h)
    return J


def stress_tensor(s):
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


def stress_voigt(t):
    return np.array([t[0, 0], t[1, 1], t[2, 2], t[0, 1], t[1, 2], t[0, 2]])


def rotate_stress(s, R):
    return stress_voigt(R @ stress_tensor(s) @ R.T)


def rotate_strain(e, R):
    """Rotate a Voigt strain with engineering shear components."""
    t = stress_tensor(np.array([e[0], e[1], e[2], e[3] / 2, e[4] / 2, e[5] / 2]))
    v = stress_voigt(R @ t @ R.T)
    v[3:] *= 2.0
    return v


def run_crs(dll, path, p0, ocr0, n_steps, total_strain=0.2, total_time=2160.0):
    """
    Constant rate of strain test of the paper (axial direction x), starting from an isotropic
    stress p0: 'oedometer' (only axial strain) or 'undrained' triaxial (isochoric).
    Returns per step: axial strain, stress, state variables.
    """
    stress = np.array([p0, p0, p0, 0.0, 0.0, 0.0])
    statev = np.zeros(5)
    de = total_strain / n_steps
    if path == "oedometer":
        dstrain = np.array([de, 0.0, 0.0, 0.0, 0.0, 0.0])
    else:
        dstrain = np.array([de, -de / 2.0, -de / 2.0, 0.0, 0.0, 0.0])
    history = []
    for i in range(n_steps):
        stress, _, statev = step(dll, stress, statev, dstrain, props(ocr0), total_time / n_steps)
        history.append(((i + 1) * de, stress.copy(), statev.copy()))
    return history


# --------------------------------------------------------------------------- #
# Initial state and elasticity                                                  #
# --------------------------------------------------------------------------- #
def test_initial_preconsolidation_pressure(dll):
    """On the first call (p_p = 0) the pre-consolidation pressure is OCR0 p_eq0 (Eq. 4)."""
    stress = np.array([100.0, 60.0, 60.0, 0.0, 0.0, 0.0])
    ocr0 = 2.0
    _, _, sv = step(dll, stress, np.zeros(5), np.zeros(6), props(ocr0), 1e-6)

    p_eq0 = equivalent_pressure(stress)
    assert sv[0] == pytest.approx(0.0, abs=1e-12)          # no creep in 1e-6 min at OCR 2
    assert sv[1] == pytest.approx(ocr0 * p_eq0, rel=1e-10)
    assert sv[2] == pytest.approx(p_eq0, rel=1e-10)
    assert sv[3] == pytest.approx(ocr0, rel=1e-10)
    assert sv[4] == 0.0


def test_elastic_isotropic_response(dll):
    """
    Heavily over-consolidated: the isotropic response follows the logarithmic elastic law
    K = p / A, integrated exactly: p = p0 exp(deps_v / A).
    """
    p0, e = 100.0, 1e-3
    stress = np.array([p0, p0, p0, 0.0, 0.0, 0.0])
    statev = np.array([0.0, 1e4, 0.0, 0.0, 0.0])
    dstrain = np.array([e, e, e, 0.0, 0.0, 0.0])

    s_new, ddsdde, sv = step(dll, stress, statev, dstrain, props(), 1e-3)

    p = p0 * np.exp(3.0 * e / A)
    assert np.allclose(s_new[:3], p, rtol=1e-10)
    assert np.allclose(s_new[3:], 0.0, atol=1e-12)
    assert sv[0] == pytest.approx(0.0, abs=1e-15)
    assert sv[1] == 1e4
    # tangent: K of the logarithmic law at the end, G at the start of the increment
    assert np.allclose(ddsdde, elastic_matrix(p / A, shear_modulus(p0)), rtol=1e-10)


def test_zero_time_increment_is_elastic(dll):
    """Without time there is no creep, even in a normally consolidated state (Eq. 5)."""
    stress = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    statev = np.array([0.0, 100.0, 0.0, 0.0, 0.0])
    dstrain = np.array([2e-3, 0.0, 0.0, 0.0, 0.0, 0.0])

    s_new, ddsdde, sv = step(dll, stress, statev, dstrain, props(phi_deg=0.0), 0.0)

    G = shear_modulus(100.0)
    s_expected = np.full(6, 0.0)
    s_expected[:3] = 100.0 * np.exp(2e-3 / A)
    s_expected[:3] += 2.0 * G * (np.array([2e-3, 0.0, 0.0]) - 2e-3 / 3.0)
    assert np.allclose(s_new, s_expected, rtol=1e-10, atol=1e-10)
    assert sv[0] == 0.0
    assert sv[1] == 100.0


def test_tangent_matches_finite_differences_when_elastic(dll):
    stress = np.array([100.0, 80.0, 60.0, 5.0, -3.0, 2.0])
    statev = np.array([0.0, 500.0, 0.0, 0.0, 0.0])
    dstrain = np.array([1e-4, -2e-5, 3e-5, 1e-4, 0.0, -5e-5])

    _, ddsdde, sv = step(dll, stress, statev, dstrain, props(), 1e-6)
    J = fd_tangent(dll, stress, statev, dstrain, props(), 1e-6)

    assert sv[0] < 1e-12
    assert np.abs(ddsdde - J).max() <= 1e-7 * np.abs(J).max()


# --------------------------------------------------------------------------- #
# Creep                                                                         #
# --------------------------------------------------------------------------- #
def test_isotropic_creep_at_constant_stress_matches_analytic(dll):
    """
    At constant isotropic stress p = p_p0 the creep strain is eps_c = C ln(1 + t / tau) (Eqs. 2-3).
    The strain increment that keeps the stress constant is the creep increment of Eq. 5, and the
    increments add up exactly for any sequence of time steps.
    """
    p0 = 100.0
    stress = np.array([p0, p0, p0, 0.0, 0.0, 0.0])
    statev = np.array([0.0, p0, 0.0, 0.0, 0.0])
    time = 0.0
    for dt in [0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0]:
        x = creep_increment(p0, statev[1], dt)
        dstrain = np.array([x / 3.0, x / 3.0, x / 3.0, 0.0, 0.0, 0.0])
        stress, _, statev = step(dll, stress, statev, dstrain, props(), dt)
        time += dt

        eps_c = C * np.log1p(time / TAU)
        assert np.allclose(stress[:3], p0, rtol=1e-9)
        assert statev[0] == pytest.approx(eps_c, rel=1e-9)
        assert statev[1] == pytest.approx(p0 * np.exp(eps_c / B), rel=1e-9)


def test_relaxation_at_constant_strain(dll):
    """
    Without strain, creep is compensated by elastic swelling: p = p0 exp(-x / A) with
    x = deps_c(p) of Eq. 5 evaluated at the end of the increment (implicit).
    """
    p0, dt = 100.0, 21.6
    stress = np.array([p0, p0, p0, 0.0, 0.0, 0.0])
    statev = np.array([0.0, p0, 0.0, 0.0, 0.0])

    s_new, _, sv = step(dll, stress, statev, np.zeros(6), props(), dt)

    x = bisect(lambda x: creep_increment(p0 * np.exp(-x / A), p0, dt) - x,
               0.0, creep_increment(p0, p0, dt))
    assert sv[0] == pytest.approx(x, rel=1e-9)
    assert np.allclose(s_new[:3], p0 * np.exp(-x / A), rtol=1e-9)

    # the stress keeps relaxing, ever more slowly
    p_prev, drop_prev = s_new[0], p0 - s_new[0]
    for _ in range(5):
        s_new, _, sv = step(dll, s_new, sv, np.zeros(6), props(), dt)
        drop = p_prev - s_new[0]
        assert 0.0 < drop < drop_prev
        p_prev, drop_prev = s_new[0], drop


def test_tangent_is_eq9(dll):
    """
    DDSDDE is D^ec of Eq. 9 at the stress after the creep update, with the pre-consolidation
    pressure of the start of the increment.
    """
    stress = np.array([100.0, 80.0, 60.0, 5.0, -3.0, 2.0])
    dstrain = np.array([1e-3, -2e-4, 3e-4, 1e-3, 0.0, -5e-4])
    dt = 21.6
    s_new, ddsdde, sv = step(dll, stress, np.zeros(5), dstrain, props(phi_deg=0.0), dt)
    assert sv[0] > 1e-3  # creep active

    p_p = equivalent_pressure(stress)  # OCR0 = 1
    p, q = p_q(s_new)
    p_eq = p + q * q / (M * M * p)
    dev = s_new.copy()
    dev[:3] -= p
    dq_dsigma = 1.5 * dev / q
    dq_dsigma[3:] *= 2.0
    a = (1.0 - q * q / (M * M * p * p)) * np.array([1, 1, 1, 0, 0, 0]) / 3.0 \
        + 2.0 * q / (M * M * p) * dq_dsigma
    x = np.log(dt / TAU) + B / C * np.log(p_eq / p_p)
    c = B / p_eq / (1.0 + np.exp(-x))  # d deps_c / d p_eq, Eq. 5
    D = elastic_matrix(p / A, shear_modulus(np.mean(stress[:3])))
    Da = D @ a
    D_ec = D - c * np.outer(Da, Da) / (1.0 - q * q / (M * M * p * p) + c * a @ Da)

    assert np.allclose(ddsdde, D_ec, rtol=1e-10, atol=1e-10 * np.abs(D).max())


# --------------------------------------------------------------------------- #
# Constant rate of strain tests of the paper                                    #
# --------------------------------------------------------------------------- #
def test_oedometer_reaches_k0_line(dll):
    """
    Fig. 2: CRS oedometer test, OCR0 = 6, p0 = 10 kPa, 20% strain in 2160 min. After the
    recompression range the stress path follows the K0 line with K0 of about 0.62.
    """
    history = run_crs(dll, "oedometer", 10.0, 6.0, 100)

    sigma_v_prev = 0.0

    for t, s, sv in history:
        assert s[0] > sigma_v_prev  # vertical stress increases monotonically
        assert sv[4] == 0.0         # far from failure
        sigma_v_prev = s[0]


    # recompression: no creep in the first steps
    assert history[4][2][0] < 1e-6
    # normally consolidated at the end, on the K0 line
    _, s_end, sv_end = history[-1]
    K0 = s_end[1] / s_end[0]
    assert 0.58 < K0 < 0.66
    assert s_end[2] == pytest.approx(s_end[1], rel=1e-10)
    assert sv_end[0] > 0.1                      # most of the strain is creep
    assert 1.0 < sv_end[3] < 1.5                # constant OCR under constant strain rate
    assert history[-10][2][3] == pytest.approx(sv_end[3], rel=1e-2)


def test_undrained_triaxial_reaches_failure(dll):
    """
    Fig. 3: CRS undrained triaxial test on a normally consolidated clay, p0 = 100 kPa. q increases
    until the Mohr-Coulomb envelope q = M_mc p is reached, and then moves down along it due to the
    creep-induced excess pore pressure (strain softening).
    """
    history = run_crs(dll, "undrained", 100.0, 1.0, 100)
    p_values = np.array([p_q(s)[0] for _, s, _ in history])
    q_values = np.array([p_q(s)[1] for _, s, _ in history])
    strain_values = [t for t, _, _ in history]

    # never above the failure envelope; the lateral stresses stay equal
    for (_, s, _), p, q in zip(history, p_values, q_values):
        assert q <= M_MC * p * (1.0 + 1e-10)
        assert s[1] == pytest.approx(s[2], rel=1e-10)
    # mean effective stress decreases monotonically (positive excess pore pressure)
    assert np.all(np.diff(p_values) < 0.0)

    i_peak = np.argmax(q_values)
    assert 50.0 < q_values[i_peak] < 60.0       # paper: about 55 kPa
    assert q_values[-1] < q_values[i_peak]      # strain softening
    assert history[-1][2][4] == 1.0             # on the failure surface
    assert q_values[-1] == pytest.approx(M_MC * p_values[-1], rel=1e-10)


def test_large_time_steps_match_small_time_steps(dll):
    """The modified procedure is accurate with large time steps (21.6 vs 2.16 min, Figs. 2-3)."""
    for path, p0, ocr0 in [("oedometer", 10.0, 6.0), ("undrained", 100.0, 1.0)]:
        coarse = run_crs(dll, path, p0, ocr0, 100)
        fine = run_crs(dll, path, p0, ocr0, 1000)
        q_coarse = np.array([p_q(s)[1] for _, s, _ in coarse])
        q_fine = np.array([p_q(s)[1] for _, s, _ in fine])[9::10]
        assert np.max(q_coarse) == pytest.approx(np.max(q_fine), rel=0.02)
        assert np.allclose(coarse[-1][1][:3], fine[-1][1][:3], rtol=0.05)


def test_undrained_strength_increases_with_strain_rate(dll):
    """The peak undrained shear strength is rate dependent: a faster test gives a higher peak."""
    peaks = []
    for total_time in [21600.0, 2160.0, 216.0]:
        history = run_crs(dll, "undrained", 100.0, 1.0, 100, total_time=total_time)
        peaks.append(max(p_q(s)[1] for _, s, _ in history))
    assert peaks[0] < peaks[1] < peaks[2]


# --------------------------------------------------------------------------- #
# General 3D states                                                              #
# --------------------------------------------------------------------------- #
def test_rotation_invariance(dll):
    """Isotropic model: rotating the stress and the strain increment rotates the result."""
    stress = np.array([100.0, 70.0, 60.0, 8.0, -4.0, 3.0])
    dstrain = np.array([3e-3, -1e-3, -1.5e-3, 2e-3, 1e-3, -1e-3])
    angles = np.radians([30.0, -20.0, 45.0])
    Rx = np.array([[1, 0, 0], [0, np.cos(angles[0]), -np.sin(angles[0])],
                   [0, np.sin(angles[0]), np.cos(angles[0])]])
    Rz = np.array([[np.cos(angles[2]), -np.sin(angles[2]), 0],
                   [np.sin(angles[2]), np.cos(angles[2]), 0], [0, 0, 1]])
    R = Rz @ Rx

    for dt in [1e-6, 21.6]:  # failure dominated and creep dominated
        s_a, _, sv_a = step(dll, stress, np.zeros(5), dstrain, props(), dt)
        s_b, _, sv_b = step(dll, rotate_stress(stress, R), np.zeros(5), rotate_strain(dstrain, R),
                            props(), dt)
        assert np.allclose(rotate_stress(s_a, R), s_b, rtol=1e-9, atol=1e-9)
        assert np.allclose(sv_a, sv_b, rtol=1e-9, atol=1e-12)


def test_mohr_coulomb_is_respected_on_general_paths(dll):
    """Along a 3D shear path the returned stresses never violate the Mohr-Coulomb criterion."""
    sin_phi = np.sin(np.radians(PHI_DEG))
    stress = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    statev = np.zeros(5)
    dstrain = np.array([2e-3, -1.5e-3, -0.5e-3, 1e-3, -2e-3, 0.5e-3])
    failed = False
    for _ in range(40):
        stress, _, statev = step(dll, stress, statev, dstrain, props(), 1.0)
        s = np.sort(np.linalg.eigvalsh(stress_tensor(stress)))[::-1]
        f13 = 0.5 * (s[0] - s[2]) - 0.5 * (s[0] + s[2]) * sin_phi
        assert f13 <= 1e-9 * s[0]
        failed = failed or statev[4] == 1.0
    assert failed
