"""
Tests for the Matsuoka-Nakai UMAT (c_models/matsuoka_nakai/matsuoka_nakai.c).

The model is linear elastic, perfectly plastic with the Matsuoka-Nakai yield surface
I1 I2 / I3 = 9 + 8 tan^2(phi) of the stresses shifted by a = c cot(phi) (Matsuoka & Nakai, 1974),
in the formulation of Lagioia & Panteghini (2016), and a plastic potential of the same form with
the dilation angle psi.

    PROPS  = [E, nu, c, phi_deg, psi_deg]
    STATEV = [state]  (0 elastic, 1 plastic on the cone, 2 plastic at the apex)
    Voigt  = [xx, yy, zz, xy, yz, xz], engineering shear strains

The first two tests drive the UMAT through the incremental driver in its own (tension positive)
convention. The other tests call the UMAT directly through the helper `step`, which converts to
compression positive stresses and strains, as is usual in soil mechanics.
"""

import os
import sys
import shutil
import subprocess

import cffi
import numpy as np
import pytest

from tests.incr_driver import IncrDriver
from tests.utils import Utils



def test_strain_controlled_compression_triaxial():
    """
    Test the strain controlled triaxial test in compression.

    :return:
    """

    G = 2000.
    nu = 0.33
    E = 2. * G * (1. + nu)

    params = {'E': E, 'poison_ratio': nu, 'cohesion': 0.1, 'friction_angle': 30., 'dilation_angle': 0.0}
    state_vars={'plastic_strain': 0.}

    # Define the original stress vector
    orig_stress_vector = np.array([-1., -1., -1., 0, 0, 0])

    project_dir = os.getcwd()
    # check operating system
    if sys.platform == 'win32':
        extension = 'dll'
    elif sys.platform == 'linux':
        extension = 'so'
    else:
        raise Exception("Unsupported operating system")
    model_loc = os.path.join(project_dir,'build_C','lib','matsuoka_nakai.'+extension)

    const_model_info = {"language": "c",
                        "file_name":model_loc,
                        "properties": list(params.values()),
                        "state_vars": list(state_vars.values())}

    # strain increment per time step in the form of [eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz]

    vertical_strain_increment = -1e-4

    strain_increment = np.array([0, vertical_strain_increment, 0, 0, 0, 0])
    stress_increment = np.zeros_like(strain_increment)

    # solve the problem
    incr_driver = IncrDriver(orig_stress_vector,
                             strain_increment,
                             stress_increment,
                             const_model_info,
                             6,
                             100)

    incr_driver.solve_triaxial_strain_controlled()

    expected_sigma_3 = orig_stress_vector[0]
    phi_rad = np.radians(params['friction_angle'])
    c = params['cohesion']

    expected_first_yield_strain = 1.75020728e-04
    expected_yield_strain_increment = 0.00005

    expected_strains = np.array([[-vertical_strain_increment*nu, vertical_strain_increment, -vertical_strain_increment*nu, 0, 0, 0],
                                 [2*-vertical_strain_increment*nu, 2*vertical_strain_increment, 2*-vertical_strain_increment*nu, 0, 0, 0],
                                 [3*-vertical_strain_increment*nu, 3*vertical_strain_increment, 3*-vertical_strain_increment*nu, 0, 0, 0],
                                 [4*-vertical_strain_increment*nu, 4*vertical_strain_increment, 4*-vertical_strain_increment*nu, 0, 0, 0],
                                 [expected_first_yield_strain, 5*vertical_strain_increment, expected_first_yield_strain, 0, 0, 0],
                                 [expected_first_yield_strain+ expected_yield_strain_increment, 6*vertical_strain_increment, expected_first_yield_strain + expected_yield_strain_increment, 0, 0, 0]])


    # expected vertical_yield_stress  is equal to the mohr coulomb vertical yield stress
    expected_vertical_yield_stress = (expected_sigma_3 * (1 + np.sin(phi_rad)) - 2 * c * np.cos(phi_rad)) / (1 - np.sin(phi_rad))

    expected_stresses = np.array([orig_stress_vector,orig_stress_vector,orig_stress_vector,
                                  orig_stress_vector, orig_stress_vector ,orig_stress_vector], dtype=float)
    expected_stresses[0,1] = expected_stresses[0,1]+ vertical_strain_increment * E
    expected_stresses[1,1] = expected_stresses[1,1] + vertical_strain_increment * E * 2
    expected_stresses[2,1] = expected_stresses[2,1] + vertical_strain_increment * E * 3
    expected_stresses[3, 1] = expected_stresses[3, 1] + vertical_strain_increment * E * 4
    expected_stresses[4, 1] = expected_vertical_yield_stress
    expected_stresses[5, 1] = expected_vertical_yield_stress

    # Check the results
    np.testing.assert_allclose(incr_driver.strains, expected_strains, rtol=1e-6)
    np.testing.assert_allclose(incr_driver.stresses, expected_stresses, rtol=1e-6)


def test_strain_controlled_tension_triaxial():
    """
    Test the strain controlled triaxial test in tension.

    :return:
    """

    G = 2000
    nu = 0.33
    E = 2 * G * (1 + nu)

    params = {'E': E, 'poison_ratio': nu, 'cohesion': 0.1, 'friction_angle': 30., 'dilation_angle': 0.0}
    state_vars={'plastic_strain': 0.}

    # Define the original stress vector
    orig_stress_vector = np.array([-1., -1., -1., 0, 0, 0])

    project_dir = os.getcwd()
    # check operating system
    if sys.platform == 'win32':
        extension = 'dll'
    elif sys.platform == 'linux':
        extension = 'so'
    else:
        raise Exception("Unsupported operating system")
    model_loc = os.path.join(project_dir,'build_C','lib','matsuoka_nakai.'+extension)

    const_model_info = {"language": "c",
                        "file_name":model_loc,
                        "properties": list(params.values()),
                        "state_vars": list(state_vars.values())}

    # strain increment per time step in the form of [eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz]
    vertical_strain_increment = 1e-4
    strain_increment = np.array([0, vertical_strain_increment, 0, 0, 0, 0])
    stress_increment = np.zeros_like(strain_increment)

    # solve the problem
    incr_driver = IncrDriver(orig_stress_vector,
                             strain_increment,
                             stress_increment,
                             const_model_info,
                             3,
                             100)

    incr_driver.solve_triaxial_strain_controlled()


    expected_sigma_3 = orig_stress_vector[0]
    phi_rad = np.radians(params['friction_angle'])
    c = params['cohesion']

    # expected vertical_yield_stress  is equal to the mohr coulomb vertical yield stress
    expected_vertical_yield_stress = (expected_sigma_3 * (1 - np.sin(phi_rad)) + 2 * c * np.cos(phi_rad)) / (1 + np.sin(phi_rad))

    expected_first_yield_strain = -7.50069093e-05
    expected_yield_strain_increment = -0.00005

    expected_strains = np.array([[-vertical_strain_increment*nu, vertical_strain_increment, -vertical_strain_increment*nu, 0, 0, 0],
                                    [expected_first_yield_strain, 2*vertical_strain_increment, expected_first_yield_strain, 0, 0, 0],
                                    [expected_first_yield_strain + expected_yield_strain_increment, 3*vertical_strain_increment, expected_first_yield_strain + expected_yield_strain_increment, 0, 0, 0]])

    expected_stresses = np.array([orig_stress_vector,orig_stress_vector,orig_stress_vector], dtype=float)
    expected_stresses[0,1] = expected_stresses[0,1]+ vertical_strain_increment * E
    expected_stresses[1,1] = expected_vertical_yield_stress
    expected_stresses[2,1] = expected_vertical_yield_stress

    # Check the results
    np.testing.assert_allclose(incr_driver.strains, expected_strains, rtol=1e-6)
    np.testing.assert_allclose(incr_driver.stresses, expected_stresses, rtol=1e-6)


# --------------------------------------------------------------------------- #
# Direct UMAT tests (compression positive)                                      #
# --------------------------------------------------------------------------- #
E_MOD, NU = 20000.0, 0.3
C, PHI_DEG, PSI_DEG = 5.0, 30.0, 10.0


def props_of(c=C, phi=PHI_DEG, psi=PSI_DEG):
    return [E_MOD, NU, c, phi, psi]


def _dll_path():
    ext = "dll" if sys.platform == "win32" else "so"
    return os.path.join(os.getcwd(), "build_C", "lib", f"matsuoka_nakai.{ext}")


@pytest.fixture(scope="module")
def dll():
    path = _dll_path()
    if os.path.exists(path):
        return path

    # Build with gcc from the model and the shared modules it uses.
    if shutil.which("gcc") is None:
        pytest.skip("matsuoka_nakai shared library not built and gcc not available")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    root = os.path.join(os.getcwd(), "c_models")
    sources = [
        os.path.join(root, "matsuoka_nakai", "matsuoka_nakai.c"),
        os.path.join(root, "globals.c"),
        os.path.join(root, "utils.c"),
        os.path.join(root, "stress_utils.c"),
        os.path.join(root, "elastic_laws", "hookes_law.c"),
        os.path.join(root, "yield_surfaces", "matsuoka_nakai_surface.c"),
    ]
    subprocess.run(["gcc", "-O2", "-shared", "-o", path, *sources, "-lm"], check=True)
    assert os.path.exists(path)
    return path


def elastic_stiffness():
    G = E_MOD / (2.0 * (1.0 + NU))
    lam = E_MOD * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] += 2.0 * G
    D[3:, 3:] = G * np.eye(3)
    return D


def step(dll, stress, dstrain, props):
    """
    UMAT call with stresses and strains in the compression-positive convention; the UMAT
    interface itself is tension positive (Abaqus convention). The tangent is the same in both.
    """
    stress_new, ddsdde, statev = Utils.run_c_umat(
        dll, -np.asarray(stress, dtype=float), np.zeros(1), np.zeros(6),
        -np.asarray(dstrain, dtype=float), props, 1)
    return -stress_new, ddsdde, statev


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


def attraction(c, phi_deg):
    """a = c cot(phi), the shift of the stresses (tensile strength under isotropic tension)."""
    return c / np.tan(np.radians(phi_deg))


def principal_stresses(s):
    """Principal stresses sigma_1 >= sigma_2 >= sigma_3 (compression positive)."""
    return np.sort(np.linalg.eigvalsh(stress_tensor(s)))[::-1]


def matsuoka_nakai_criterion(s, c, phi_deg):
    """
    I1 I2 / I3 / (9 + 8 tan^2 phi) - 1 of the principal stresses shifted by a = c cot(phi):
    zero on the Matsuoka-Nakai surface, positive outside.
    """
    sig = principal_stresses(s) + attraction(c, phi_deg)
    i1 = sig.sum()
    i2 = sig[0] * sig[1] + sig[1] * sig[2] + sig[2] * sig[0]
    i3 = sig.prod()
    return i1 * i2 / i3 / (9.0 + 8.0 * np.tan(np.radians(phi_deg)) ** 2) - 1.0


def mean_stress(s):
    return (s[0] + s[1] + s[2]) / 3.0


def deviator_q(s):
    p = mean_stress(s)
    j2 = 0.5 * ((s[0] - p) ** 2 + (s[1] - p) ** 2 + (s[2] - p) ** 2) + s[3] ** 2 + s[4] ** 2 + s[5] ** 2
    return np.sqrt(3.0 * j2)


def mixed_step(dll, stress, strain, dstrain, groups, targets, props, tol=1e-10):
    """
    Strain increment in which every group of components (sharing one strain increment) is adjusted
    by a Newton iteration such that the stress of its first component equals the target.
    """
    d = np.array(dstrain, dtype=float)
    u = np.array([d[g[0]] for g in groups])
    first = [g[0] for g in groups]

    def apply(values):
        out = d.copy()
        for g, v in zip(groups, values):
            out[list(g)] = v
        return out

    for _ in range(50):
        s, ddsdde, _ = step(dll, stress, apply(u), props)
        r = s[first] - targets
        if np.max(np.abs(r)) < tol * (1.0 + np.max(np.abs(targets))):
            break
        # consistent tangent of the stress update, summed over the components of each group
        jac = np.array([[ddsdde[i, list(g)].sum() for g in groups] for i in first])
        u -= np.linalg.solve(jac, r)

    d = apply(u)
    s, _, statev = step(dll, stress, d, props)
    return s, statev, strain + d


def drained_triaxial(dll, props, cell, axial_strain, n_steps, extension=False):
    """
    Drained triaxial test from an isotropic state at constant cell pressure: compression along x
    (sigma_y = sigma_z = cell), or extension along z (sigma_x = sigma_y = cell).
    Returns the strains, stresses and states of all steps.
    """
    stress, strain = np.array([cell, cell, cell, 0.0, 0.0, 0.0]), np.zeros(6)
    if extension:
        dstrain, groups = np.array([0.0, 0.0, -axial_strain / n_steps, 0.0, 0.0, 0.0]), [[0, 1]]
    else:
        dstrain, groups = np.array([axial_strain / n_steps, 0.0, 0.0, 0.0, 0.0, 0.0]), [[1, 2]]
    strains, stresses, states = [strain], [stress], [0.0]
    for _ in range(n_steps):
        stress, statev, strain = mixed_step(dll, stress, strain, dstrain, groups, np.array([cell]), props)
        strains.append(strain)
        stresses.append(stress)
        states.append(statev[0])
    return np.array(strains), np.array(stresses), np.array(states)


def test_returned_stress_satisfies_matsuoka_nakai_criterion_at_all_lode_angles(dll):
    """
    The yield function of Lagioia & Panteghini (2016) is the Matsuoka-Nakai criterion
    I1 I2 / I3 = 9 + 8 tan^2(phi) of the shifted stresses: large deviatoric strain increments in
    all directions of the deviatoric plane are returned onto that surface, also with cohesion.
    """
    stress0 = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    for c, phi_deg in [(0.0, 20.0), (10.0, 40.0), (5.0, 30.0)]:
        props = props_of(c=c, phi=phi_deg, psi=phi_deg / 3.0)
        for lode in np.linspace(-np.pi / 6.0, np.pi / 6.0, 13):
            # deviatoric principal strain direction: lode = -pi/6 is triaxial compression along x,
            # +pi/6 triaxial extension along z
            direction = np.array([np.sin(lode + 2 * np.pi / 3), np.sin(lode), np.sin(lode - 2 * np.pi / 3)])
            s, _, statev = step(dll, stress0, 0.02 * np.r_[direction, 0.0, 0.0, 0.0], props)
            assert statev[0] == 1.0
            assert matsuoka_nakai_criterion(s, c, phi_deg) == pytest.approx(0.0, abs=1e-10)


@pytest.mark.parametrize("extension", [False, True])
@pytest.mark.parametrize("c", [0.0, 10.0])
def test_triaxial_strength_equals_mohr_coulomb(dll, extension, c):
    """
    The Matsuoka-Nakai surface passes through the corners of the Mohr-Coulomb pyramid: in drained
    triaxial compression and extension the strength is (sigma_1 + a) / (sigma_3 + a) =
    (1 + sin phi) / (1 - sin phi).
    """
    props = props_of(c=c)
    _, stresses, states = drained_triaxial(dll, props, 100.0, 0.04, 40, extension)

    s1, s2, s3 = principal_stresses(stresses[-1])
    a = attraction(c, PHI_DEG)
    sin_phi = np.sin(np.radians(PHI_DEG))
    assert states[-1] == 1.0
    assert (s1 + a) / (s3 + a) == pytest.approx((1.0 + sin_phi) / (1.0 - sin_phi), rel=1e-9)
    if extension:
        assert s1 == pytest.approx(s2, rel=1e-9)
        assert s1 == pytest.approx(100.0, rel=1e-9)
    else:
        assert s2 == pytest.approx(s3, rel=1e-9)
        assert s3 == pytest.approx(100.0, rel=1e-9)


@pytest.mark.parametrize("extension", [False, True])
@pytest.mark.parametrize("psi_deg", [-10.0, 0.0, 10.0, 25.0])
def test_dilatancy_at_triaxial_failure(dll, extension, psi_deg):
    """
    At failure in drained triaxial tests the stress is stationary, so the strain increments are
    plastic. At the compression and extension meridians the Matsuoka-Nakai plastic potential gives
    the dilatancy of the Mohr-Coulomb potential at its corners (Koiter's rule):
    -d eps_v / d(eps_1 - eps_3) = 4 sin(psi) / (3 -/+ sin(psi)), also for contractant flow
    (psi < 0).
    """
    props = props_of(psi=psi_deg)
    strains, stresses, _ = drained_triaxial(dll, props, 100.0, 0.06, 60, extension)

    ev = strains[:, :3].sum(axis=1)
    gamma = strains[:, 0] - (strains[:, 2] if extension else strains[:, 1])
    ratio = -(ev[-1] - ev[-11]) / (gamma[-1] - gamma[-11])
    sin_psi = np.sin(np.radians(psi_deg))
    expected = 4.0 * sin_psi / (3.0 + sin_psi if extension else 3.0 - sin_psi)

    assert np.allclose(stresses[-1], stresses[-11], rtol=1e-9, atol=1e-9)
    assert ratio == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("psi_deg", [PHI_DEG, 0.0])
def test_plane_strain_intermediate_stress_at_failure(dll, psi_deg):
    """
    Plane strain (eps_y = 0, sigma_z = cell): at failure the plastic strain in y vanishes. With
    associated flow (psi = phi) this gives sigma_2 + a = sqrt((sigma_1 + a)(sigma_3 + a)), the
    classical result of Matsuoka & Nakai; with psi = 0 the plastic potential is the von Mises
    cylinder and sigma_2 = (sigma_1 + sigma_3) / 2.
    """
    props = props_of(psi=psi_deg)
    cell = 100.0
    stress, strain = np.array([cell, cell, cell, 0.0, 0.0, 0.0]), np.zeros(6)
    # the stress approaches the stationary state at failure exponentially in the strain
    for _ in range(120):
        stress, statev, strain = mixed_step(dll, stress, strain, [2.5e-3, 0, 0, 0, 0, 0], [[2]],
                                            np.array([cell]), props)

    s1, s2, s3 = stress[0], stress[1], stress[2]
    a = attraction(C, PHI_DEG)
    assert statev[0] == 1.0
    assert s1 > s2 > s3
    assert matsuoka_nakai_criterion(stress, C, PHI_DEG) == pytest.approx(0.0, abs=1e-10)
    if psi_deg == PHI_DEG:
        assert s2 + a == pytest.approx(np.sqrt((s1 + a) * (s3 + a)), rel=1e-8)
    else:
        assert s2 == pytest.approx(0.5 * (s1 + s3), rel=1e-8)


@pytest.mark.parametrize("psi_deg", [0.0, 10.0])
def test_undrained_triaxial_compression(dll, psi_deg):
    """
    Undrained (isochoric) triaxial compression: without dilatancy the mean effective stress stays
    at its initial value and the path ends at q = 6 sin(phi) / (3 - sin(phi)) (p + a); with
    dilatancy p increases once the surface is reached.
    """
    props = props_of(psi=psi_deg)
    stress = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    d = 5e-4
    p, q = [mean_stress(stress)], [deviator_q(stress)]
    for _ in range(60):
        stress, _, statev = step(dll, stress, [d, -d / 2, -d / 2, 0, 0, 0], props)
        p.append(mean_stress(stress))
        q.append(deviator_q(stress))

    sin_phi = np.sin(np.radians(PHI_DEG))
    m_q = 6.0 * sin_phi / (3.0 - sin_phi)
    assert statev[0] == 1.0
    assert q[-1] == pytest.approx(m_q * (p[-1] + attraction(C, PHI_DEG)), rel=1e-9)
    if psi_deg == 0.0:
        assert np.allclose(p, 100.0, rtol=1e-12)
    else:
        assert np.all(np.diff(p) >= -1e-9)
        assert p[-1] > 150.0


@pytest.mark.parametrize("n_steps", [100, 20])
def test_undrained_contractant_flow_reaches_the_apex(dll, n_steps):
    """
    Undrained triaxial compression with contractant flow (psi < 0, c = 0): once the failure line is
    reached, the mean effective stress decreases along it (static liquefaction) until the stress
    reaches the apex, sigma = 0, where it stays under further shearing.
    """
    props = props_of(c=0.0, psi=-10.0)
    stress = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    d = 0.05 / n_steps
    p, states = [mean_stress(stress)], []
    for _ in range(n_steps):
        stress, _, statev = step(dll, stress, [d, -d / 2, -d / 2, 0, 0, 0], props)
        p.append(mean_stress(stress))
        states.append(statev[0])
        if statev[0] == 1.0:
            assert matsuoka_nakai_criterion(stress, 0.0, PHI_DEG) == pytest.approx(0.0, abs=1e-10)

    assert np.all(np.diff(p) <= 1e-9)
    assert 1.0 in states
    assert states[-1] == 2.0
    assert np.allclose(stress, 0.0, atol=1e-10)


def test_zero_friction_angle_gives_von_mises(dll):
    """
    For phi = 0 the Matsuoka-Nakai surface is the von Mises cylinder through the corners of the
    Tresca hexagon: q = 2 c in triaxial compression, triaxial extension and simple shear.
    """
    c = 20.0
    props = props_of(c=c, phi=0.0, psi=0.0)
    stress0 = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    for dstrain in ([0.01, -0.005, -0.005, 0, 0, 0], [0.005, 0.005, -0.01, 0, 0, 0], [0, 0, 0, 0.02, 0, 0]):
        s, _, statev = step(dll, stress0, dstrain, props)
        assert statev[0] == 1.0
        assert deviator_q(s) == pytest.approx(2.0 * c, rel=1e-10)
        assert mean_stress(s) == pytest.approx(100.0, rel=1e-12)


@pytest.mark.parametrize("psi_deg", [0.0, 10.0, 30.0])
@pytest.mark.parametrize("c", [0.0, 5.0])
def test_tangent_matches_numerical_derivative(dll, psi_deg, c):
    """
    DDSDDE is the consistent tangent of the implicit return mapping: it equals the numerical
    derivative of the stress update with respect to the strain increment, for a large plastic
    increment in a general direction.
    """
    props = props_of(c=c, psi=psi_deg)
    stress = np.array([120.0, 90.0, 100.0, 10.0, -5.0, 8.0])
    dstrain = np.array([1.2e-2, -6e-3, 3e-3, 9e-3, -3e-3, 6e-3])
    _, ddsdde, statev = step(dll, stress, dstrain, props)
    assert statev[0] == 1.0

    h = 1e-7
    numerical = np.zeros((6, 6))
    for j in range(6):
        e = np.zeros(6)
        e[j] = h
        numerical[:, j] = (step(dll, stress, dstrain + e, props)[0] - step(dll, stress, dstrain - e, props)[0]) / (2 * h)

    assert np.linalg.norm(ddsdde - numerical) / np.linalg.norm(numerical) < 1e-6
    asymmetry = np.linalg.norm(ddsdde - ddsdde.T) / np.linalg.norm(ddsdde)
    if psi_deg == PHI_DEG:
        assert asymmetry < 1e-12  # associated flow
    else:
        assert asymmetry > 1e-3


def test_principal_frame_invariance(dll):
    """
    The model is isotropic: the same loading path applied in a rotated frame gives the rotated
    stress. (This failed before the third invariant J3 = det(s) was corrected for stress states
    with shear components.)
    """
    a, b, g = 0.4, -0.7, 1.1
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    Rx = np.array([[1, 0, 0], [0, np.cos(g), -np.sin(g)], [0, np.sin(g), np.cos(g)]])
    rotations = [Rz @ Ry @ Rx, np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])]
    props = props_of()
    dstrain_local = np.array([6e-4, -2.5e-4, -1e-4, 0.0, 0.0, 0.0])

    def run(R, n=20):
        stress = rotate_stress(np.array([100.0, 80.0, 60.0, 0.0, 0.0, 0.0]), R)
        for _ in range(n):
            stress, _, statev = step(dll, stress, rotate_strain(dstrain_local, R), props)
        return rotate_stress(stress, R.T), statev

    stress_ref, statev_ref = run(np.eye(3))
    assert statev_ref[0] == 1.0
    for R in rotations:
        stress_rot, statev_rot = run(R)
        assert np.allclose(stress_rot, stress_ref, rtol=1e-9, atol=1e-8)
        assert statev_rot[0] == statev_ref[0]


@pytest.mark.parametrize("c, psi_deg", [(10.0, 0.0), (10.0, 10.0), (0.0, 10.0)])
def test_tension_beyond_apex_returns_to_apex(dll, c, psi_deg):
    """
    Trial stresses in tension beyond the apex of the cone are returned to the apex,
    sigma = -a (isotropic tension c cot(phi)), with a fraction 1e-2 of the elastic stiffness as
    tangent. Compression from the apex is elastic again.
    """
    props = props_of(c=c, psi=psi_deg)
    a = attraction(c, PHI_DEG)
    stress0 = np.array([50.0, 50.0, 50.0, 0.0, 0.0, 0.0])
    for dstrain in ([-1e-2, -1e-2, -1e-2, 0, 0, 0], [-1e-2, -5e-3, -8e-3, 2e-3, 0, -1e-3]):
        s, ddsdde, statev = step(dll, stress0, dstrain, props)
        assert statev[0] == 2.0
        assert np.allclose(s, [-a, -a, -a, 0.0, 0.0, 0.0], rtol=0.0, atol=1e-10 * (a + 1.0))
        assert np.allclose(ddsdde, 1e-2 * elastic_stiffness(), rtol=1e-12)

    reload = np.array([1e-4, 1e-4, 1e-4, 0.0, 0.0, 0.0])
    s_new, ddsdde, statev = step(dll, s, reload, props)
    assert statev[0] == 0.0
    assert np.allclose(s_new - s, elastic_stiffness() @ reload, rtol=1e-10)
    assert np.allclose(ddsdde, elastic_stiffness(), rtol=1e-12)


@pytest.mark.parametrize("psi_deg", [PSI_DEG, -10.0])
def test_large_increments_return_to_the_surface(dll, psi_deg):
    """
    Large strain increments in arbitrary directions (up to the order of 10 % strain) end on the
    yield surface, at the apex, or elastically inside the surface, without a failure of the local
    iteration (which would leave the stress unchanged), also with contractant flow.
    """
    rng = np.random.default_rng(2024)
    props = props_of(psi=psi_deg)
    a = attraction(C, PHI_DEG)
    n_plastic = 0
    for _ in range(200):
        stress0 = np.r_[np.full(3, rng.uniform(10.0, 300.0)), 0.0, 0.0, 0.0] + rng.normal(size=6) * 20.0
        dstrain = rng.normal(size=6) * 10.0 ** rng.uniform(-4.0, -1.0)
        s, _, statev = step(dll, stress0, dstrain, props)
        assert np.all(np.isfinite(s))
        if statev[0] == 1.0:
            n_plastic += 1
            assert matsuoka_nakai_criterion(s, C, PHI_DEG) == pytest.approx(0.0, abs=1e-9)
        elif statev[0] == 2.0:
            assert np.allclose(s, [-a, -a, -a, 0.0, 0.0, 0.0], atol=1e-9)
        else:
            assert np.allclose(s, stress0 + elastic_stiffness() @ dstrain, rtol=1e-12, atol=1e-9)
            assert matsuoka_nakai_criterion(s, C, PHI_DEG) <= 1e-9
    assert n_plastic >= 50


def test_energy_terms(dll):
    """
    SSE is updated with the change of the elastic strain energy 1/2 sigma : Ce^-1 sigma and SPD
    with the plastic work sigma : d eps_p of the increment.
    """
    ffi = cffi.FFI()
    ffi.cdef("""
    void umat(double* STRESS, double* STATEV, double* DDSDDE, double* SSE, double* SPD, double* SCD,
              double* RPL, double* DDSDDT, double* DRPLDE, double* DRPLDT, double* STRAN, double* DSTRAN,
              double* TIME, double* DTIME, double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
              char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV, double* PROPS, int* NPROPS,
              double* COORDS, double* DROT, double* PNEWDT, double* CELENT, double* DFGRD0,
              double* DFGRD1, int* NOEL, int* NPT, int* LAYER, int* KSPT, int* KSTEP, int* KINC);
    """)
    lib = ffi.dlopen(dll)
    props = props_of()
    compliance = np.linalg.inv(elastic_stiffness())
    stress0 = -np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])  # tension positive
    for dstrain in (np.array([-1e-4, 0.0, 0.0, 0.0, 0.0, 0.0]),  # elastic
                    np.array([-1e-2, 4e-3, 4e-3, 2e-3, 0.0, 0.0])):  # plastic
        c_stress = ffi.new("double[]", list(stress0))
        sse, spd = ffi.new("double*", 1.0), ffi.new("double*", 2.0)
        ints = [ffi.new("int*", v) for v in (3, 3, 6, 1, 5, 1, 1, 1, 1)]
        lib.umat(c_stress, ffi.new("double[]", [0.0]), ffi.new("double[]", 36), sse, spd, ffi.new("double*", 0.0),
                 ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL, ffi.new("double[]", 6), ffi.new("double[]", list(dstrain)),
                 ffi.new("double[]", [0.0, 0.0]), ffi.new("double*", 1.0), ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL,
                 ffi.new("char[]", b"MN".ljust(80)), ints[0], ints[1], ints[2], ints[3],
                 ffi.new("double[]", props), ints[4], ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL,
                 ints[5], ints[6], ffi.NULL, ffi.NULL, ints[7], ints[8])
        s = np.array(list(c_stress))
        d_eps_p = dstrain - compliance @ (s - stress0)
        assert sse[0] - 1.0 == pytest.approx(0.5 * s @ compliance @ s - 0.5 * stress0 @ compliance @ stress0,
                                             rel=1e-12)
        assert spd[0] - 2.0 == pytest.approx(s @ d_eps_p, rel=1e-10, abs=1e-14)
    assert spd[0] - 2.0 > 0.0


def test_unloading_from_failure_is_elastic(dll):
    """After failure, reversing the strain increment unloads elastically with the elastic tangent."""
    props = props_of()
    stress = np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0])
    dstrain = np.array([2e-3, -1e-3, -1e-3, 0.0, 0.0, 0.0])
    for _ in range(10):
        stress, _, statev = step(dll, stress, dstrain, props)
    assert statev[0] == 1.0

    s_new, ddsdde, statev = step(dll, stress, -0.1 * dstrain, props)
    assert statev[0] == 0.0
    assert np.allclose(s_new - stress, elastic_stiffness() @ (-0.1 * dstrain), rtol=1e-10)
    assert np.allclose(ddsdde, elastic_stiffness(), rtol=1e-12)


def test_stress_update_is_continuous_in_strain_increment(dll):
    """
    The integrated stress is a continuous function of the strain increment, also across the
    transition from elastic to plastic increments, so that a global Newton iteration converges.
    """
    props = props_of()
    stress = np.array([180.0, 90.0, 100.0, 10.0, 0.0, -5.0])
    dstrain = np.array([2e-3, -5e-4, -1e-3, 1e-3, 0.0, 5e-4])
    factors = np.linspace(0.0, 4.0, 801)  # plastic from about 2.1
    response = np.array([step(dll, stress, t * dstrain, props)[0] for t in factors])
    states = np.array([step(dll, stress, t * dstrain, props)[2][0] for t in factors[::100]])
    changes = np.linalg.norm(np.diff(response, axis=0), axis=1)

    assert states[0] == 0.0 and states[-1] == 1.0
    # a jump shows up as one change much larger than both of its neighbours (the slope itself
    # changes at the elastic-plastic transition)
    assert np.max(changes[1:-1] / np.maximum(changes[:-2], changes[2:])) < 1.05