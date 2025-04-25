import os
import sys
import numpy as np

from tests.incr_driver import IncrDriver



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
        model_loc = rf"{project_dir}\build_C\lib\matsuoka_nakai.dll"
    elif sys.platform == 'linux':
        model_loc = f"{project_dir}/build_C/lib/matsuoka_nakai.so"
    else:
        raise Exception("Unsupported operating system")

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
        model_loc = f"{project_dir}/build_C/lib/matsuoka_nakai.dll"
    elif sys.platform == 'linux':
        model_loc = f"{project_dir}/build_C/lib/matsuoka_nakai.so"
    else:
        raise Exception("Unsupported operating system")

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