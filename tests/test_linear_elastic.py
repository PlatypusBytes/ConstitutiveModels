import os
import numpy as np

from tests.incr_driver import IncrDriver


def test_strain_controlled_compression_triaxial():
    """
    Test the strain controlled triaxial test with a tensile cutoff in the vertical direction. Where the sample is
    compressed, thus the tensile cutoff is not active, the result should be the same as the linear elastic case.

    :return:
    """

    G = 2000
    nu = 0.33
    E = 2 * G * (1 + nu)

    params = {'E': E, 'poison_ratio': nu, 'tensile_cutoff': 1}
    state_vars={'is_plastic': 0}

    # Define the original stress vector
    orig_stress_vector = np.array([-1., -1., -1., 0, 0, 0])

    project_dir = os.getcwd()
    const_model_info = {"language": "c",
                        "file_name": f"{project_dir}\c_models\linear_elastic_with_vertical_tens_cutoff\linear_elastic_with_vertical_tens_cutoff.dll",
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
                             3,
                             5)

    incr_driver.solve_triaxial_strain_controlled()

    expected_strains = np.array([[-vertical_strain_increment*nu, vertical_strain_increment, -vertical_strain_increment*nu, 0, 0, 0],
                                    [2*-vertical_strain_increment*nu, 2*vertical_strain_increment, 2*-vertical_strain_increment*nu, 0, 0, 0],
                                    [3*-vertical_strain_increment*nu, 3*vertical_strain_increment, 3*-vertical_strain_increment*nu, 0, 0, 0]])

    expected_stresses = np.array([orig_stress_vector,orig_stress_vector,orig_stress_vector], dtype=float)
    expected_stresses[0,1] = expected_stresses[0,1]+ vertical_strain_increment * E
    expected_stresses[1,1] = expected_stresses[1,1] + vertical_strain_increment * E * 2
    expected_stresses[2,1] = expected_stresses[2,1] + vertical_strain_increment * E * 3

    # Check the results
    np.testing.assert_allclose(incr_driver.strains, expected_strains, rtol=1e-6)
    np.testing.assert_allclose(incr_driver.stresses, expected_stresses, rtol=1e-6)


def test_strain_controlled_tension_triaxial():
    """
    Test the strain controlled triaxial test with a tensile cutoff in the vertical direction. Where enough tension is
    applied such that the tensile cutoff is active.

    :return:
    """

    G = 2000
    nu = 0.33
    E = 2 * G * (1 + nu)

    tensile_cutoff = 0.5
    params = {'G': G, 'poison_ratio': nu}
    state_vars={}

    # Define the original stress vector
    orig_stress_vector = np.array([-1., -1., -1., 0, 0, 0])

    project_dir = os.getcwd()
    const_model_info = {"language": "fortran",
                        "file_name": rf"{project_dir}\fortran_models\linear_elastic\linear_elastic.dll",
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
                             5)

    incr_driver.solve_triaxial_strain_controlled()

    expected_strains = np.array([[-vertical_strain_increment*nu, vertical_strain_increment, -vertical_strain_increment*nu, 0, 0, 0],
                                    [2*-vertical_strain_increment*nu, 2*vertical_strain_increment, 2*-vertical_strain_increment*nu, 0, 0, 0],
                                    [3*-vertical_strain_increment*nu, 3*vertical_strain_increment, 3*-vertical_strain_increment*nu, 0, 0, 0]])

    expected_stresses = np.array([orig_stress_vector,orig_stress_vector,orig_stress_vector], dtype=float)
    expected_stresses[0,1] = expected_stresses[0,1]+ vertical_strain_increment * E
    expected_stresses[1,1] = expected_stresses[1,1] + vertical_strain_increment * E * 2
    expected_stresses[2,1] = expected_stresses[2,1] + vertical_strain_increment * E * 3

    # Check the results
    np.testing.assert_allclose(incr_driver.strains, expected_strains, rtol=1e-6)
    np.testing.assert_allclose(incr_driver.stresses, expected_stresses, rtol=1e-6)