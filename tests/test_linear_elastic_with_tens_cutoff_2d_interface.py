import os
import sys
import numpy as np

from tests.incr_driver import IncrDriver


def test_strain_controlled_compression_oedometer():
    """
    Test the strain controlled compression oedometer  with a tensile cutoff in the vertical direction. The sample
    is compressed, thus the tensile cutoff is not active, the result should be the same as the linear elastic case.

    :return:
    """

    G = 2000
    nu = 0.33
    E = 2 * G * (1 + nu)

    M = E * (1-nu)/((1+nu)*(1-2*nu)) # P wave modulus

    params = {'E': E, 'poison_ratio': nu, 'tensile_cutoff': 1}
    state_vars={'is_plastic': 0}

    # Define the original stress vector
    orig_stress_vector = np.array([-1, 0])

    project_dir = os.getcwd()
    # check operating system
    if sys.platform == 'win32':
        extension = 'dll'
    elif sys.platform == 'linux':
        extension = 'so'
    else:
        raise Exception("Unsupported operating system")
    model_loc = os.path.join(project_dir,'build_C','lib','linear_elastic_with_vertical_tens_cutoff_2d_interface.'+extension)

    const_model_info = {"language": "c",
                        "file_name":model_loc,
                        "properties": list(params.values()),
                        "state_vars": list(state_vars.values())}

    # strain increment per time step in the form of [eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz]

    vertical_strain_increment = -1e-4

    strain_increment = np.array([vertical_strain_increment, 0])
    stress_increment = np.zeros_like(strain_increment)

    # solve the problem
    incr_driver = IncrDriver(orig_stress_vector,
                             strain_increment,
                             stress_increment,
                             const_model_info,
                             3,
                             5, ndim=2, n_direct_stress_components=1, n_shear_components=1,
                             vertical_axis_index=0)

    incr_driver.solve_oedometer_strain_controlled()

    expected_strains = np.array([[vertical_strain_increment, 0],
                                    [2*vertical_strain_increment, 0],
                                    [3*vertical_strain_increment, 0]])

    expected_stresses = np.array([orig_stress_vector,orig_stress_vector,orig_stress_vector], dtype=float)
    expected_stresses[0,0] = expected_stresses[0,0]+ vertical_strain_increment * M
    expected_stresses[1,0] = expected_stresses[1,0] + vertical_strain_increment * M * 2
    expected_stresses[2,0] = expected_stresses[2,0] + vertical_strain_increment * M * 3

    # Check the results
    np.testing.assert_allclose(incr_driver.strains, expected_strains, rtol=1e-6)
    np.testing.assert_allclose(incr_driver.stresses, expected_stresses, rtol=1e-6)


def test_stress_controlled_extension_below_cutoff():
    """
    Test a stress controlled test with a tensile cutoff in the vertical direction. The sample is extended, but less than
    the tensile cutoff, thus the tensile cutoff is not active, the result should be the same as the linear elastic case.

    :return:
    """

    E = 30000000
    nu = 0.0

    M = E * (1-nu)/((1+nu)*(1-2*nu)) # P wave modulus
    G = E / (2 * (1 + nu))

    params = {'E': E, 'poison_ratio': nu, 'tensile_cutoff': 10000000}
    state_vars={'is_plastic': 0}

    # Define the original stress vector
    orig_stress_vector = np.array([0, 0])

    project_dir = os.getcwd()
    # check operating system
    if sys.platform == 'win32':
        extension = 'dll'
    elif sys.platform == 'linux':
        extension = 'so'
    else:
        raise Exception("Unsupported operating system")
    model_loc = os.path.join(project_dir,'build_C','lib','linear_elastic_with_vertical_tens_cutoff_2d_interface.'+extension)

    const_model_info = {"language": "c",
                        "file_name":model_loc,
                        "properties": list(params.values()),
                        "state_vars": list(state_vars.values())}

    # stress increment per time step in the form of [sigma_zz, sigma_xz]
    strain_increment = np.array([0, 0])
    stress_increment = np.array([333,667])

    # solve the problem
    incr_driver = IncrDriver(orig_stress_vector,
                             strain_increment,
                             stress_increment,
                             const_model_info,
                             3,
                             5, ndim=2, n_direct_stress_components=1, n_shear_components=1,
                             vertical_axis_index=0)

    incr_driver.solve([1,1])

    expected_strains = np.array([-stress_increment / np.array([M,G]),
                                 -2 * stress_increment / np.array([M,G]),
                                 -3* stress_increment / np.array([M,G])])

    expected_stresses = np.array([-stress_increment,-stress_increment*2,-stress_increment*3], dtype=float)


    # Check the results
    np.testing.assert_allclose(incr_driver.strains, expected_strains, rtol=1e-6)
    np.testing.assert_allclose(incr_driver.stresses, expected_stresses, rtol=1e-6)



def test_strain_controlled_tension_oedometer():
    """
    Test the strain controlled oedometer test with a tensile cutoff in the vertical direction. Enough tension is
    applied such that the tensile cutoff is active.

    :return:
    """

    G = 2000
    nu = 0.33
    E = 2 * G * (1 + nu)

    M = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))  # P wave modulus

    tensile_cutoff = 0.6
    params = {'E': E, 'poison_ratio': nu, 'tensile_cutoff': tensile_cutoff}
    state_vars={'is_plastic': 0}

    # Define the original stress vector
    orig_stress_vector = np.array([-1.0, 0.0])

    project_dir = os.getcwd()
    # check operating system
    if sys.platform == 'win32':
        extension = 'dll'
    elif sys.platform == 'linux':
        extension = 'so'
    else:
        raise Exception("Unsupported operating system")
    model_loc = os.path.join(project_dir,'build_C','lib','linear_elastic_with_vertical_tens_cutoff_2d_interface.'+extension)

    const_model_info = {"language": "c",
                        "file_name":model_loc,
                        "properties": list(params.values()),
                        "state_vars": list(state_vars.values())}

    # strain increment per time step in the form of [eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz]
    vertical_strain_increment = 1e-4
    strain_increment = np.array([vertical_strain_increment, 0.0])
    stress_increment = np.zeros_like(strain_increment)

    # solve the problem
    incr_driver = IncrDriver(orig_stress_vector,
                             strain_increment,
                             stress_increment,
                             const_model_info,
                             3,
                             5, ndim=2, n_direct_stress_components=1, n_shear_components=1,
                             vertical_axis_index=0)

    incr_driver.solve_oedometer_strain_controlled()

    expected_strains = np.array([[vertical_strain_increment, 0.0],
                                    [2*vertical_strain_increment, 0.0],
                                    [3*vertical_strain_increment, 0.0]])

    expected_stresses = np.array([orig_stress_vector,orig_stress_vector,orig_stress_vector], dtype=float)
    expected_stresses[0,0] = expected_stresses[0,0]+ vertical_strain_increment * M
    expected_stresses[1,0] = expected_stresses[1,0] + vertical_strain_increment * M * 2
    expected_stresses[2,0] = tensile_cutoff

    # Check the results
    np.testing.assert_allclose(incr_driver.strains, expected_strains, rtol=1e-6)
    np.testing.assert_allclose(incr_driver.stresses, expected_stresses, rtol=1e-6)