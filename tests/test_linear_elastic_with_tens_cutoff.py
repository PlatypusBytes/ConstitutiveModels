import os

import numpy as np


from tests.incr_driver import IncrDriver


def test_strain_controlled_pressure_triaxial():

    E = 2000
    nu = 0.33
    params = {'E': E, 'poison_ratio': nu, 'tensile_cutoff': 0}
    state_vars={'is_plastic': 0}

    # Define the original stress vector
    orig_stress_vector = np.array([-1, -1, -1, 0, 0, 0])

    project_dir = os.getcwd()
    const_model_info = {"language": "c",
                        "file_name": f"{project_dir}\c_models\linear_elastic_with_vertical_tens_cutoff\linear_elastic_with_vertical_tens_cutoff.dll",
                        "properties": list(params.values()),
                        "state_vars": list(state_vars.values())}

    # strain increment per time step in the form of [eps_x, eps_y, eps_z, gamma_xy, gamma_xz, gamma_yz]
    strain_increment = np.array([0, -0.0001, 0, 0, 0, 0])
    stress_increment = np.zeros_like(strain_increment)

    # solve the problem
    incr_driver = IncrDriver(orig_stress_vector,
                             strain_increment,
                             stress_increment,
                             const_model_info,
                             10,
                             100)

    incr_driver.solve_triaxial_strain_controlled()

    print(incr_driver.strains)


    a=1+1
if __name__ == "__main__":
    test_strain_controlled_pressure_triaxial()