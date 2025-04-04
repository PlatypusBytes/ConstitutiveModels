import os

import numpy as np


from tests.incr_driver import IncrDriver


def test_pressure():

    E = 2000
    nu = 0.33
    params = {'E': E, 'poison_ratio': nu, 'tensile_cutoff': 0}

    # Define the control type for the stress vector, where 1 indicates a controlled variable and 0 indicates an uncontrolled variable
    control_type = [1, 0, 1, 1, 1, 1]
    stress_vector = np.array([-1, -1, -1, 0, 0, 0])
    orig_stress_vector = np.copy(stress_vector)

    project_dir = os.getcwd()
    const_model_info = {"language": "c",
                        "file_name": f"{project_dir}\c_models\linear_elastic_with_vertical_tens_cutoff\linear_elastic_with_vertical_tens_cutoff.dll",
                        "properties": list(params.values()),}

    delta_strain = np.array([0, 0.0001, 0, 0, 0, 0])
    incr_driver = IncrDriver(orig_stress_vector,delta_strain,control_type,const_model_info,10,100 )
    incr_driver.solve()

    print(incr_driver.strains)

    a=1+1
if __name__ == "__main__":
    test_pressure()