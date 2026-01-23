import os
import sys
import pytest
import itertools
import numpy as np

from tests.incr_driver import IncrDriver


PLOTS = False

# Define parameter combinations
amplitudes = [0.002, -0.02]
frequencies = [1, 5]
E_values = [10e6, 200e6]
eta_values = [1e5, 10e5]

# Create all combinations
test_params = list(itertools.product(amplitudes, frequencies, E_values, eta_values))

@pytest.mark.parametrize("amplitude,frequency,E,eta", test_params)
def test_strain_controlled_cyclic_triaxial(amplitude, frequency, E, eta):
    """
    Test the strain controlled triaxial test with cyclic loading using the Kelvin-Voigt viscoelastic model.

    :return:
    """

    # Material parameters
    nu = 0.3  # Poisson's Ratio

    # Strain parameters
    omega = 2.0 * np.pi * frequency
    total_time = 2.0
    delta_t = 0.005

    time = np.arange(0, total_time, delta_t)
    num_steps = int(total_time / delta_t)
    strain = np.zeros((num_steps, 6))
    strain[:, 1] = amplitude * np.sin(omega * np.linspace(0, total_time, num_steps))
    strain_increment = np.diff(strain, axis=0, prepend=np.zeros((1, 6)))

    params = {'E': E, 'poisson_ratio': nu, 'eta': eta}
    state_vars={'elastic_stress': np.zeros(6)}

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
    model_loc = os.path.join(project_dir,'build_C','lib','viscoelastic_kelvin_voigt.'+extension)

    const_model_info = {"language": "c",
                        "file_name":model_loc,
                        "properties": list(params.values()),
                        "state_vars": list(state_vars.values())[0]}

    # stress increment
    stress_increment = np.zeros_like(strain_increment)

    # solve the problem
    incr_driver = IncrDriver(orig_stress_vector,
                             strain_increment,
                             stress_increment,
                             const_model_info,
                             len(time),
                             5,
                             delta_t,
                             )

    incr_driver.solve_triaxial_strain_controlled()


    strain_, stress_ = kelvin_voigt_triax(orig_stress_vector[:3], strain[:, :3], E, nu, eta, time)

    # Check the results
    np.testing.assert_allclose(incr_driver.strains[:, 1], strain_, atol=1e-16, rtol=1e-2)
    np.testing.assert_allclose(incr_driver.stresses[:, 1], stress_, atol=1e-16, rtol=1e-2)

    if PLOTS:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(8, 6))
        ax[0].plot(incr_driver.strains[:, 1], incr_driver.stresses[:, 1], label='UMAT Result')
        ax[0].plot(strain_, stress_, label='Analytical Solution')
        ax[0].set_xlabel('Axial Strain [-]')
        ax[0].set_ylabel('Axial Stress [Pa]')
        ax[0].grid()
        ax[0].legend()
        ax[1].plot(time, incr_driver.strains[:, 1], label='UMAT Strain')
        ax[1].plot(time, strain_, label='Analytical Strain')
        ax[1].set_xlabel('Time [s]')
        ax[1].set_ylabel('Axial Strain [-]')
        ax[1].grid()
        ax[1].legend()
        plt.show()


def kelvin_voigt_triax(ini_stress, strain, E, nu, eta, time):
    """
    Auxiliary function to compute the analytical solution for triaxial test with Kelvin-Voigt model.

    Args:
        ini_stress (array): Initial stress state [sigma_xx, sigma_yy, sigma_zz]
        strain (array): Strain over time steps (n_steps x 3)
        E (float): Young's Modulus
        nu (float): Poisson's Ratio
        eta (float): Viscosity coefficient
        time (array): Time array
    Returns:
        eps_xx (array): Axial strain over time
        sigma_xx (array): Axial stress over time
    """

    # D matrix for isotropic linear elasticity in 3D
    fct = E / ((1 + nu) * (1 - 2 * nu))
    D = fct * np.array([[1 - nu, nu,     nu    ],
                        [nu,     1 - nu, nu    ],
                        [nu,     nu,     1 - nu]
                        ])

    n_steps = len(time)
    t_step = time[1] - time[0]

    eps_xx = np.zeros(n_steps)
    eps_yy = np.zeros(n_steps)
    eps_zz = np.zeros(n_steps)

    sigma_xx = np.zeros(n_steps)
    sigma_yy = np.zeros(n_steps)
    sigma_zz = np.zeros(n_steps)

    # Initial conditions
    eps = np.zeros(3)
    # sigma = np.array(ini_stress)  # initial stress

    for i, t in enumerate(time):
        # Prescribed axial cyclic strain
        eps_xx[i] = eps[0] = strain[i, 1]

        # Lateral strain unknowns: eps_yy, eps_zz
        A = D[1:,1:] + (eta/t_step) * np.eye(2)

        # initialize b
        if i == 0:
            b = np.array([ini_stress[1] - D[1,0]*eps[0],
                          ini_stress[2] - D[2,0]*eps[0]
                          ])
            sigma = D @ eps
        else:
            b = np.array([
                ini_stress[1] - D[1,0]*eps[0] + (eta/t_step)*eps_yy[i-1],
                ini_stress[2] - D[2,0]*eps[0] + (eta/t_step)*eps_zz[i-1]
            ])
            sigma = D @ eps + eta * np.array([(eps[0] - eps_xx[i-1])/t_step,
                                              (eps[1] - eps_yy[i-1])/t_step,
                                              (eps[2] - eps_zz[i-1])/t_step
                                              ])

        eps_lat = np.linalg.solve(A, b)
        eps[1] = eps_yy[i] = eps_lat[0]
        eps[2] = eps_zz[i] = eps_lat[1]

        # Axial stress including viscous contribution
        sigma = D @ eps + eta * np.array([
            (eps[0] - eps_xx[i-1])/t_step,
            (eps[1] - eps_yy[i-1])/t_step,
            (eps[2] - eps_zz[i-1])/t_step
        ])

        sigma_xx[i] = sigma[0]
        sigma_yy[i] = sigma[1]
        sigma_zz[i] = sigma[2]

    return eps_xx, sigma_xx
