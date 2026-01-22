import os
import sys
from time import time
import numpy as np

from tests.incr_driver import IncrDriver


PLOTS = False


def test_strain_controlled_compression_triaxial():
    """
    Test the strain controlled triaxial test with a tensile cutoff in the vertical direction. Where enough tension is
    applied such that the tensile cutoff is active.

    :return:
    """


    # Material parameters
    E = 10e6  # Young's Modulus in Pa
    nu = 0.3  # Poisson's Ratio
    eta = 5e5  # Viscosity in Pa.s
    # Strain parameters
    amplitude = 0.0002
    frequency = 1.0
    omega = 2.0 * np.pi * frequency
    total_time = 2.0
    delta_t = 0.005

    time = np.arange(0, total_time, delta_t)
    num_steps = int(total_time / delta_t)
    strain = np.zeros((num_steps, 6))
    strain[:, 1] = amplitude * np.sin(omega * np.linspace(0, total_time, num_steps))
    strain_increment = np.diff(strain, axis=0, prepend=np.zeros((1, 6)))


    params = {'E': E, 'poisson_ratio': nu, 'eta': eta}
    state_vars={'strain': np.zeros(6)}

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
                             delta_t,
                             5)

    incr_driver.solve_triaxial_strain_controlled()


    strain, stress = kelvin_voigt_triax(orig_stress_vector[:3], strain, E, nu, eta, time)

    # import matplotlib.pyplot as plt
    # plt.plot(incr_driver.strains[:, 1], incr_driver.stresses[:, 1], label='UMAT Result')
    # plt.plot(strain, stress, label='Analytical Solution')
    # plt.grid()
    # plt.legend()
    # plt.show()

    # Check the results
    np.testing.assert_allclose(incr_driver.strains[:, 1], strain[:], atol=1e-16, rtol=1e-2)
    np.testing.assert_allclose(incr_driver.stresses[:, 1], stress[:], atol=1e-16, rtol=1e-2)

    if PLOTS:
        import matplotlib.pyplot as plt
        plt.plot(incr_driver.strains[:, 1], incr_driver.stresses[:, 1], label='UMAT Result')
        plt.plot(strain, stress, label='Analytical Solution')
        plt.xlabel('Axial Strain')
        plt.ylabel('Axial Stress (Pa)')
        plt.grid()
        plt.legend()
        plt.show()


    # plt.plot(time, incr_driver.strains[:, 0])
    # plt.plot(time, incr_driver.strains[:, 1])
    # plt.plot(time, incr_driver.strains[:, 2])
    # plt.plot(time, incr_driver.strains[:, 3])
    # plt.show()
    # print(1)
    # expected_strains = np.array([[-vertical_strain_increment*nu, vertical_strain_increment, -vertical_strain_increment*nu, 0, 0, 0],
    #                                 [2*-vertical_strain_increment*nu, 2*vertical_strain_increment, 2*-vertical_strain_increment*nu, 0, 0, 0],
    #                                 [2*-vertical_strain_increment*nu, 3*vertical_strain_increment, 2*-vertical_strain_increment*nu, 0, 0, 0]])

    # expected_stresses = np.array([orig_stress_vector,orig_stress_vector,orig_stress_vector], dtype=float)
    # expected_stresses[0,1] = expected_stresses[0,1]+ vertical_strain_increment * E
    # expected_stresses[1,1] = expected_stresses[1,1] + vertical_strain_increment * E * 2
    # expected_stresses[2,1] = tensile_cutoff

    # # Check the results
    # np.testing.assert_allclose(incr_driver.strains, expected_strains, rtol=1e-6)
    # np.testing.assert_allclose(incr_driver.stresses, expected_stresses, rtol=1e-6)


def kelvin_voigt_triax(ini_stress, strain_increment, E, nu, eta, time):

    D = (E / ((1 + nu) * (1 - 2 * nu))) * np.array([
        [1 - nu,  nu,       nu      ],
        [nu,      1 - nu,   nu      ],
        [nu,      nu,       1 - nu  ]
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
    sigma = np.array(ini_stress)  # initial stress

    for i, t in enumerate(time):
        # Prescribed axial cyclic strain
        eps_xx[i] = eps[0] = strain_increment[i, 1]

        # Compute lateral strains (eps_yy, eps_zz) from triaxial constraint
        # Solve: sigma_lat = sigma_conf = E*eps_lat + eta*d_eps_lat/dt + coupling
        # Here we do simple explicit step for viscous part

        # Define A*eps_lat = sigma_lat - E coupling
        # Lateral strain unknowns: eps_yy, eps_zz
        A = D[1:,1:] + (eta/t_step) * np.eye(2)
        # b = np.array([
        #     ini_stress[1] - D[1,0]*eps[0] + (eta/t_step)*eps_yy[i-1] if i>0 else 0,
        #     ini_stress[2] - D[2,0]*eps[0] + (eta/t_step)*eps_zz[i-1] if i>0 else 0
        # ])


        if i == 0:
            b = np.array([
                ini_stress[1] - D[1,0]*eps[0],
                ini_stress[2] - D[2,0]*eps[0]
            ])
            sigma = D @ eps
        else:
            b = np.array([
                ini_stress[1] - D[1,0]*eps[0] + (eta/t_step)*eps_yy[i-1],
                ini_stress[2] - D[2,0]*eps[0] + (eta/t_step)*eps_zz[i-1]
            ])
            sigma = D @ eps + eta * np.array([
                (eps[0] - eps_xx[i-1])/t_step,
                (eps[1] - eps_yy[i-1])/t_step,
                (eps[2] - eps_zz[i-1])/t_step
            ])

        eps_lat = np.linalg.solve(A, b)
        eps[1] = eps_yy[i] = eps_lat[0]
        eps[2] = eps_zz[i] = eps_lat[1]

        # Axial stress including viscous contribution
        sigma = D @ eps + eta * np.array([
            (eps[0] - eps_xx[i-1])/t_step if i>0 else 0,
            (eps[1] - eps_yy[i-1])/t_step if i>0 else 0,
            (eps[2] - eps_zz[i-1])/t_step if i>0 else 0
        ])

        sigma_xx[i] = sigma[0]
        sigma_yy[i] = sigma[1]
        sigma_zz[i] = sigma[2]

    return eps_xx, sigma_xx
