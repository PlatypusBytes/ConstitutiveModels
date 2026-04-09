import numpy as np
from python_prototypes.hardening_soil import  HardeningSoilWithCap
from python_prototypes.stress_utils import StressUtils


def matrix_to_voigt(m):
    return np.array([m[0,0], m[1,1], m[2,2], m[1,2], m[0,2], m[0,1]])

def voigt_to_matrix(v):
    m = np.zeros((3,3))
    m[0,0] = v[0]
    m[1,1] = v[1]
    m[2,2] = v[2]
    m[1,2] = m[2,1] = v[3]
    m[0,2] = m[2,0] = v[4]
    m[0,1] = m[1,0] = v[5]
    return m


def example_triaxial_test():
    from tests.incr_driver import IncrDriver
    cell_pressure = 100.0

    sigma0 = np.array([cell_pressure,cell_pressure,cell_pressure,0,0,0])


    initial_stress = sigma0
    stress_increment = np.zeros(6)
    n_time_steps = 50
    max_iterations = 50
    voigt_size = 6
    ndim = 3



    G = 2000.
    nu = 0.33
    E = 2. * G * (1. + nu)


    phi = 34.0
    k0nc = 1 - np.sin(np.radians(phi))
    k0nc = 0.4
    M = 6 * (1 - k0nc) / (2 + k0nc)
    M2 = 9 / 2 * (1 - k0nc) / (1 + 2 * k0nc)
    M = 1.47
    #
    # M = 1/np.sqrt(3) * (6*np.sin(np.radians(phi)) / (3 - np.sin(np.radians(phi))))


    bulk_stiffness_ratio = 1.84
    params = {
        'E50_ref': 30000,
        'Eur_ref': 3 * 30000,
        'm': 0.55,
        'c': 0.0,
        'phi': 42.0,
        # 'phi': 36.0,
        'psi': 16,
        'pref': 100,
        'Rf': 0.9,
        'nu': 0.25,
        'M': M,  # cap shape parameter
        'k0nc': k0nc,
        'K_ratio': bulk_stiffness_ratio
    }

    model = HardeningSoilWithCap(params)

    p0 = StressUtils.p(sigma0)
    model.set_initial_state(sigma0, p_p0=p0)  # lightly overconsolidated

    deps_axial = 0.10 / n_time_steps  #

    # compute the initial elastic matrix
    ddsdde = model.ddsde(sigma0)

    # initialize stress and strain vectors
    strain_vector = np.zeros(voigt_size)
    stress_vector = np.copy(initial_stress)

    # stress_vector = np.copy(matrix_to_voigt(initial_stress))

    stresses = []
    strains = []
    stiffnesses = []

    strain_increment = np.array([deps_axial, 0, 0, 0.0, 0.0, 0.0])
    # strain_increment = np.array([deps_axial, 0, 0, 0.0, 0.0, 0.0])
    control_type = np.zeros(voigt_size)

    # triaxial strain controlled on axial direction, stress controlled on lateral directions
    control_type[:3] = 1
    control_type[0] = 0


    # loop over time steps
    for t in range(n_time_steps):
        print(f"Time step {t+1}/{n_time_steps}")

        if t >= n_time_steps - 1:
            print('test')

        delta_strain = np.copy(strain_increment)
        correction_delta_strain = np.zeros(6)
        approx_delta_stress = np.zeros(6)

        model.update_variables(model.sigma)

        # save the stress vector of the previous time step
        old_d_stress_vector = np.copy(stress_vector)

        # save the state variable vector of the previous time step
        prev_state_variables = np.copy([model.p_p, model.gamma_p])

        # loop over the maximum number of non-linear iterations
        for i in range(max_iterations + 1):

            # get stress increment at stress controlled components, else zero
            stress_increment = np.where(control_type, stress_increment, 0.0)

            # calculate the undesired stress
            u_d_stress = stress_increment + approx_delta_stress

            # update the delta strain using the elastic matrix and the undesired stress
            IncrDriver.calculate_iteration(correction_delta_strain, u_d_stress, control_type, ddsdde)

            # correct delta strain
            delta_strain = delta_strain - correction_delta_strain

            # run constitutive model
            model.integrate(np.copy(delta_strain))
            stress_updated = np.copy(model.sigma)
            ddsdde = model.ddsde(model.sigma)

            if i < max_iterations:

                if np.any(np.isnan(stress_updated)):
                    raise ValueError(f"NaN detected in stress at iteration {i} of time step {t}.")
                approx_delta_stress = stress_updated - old_d_stress_vector

                # reset stress vector and state variables
                stress_vector = np.copy(old_d_stress_vector)
                model.sigma = np.copy(stress_vector)
                state_variables = np.copy(prev_state_variables)
                model.p_p = state_variables[0]
                model.gamma_p = state_variables[1]
            else:
                strain_vector = strain_vector + delta_strain
                stress_vector = np.copy(stress_updated)
                model.sigma = np.copy(stress_vector)

                stiffness = ddsdde[0,0]  # axial stiffness
                stiffnesses.append(stiffness)

        stresses.append(stress_vector)
        np_stresses = np.array(stresses)

        strains.append(strain_vector)
        np_strains = np.array(strains)



    return np_strains, np_stresses, np.array(stiffnesses)

if __name__ == "__main__":



    import cProfile
    #
    # cProfile.run('example_triaxial_test()', "hs_profile")

    # import pstats

    # p = pstats.Stats("hs_profile")
    # p.strip_dirs().sort_stats("time").print_stats(20)


    #
    np_strains, np_stresses, np_stiffnesses= example_triaxial_test()

    q = np.sqrt(3/2 * ((np_stresses[:,0] - np_stresses[:,1])**2 + (np_stresses[:,1] - np_stresses[:,2])**2 + (np_stresses[:,2] - np_stresses[:,0])**2))
    p = (np_stresses[:,0] + np_stresses[:,1] + np_stresses[:,2]) / 3

    import matplotlib.pyplot as plt

    plt.figure(figsize=(12, 4))

    plt.subplot(1, 5, 1)
    # plt.plot(axial_strain, deviator_stress, 'b-')
    plt.plot(np_strains[:,0], np_stresses[:,0]/np_stresses[:,2], 'b-')
    # plt.ylim([1,6])
    plt.xlabel('Axial strain')
    plt.ylabel('sigma1/sigma3 (-)')
    plt.title('Stress‑strain')

    plt.subplot(1, 5, 2)
    plt.plot(p, q,marker='o')
    plt.xlabel('mean stress p (kPa)')
    plt.ylabel('Deviator stress q (kPa)')

    plt.subplot(1, 5, 3)
    plt.plot(np_stresses[:,0], label='sigma1')
    plt.plot(np_stresses[:,1], label='sigma2')
    plt.plot(np_stresses[:,2], label='sigma3')
    plt.xlabel('Time step')
    plt.ylabel('Stress (kPa)')
    plt.title('Stress components')
    plt.legend()

    plt.subplot(1, 5, 4)
    plt.plot(np_strains[:,0],label='eps1')
    plt.plot(np_strains[:,1],label='eps2')
    plt.plot(np_strains[:,2],label='eps3')
    plt.xlabel('Time step')
    plt.ylabel('Strain')
    plt.title('Strain components')
    plt.legend()

    plt.subplot(1, 5, 5)
    plt.plot(np_strains[:,0], np_stiffnesses, 'b-')
    plt.xlabel('Axial strain')
    plt.ylabel('Axial stiffness (kPa)')
    plt.title('Axial stiffness vs axial strain')

    plt.show()