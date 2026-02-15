import numpy as np
import matplotlib.pyplot as plt


def drucker_prager_yield(sigma, cohesion, friction_angle):
    """Check Drucker-Prager yield condition."""
    I1 = np.trace(sigma)
    J2 = 0.5 * np.sum((sigma - np.eye(3) * I1 / 3) ** 2)
    return np.sqrt(J2) + friction_angle * I1 - cohesion


def voigt_map(stress):
    """Convert stress tensor to Voigt notation."""
    return np.array([stress[0, 0], stress[1, 1], stress[2, 2], stress[0, 1], stress[1, 2], stress[0, 2]])


def inverse_voigt_map(voigt_stress):
    """Convert Voigt notation back to stress tensor."""
    return np.array([[voigt_stress[0], voigt_stress[3], voigt_stress[5]],
                     [voigt_stress[3], voigt_stress[1], voigt_stress[4]],
                     [voigt_stress[5], voigt_stress[4], voigt_stress[2]]])


def triaxial_test(E, nu, cohesion, friction_angle, p_confining, axial_strains):
    """Simulates a triaxial test on an elastic-plastic material."""
    stress_history = []

    # Elastic stiffness matrix in Voigt notation
    C = E / ((1 + nu) * (1 - 2 * nu)) * np.array([
        [1 - nu, nu, nu, 0, 0, 0],
        [nu, 1 - nu, nu, 0, 0, 0],
        [nu, nu, 1 - nu, 0, 0, 0],
        [0, 0, 0, (1 - 2 * nu) / 2, 0, 0],
        [0, 0, 0, 0, (1 - 2 * nu) / 2, 0],
        [0, 0, 0, 0, 0, (1 - 2 * nu) / 2]
    ])

    sigma = np.array([[p_confining, 0, 0],
                      [0, p_confining, 0],
                      [0, 0, p_confining]])  # Initial confining stress

    for axial_strain in axial_strains:
        strain = np.array([axial_strain, -nu * axial_strain, -nu * axial_strain, 0, 0, 0])  # Voigt strain vector

        trial_stress_voigt = voigt_map(sigma) + C @ strain
        trial_stress = inverse_voigt_map(trial_stress_voigt)

        if drucker_prager_yield(trial_stress, cohesion, friction_angle) > 0:
            # Plastic correction using return mapping
            flow_direction = trial_stress - np.eye(3) * np.trace(trial_stress) / 3
            norm_flow = np.linalg.norm(flow_direction)

            if norm_flow > 1e-8:  # Avoid division by zero
                flow_direction /= norm_flow
                delta_lambda = drucker_prager_yield(trial_stress, cohesion, friction_angle) / E
                plastic_correction = delta_lambda * flow_direction
                trial_stress -= plastic_correction

        stress_history.append(trial_stress[0, 0])  # Axial stress
        sigma = trial_stress  # Update stress state

    return stress_history


# Parameters
E = 200e3  # Young's modulus (MPa)
nu = 0.3  # Poisson's ratio
cohesion = 20  # Cohesion (MPa)
friction_angle = 0.3  # Friction coefficient
p_confining = 50  # Confining pressure (MPa)
axial_strains = np.linspace(0, 0.02, 50)  # 2% axial strain

# Run simulation
axial_stresses = triaxial_test(E, nu, cohesion, friction_angle, p_confining, axial_strains)

# Plot results
plt.plot(axial_strains, axial_stresses, label='Triaxial Test')
plt.xlabel('Axial Strain')
plt.ylabel('Axial Stress (MPa)')
plt.title('Triaxial Test Simulation')
plt.legend()
plt.grid()
plt.show()
