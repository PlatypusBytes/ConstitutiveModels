"""
Drained triaxial compression test with the Python Hardening Soil prototype. Same test, parameters
and figure as run_hardening_soil_c.py, which drives the compiled C model.

Convention: compression positive, Voigt [xx, yy, zz, xy, yz, xz].

The lateral (confining) stress is held constant with a secant loop on the lateral strain.
"""

import time

import numpy as np

from python_prototypes.hardening_soil import HardeningSoil


def trial_increment(model, ds):
    """
    Integrates ds from the current state without committing it.

    :return: converged flag and the new state (stress, state variables, volumetric strain, tangent)
    """
    start = (model.sigma.copy(), model.state_variables, model.eps_v)
    converged = model.integrate(ds)
    state = (model.sigma.copy(), model.state_variables, model.eps_v, model.ddsdde.copy())
    model.sigma, model.state_variables, model.eps_v = start
    return converged, state


def commit_increment(model, state):
    model.sigma, model.state_variables, model.eps_v, model.ddsdde = state


def run_triaxial(cell_pressure=600.0, total_axial_strain=0.10, n_time_steps=50):
    # E50_ref, Eur_ref, m, c, phi, psi, p_ref, Rf, nu, M_cap, K_ratio, e0, e_cv
    props = [30000, 3 * 30000, 0.55, 0.0, 42, 16, 100.0, 0.85, 0.25, 1.5, 1.84, 0.63, 1]
    model = HardeningSoil.from_props(props)
    model.set_initial_state(np.array([cell_pressure, cell_pressure, cell_pressure, 0.0, 0.0, 0.0]))

    strain = np.zeros(6)
    d_axial = total_axial_strain / n_time_steps

    stresses, strains, stiffnesses = [], [], []
    tim = time.time()

    e = 0.0  # lateral strain increment (warm-started between steps)
    for t in range(n_time_steps):
        # Secant iteration on the (symmetric) lateral strain so that the lateral stress returns to
        # `cell_pressure`. The first iteration uses the derivative D[1,1] + D[1,2] of the tangent.
        e_prev, r_prev = None, None
        for _ in range(80):
            ds = np.array([d_axial, e, e, 0.0, 0.0, 0.0])
            converged, state = trial_increment(model, ds)
            r = state[0][1] - cell_pressure
            if abs(r) < 1e-6:
                break
            if r_prev is not None and r != r_prev:
                deriv = (r - r_prev) / (e - e_prev)
            else:
                deriv = state[3][1, 1] + state[3][1, 2]
            if deriv < 1e-9:
                break
            e_prev, r_prev = e, r
            e -= r / deriv

        if not converged:
            raise RuntimeError(f"Integration did not converge at time step {t}.")
        commit_increment(model, state)
        strain = strain + ds

        if np.any(np.isnan(model.sigma)):
            raise ValueError(f"NaN detected in stress at time step {t}.")

        stresses.append(model.sigma.copy())
        strains.append(strain.copy())
        stiffnesses.append(model.ddsdde[0, 0])

    print(f"Completed {n_time_steps} steps in {time.time() - tim:.2f} seconds.")

    return np.array(strains), np.array(stresses), np.array(stiffnesses)


if __name__ == "__main__":
    np_strains, np_stresses, np_stiffnesses = run_triaxial()

    q = np.sqrt(3 / 2 * ((np_stresses[:, 0] - np_stresses[:, 1]) ** 2 +
                         (np_stresses[:, 1] - np_stresses[:, 2]) ** 2 +
                         (np_stresses[:, 2] - np_stresses[:, 0]) ** 2))
    p = (np_stresses[:, 0] + np_stresses[:, 1] + np_stresses[:, 2]) / 3

    import matplotlib.pyplot as plt

    plt.figure(figsize=(16, 4))

    plt.subplot(1, 5, 1)
    plt.plot(np_strains[:, 0], np_stresses[:, 0] / np_stresses[:, 2], 'b-')
    plt.xlabel('Axial strain')
    plt.ylabel('sigma1/sigma3 (-)')
    plt.title('Stress-strain')

    plt.subplot(1, 5, 2)
    plt.plot(p, q, marker='o')
    plt.xlabel('mean stress p (kPa)')
    plt.ylabel('Deviator stress q (kPa)')
    plt.title('p-q path')

    plt.subplot(1, 5, 3)
    plt.plot(np_stresses[:, 0], label='sigma1')
    plt.plot(np_stresses[:, 1], label='sigma2')
    plt.plot(np_stresses[:, 2], label='sigma3')
    plt.xlabel('Time step')
    plt.ylabel('Stress (kPa)')
    plt.title('Stress components')
    plt.legend()

    plt.subplot(1, 5, 4)
    plt.plot(np_strains[:, 0], label='eps1')
    plt.plot(np_strains[:, 1], label='eps2')
    plt.plot(np_strains[:, 2], label='eps3')
    plt.xlabel('Time step')
    plt.ylabel('Strain')
    plt.title('Strain components')
    plt.legend()

    plt.subplot(1, 5, 5)
    plt.plot(np_strains[:, 0], np_stiffnesses, 'b-')
    plt.xlabel('Axial strain')
    plt.ylabel('Axial stiffness (kPa)')
    plt.title('Axial stiffness vs axial strain')

    plt.tight_layout()
    plt.show()
