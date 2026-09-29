"""
Run a drained triaxial compression test with the C Hardening Soil UMAT and
reproduce the same 5-panel figure as run_hardening_soil.py, but driven by the
compiled C model (build_C/lib/hardening_soil.dll).

Convention: compression positive, Voigt [xx, yy, zz, xy, yz, xz]. The UMAT interface is
tension positive (Abaqus convention); run_umat converts.
PROPS = [E50_ref, Eur_ref, m, c, phi_deg, psi_deg, p_ref, Rf, nu, M_cap, K_ratio]
STATEV = [gamma_p, p_p, at_failure]

The lateral (confining) stress is held constant with a small secant loop on the
lateral strain, which is far more robust than the fixed-point stress control and
lets the model harden all the way to Mohr-Coulomb failure.
"""

import os
import sys
import shutil
import subprocess

import numpy as np

from tests.utils import Utils


def dll_path():
    ext = "dll" if sys.platform == "win32" else "so"
    return os.path.join(os.getcwd(), "build_C", "lib", f"hardening_soil.{ext}")


def ensure_dll():
    path = dll_path()
    if os.path.exists(path):
        return path
    if shutil.which("gcc") is None:
        raise RuntimeError("hardening_soil library not built and gcc not available")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    root = os.getcwd()
    sources = [
        os.path.join(root, "c_models", "hardening_soil", "hardening_soil.c"),
        os.path.join(root, "c_models", "globals.c"),
        os.path.join(root, "c_models", "utils.c"),
        os.path.join(root, "c_models", "stress_utils.c"),
        os.path.join(root, "c_models", "elastic_laws", "hookes_law.c"),
    ]
    subprocess.run(["gcc", "-O2", "-shared", "-o", path, *sources, "-lm"], check=True)
    return path


def run_umat(dll, stress, statev, strain, dstrain, props, time_step):
    """UMAT call with compression-positive stresses and strains (the UMAT is tension positive)."""
    stress_new, ddsdde, statev_new = Utils.run_c_umat(dll, -stress, statev, -strain, -dstrain, props,
                                                      time_step)
    return -stress_new, ddsdde, statev_new


def run_triaxial_c(cell_pressure=600.0, total_axial_strain=0.10, n_time_steps=50):
    dll = ensure_dll()

    # E50_ref, Eur_ref, m, c, phi, psi, p_ref, Rf, nu, M_cap, K_ratio e0, ecv
    props = [30000, 3 * 30000, 0.55, 0.0, 42, 16, 100.0, 0.85, 0.25, 1.5, 1.84, 0.63, 1]

    stress = np.array([cell_pressure, cell_pressure, cell_pressure, 0.0, 0.0, 0.0])
    statev = np.array([0.0, 0.0, 0.0])
    strain = np.zeros(6)

    d_axial = total_axial_strain / n_time_steps

    stresses, strains, stiffnesses = [], [], []

    import time

    tim = time.time()

    e = 0.0  # lateral strain increment (warm-started between steps)
    for t in range(n_time_steps):
        # Secant iteration on the (symmetric) lateral strain so that the lateral
        # stress returns to `cell_pressure`. The first iteration uses the
        # derivative D[1,1] + D[1,2] of the tangent returned by the UMAT.
        e_prev, r_prev = None, None
        for _ in range(80):
            ds = np.array([d_axial, e, e, 0.0, 0.0, 0.0])
            s_try, ddsdde, statev_try = run_umat(dll, stress, statev.copy(), strain, ds, props, t)
            r = s_try[1] - cell_pressure
            if abs(r) < 1e-6:
                break
            if r_prev is not None and r != r_prev:
                deriv = (r - r_prev) / (e - e_prev)
            else:
                deriv = ddsdde[1, 1] + ddsdde[1, 2]
            if deriv < 1e-9:
                break
            e_prev, r_prev = e, r
            e -= r / deriv

        stress, statev = s_try, statev_try
        strain = strain + ds

        if np.any(np.isnan(stress)):
            raise ValueError(f"NaN detected in stress at time step {t}.")

        stresses.append(stress.copy())
        strains.append(strain.copy())
        stiffnesses.append(ddsdde[0, 0])

    print(f"Completed {n_time_steps} steps in {time.time() - tim:.2f} seconds.")

    return np.array(strains), np.array(stresses), np.array(stiffnesses)


if __name__ == "__main__":
    np_strains, np_stresses, np_stiffnesses = run_triaxial_c()

    q = np.sqrt(3 / 2 * ((np_stresses[:, 0] - np_stresses[:, 1]) ** 2 +
                         (np_stresses[:, 1] - np_stresses[:, 2]) ** 2 +
                         (np_stresses[:, 2] - np_stresses[:, 0]) ** 2))
    p = (np_stresses[:, 0] + np_stresses[:, 1] + np_stresses[:, 2]) / 3

    import matplotlib
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
    out = os.path.join(os.getcwd(), "hardening_soil_c_triaxial.png")
    plt.savefig(out, dpi=120)
    print(f"Saved figure to {out}")
    # Also show if a GUI backend is available.
    if matplotlib.get_backend().lower() not in ("agg", "pdf", "svg", "ps"):
        plt.show()



