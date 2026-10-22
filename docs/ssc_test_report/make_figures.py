"""
Generates the figures of the Soft Soil Creep test report by running the constant rate of strain
(CRS) tests of tests/test_soft_soil_creep.py:

    python docs/ssc_test_report/make_figures.py
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # docs
import figure_utils as fu  # noqa: E402
from figure_utils import ENVELOPE, LINESTYLES, SERIES  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, fu.ROOT)
from tests.test_soft_soil_creep import run_crs, p_q, M_MC  # noqa: E402

save = fu.figure_saver(os.path.join(HERE, "figures"))


def curves(history):
    eps = np.array([t for t, _, _ in history])
    pq = np.array([p_q(s) for _, s, _ in history])
    return eps, pq[:, 0], pq[:, 1]


def strain_q_and_p_q(cases, name, p0, envelope=False):
    """
    cases: list of (label, history) starting from the isotropic stress p0.
    Writes <name>_eps_q.pdf and <name>_p_q.pdf.
    """
    fig, ax = plt.subplots()
    for (label, h), color, ls in zip(cases, SERIES, LINESTYLES):
        eps, _, q = curves(h)
        ax.semilogy(100 * eps, q, color=color, ls=ls, label=label)
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [%]")
    ax.set_ylabel(r"deviator stress $q$ [kPa]")
    ax.legend()
    save(fig, f"{name}_eps_q")

    fig, ax = plt.subplots()
    p_max = 0.0
    for (label, h), color, ls in zip(cases, SERIES, LINESTYLES):
        _, p, q = curves(h)
        p, q = np.insert(p, 0, p0), np.insert(q, 0, 0.0)
        ax.plot(p, q, color=color, ls=ls, label=label)
        p_max = max(p_max, p.max())
    if envelope:
        pp = np.linspace(0, 1.05 * p_max, 2)
        ax.plot(pp, M_MC * pp, **ENVELOPE, label=rf"MC envelope $q = {M_MC:.2f}\,p$")
    ax.set_xlabel(r"mean effective stress $p'$ [kPa]")
    ax.set_ylabel(r"deviator stress $q$ [kPa]")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend()
    save(fig, f"{name}_p_q")


def main():
    fu.apply_style()
    dll = fu.library_path("soft_soil_creep")
    if not os.path.exists(dll):
        sys.exit(f"shared library not found: {dll} (run the pytest suite once to build it)")

    # Fig. 2 of the paper: CRS oedometer, OCR0 = 6, p0 = 10 kPa
    strain_q_and_p_q([(r"100 steps ($\Delta t$ = 21.6 min)", run_crs(dll, "oedometer", 10.0, 6.0, 100)),
                      (r"1000 steps ($\Delta t$ = 2.16 min)", run_crs(dll, "oedometer", 10.0, 6.0, 1000))],
                     "oedometer", 10.0)

    # Fig. 3 of the paper: CRS undrained triaxial, NC, p0 = 100 kPa
    strain_q_and_p_q([(r"100 steps ($\Delta t$ = 21.6 min)", run_crs(dll, "undrained", 100.0, 1.0, 100)),
                      (r"1000 steps ($\Delta t$ = 2.16 min)", run_crs(dll, "undrained", 100.0, 1.0, 1000))],
                     "undrained", 100.0, envelope=True)

    # rate dependence of the undrained test
    strain_q_and_p_q([(f"$t_{{end}}$ = {T:g} min", run_crs(dll, "undrained", 100.0, 1.0, 100, total_time=T))
                      for T in (21600.0, 2160.0, 216.0)],
                     "rate", 100.0, envelope=True)

    # key numbers for the report
    h = run_crs(dll, "oedometer", 10.0, 6.0, 100)
    s, sv = h[-1][1], h[-1][2]
    print(f"oedometer: K0 = {s[1] / s[0]:.3f}, eps_c = {sv[0]:.4f}, OCR = {sv[3]:.3f}")
    for T in (21600.0, 2160.0, 216.0):
        _, p, q = curves(run_crs(dll, "undrained", 100.0, 1.0, 100, total_time=T))
        i = np.argmax(q)
        print(f"undrained t_end={T:g}: q_peak = {q[i]:.2f} kPa at p = {p[i]:.2f}, end q = {q[-1]:.2f}, p = {p[-1]:.2f}")


if __name__ == "__main__":
    main()
