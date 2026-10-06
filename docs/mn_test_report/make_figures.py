"""
Figures for the Matsuoka-Nakai test report (docs/mn_test_report/mn_test_report.tex).

Drives the Matsuoka-Nakai UMAT (build_C/lib/matsuoka_nakai.dll or .so) along element test paths
and compares the results with closed-form solutions of the Matsuoka-Nakai criterion
(Matsuoka & Nakai, 1974) and of the plastic flow at failure.

    python docs/mn_test_report/make_figures.py

Stresses and strains are compression positive, as in the tests; the UMAT interface is tension
positive and `step` converts. The figures are written to docs/mn_test_report/figures and the values
quoted in the report, including the verification summary (Table 1), are printed.
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # docs
import figure_utils as fu  # noqa: E402
from figure_utils import AXIS, ENVELOPE, INK2, LINESTYLES, REF, SERIES, History, marker, p_q, summary  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

FIG_DIR = os.path.join(HERE, "figures")
DLL = fu.library_path("matsuoka_nakai")
UMAT = fu.Umat(DLL, "MN")
save = fu.figure_saver(FIG_DIR)


def step(stress, dstrain, props, umat=None):
    """UMAT call in the compression-positive convention (the tangent is the same in both); the
    single state variable is the return type."""
    s, statev, ddsdde, _, _ = (umat or UMAT)(-np.asarray(stress, float), [0.0], np.zeros(6),
                                             -np.asarray(dstrain, float), props)
    return -s, ddsdde, float(statev[0])


# --------------------------------------------------------------------------- #
# Parameters                                                                    #
# --------------------------------------------------------------------------- #
BASE = dict(E=20000.0, nu=0.3, c=5.0, phi=30.0, psi=10.0)


def props_of(**overrides):
    v = dict(BASE, **overrides)
    return [v["E"], v["nu"], v["c"], v["phi"], v["psi"]]


def attraction(c, phi):
    return c / np.tan(np.radians(phi))


# --------------------------------------------------------------------------- #
# Closed-form solutions                                                         #
# --------------------------------------------------------------------------- #
def tensor(s):
    return np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])


def voigt(t):
    return np.array([t[0, 0], t[1, 1], t[2, 2], t[0, 1], t[1, 2], t[0, 2]])


def principal(s):
    """Principal stresses sigma_1 >= sigma_2 >= sigma_3 (compression positive)."""
    return np.sort(np.linalg.eigvalsh(tensor(s)))[::-1]


def mn_k(phi):
    return 9.0 + 8.0 * np.tan(np.radians(phi)) ** 2


def mn_ratio(sig):
    """I1 I2 / I3 of (shifted, compression positive) principal stresses."""
    sig = np.asarray(sig, float)
    return sig.sum() * (sig[0] * sig[1] + sig[1] * sig[2] + sig[2] * sig[0]) / sig.prod()


def mn_criterion(s, c, phi):
    """I1 I2 / I3 / k - 1 of the shifted principal stresses: zero on the Matsuoka-Nakai surface."""
    return mn_ratio(principal(s) + attraction(c, phi)) / mn_k(phi) - 1.0


def bisect(fun, lo, hi, n=200):
    f_lo = fun(lo)
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        f_mid = fun(mid)
        if np.sign(f_mid) == np.sign(f_lo):
            lo, f_lo = mid, f_mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def mn_strength_ratio(phi, b):
    """sigma_1* / sigma_3* at failure for a given b = (sigma_2 - sigma_3) / (sigma_1 - sigma_3)."""
    return bisect(lambda r: mn_ratio([r, 1.0 + b * (r - 1.0), 1.0]) - mn_k(phi), 1.0 + 1e-12, 1e3)


def sin_mobilised(s1, s3, a):
    return (s1 - s3) / (s1 + s3 + 2.0 * a)


def lagioia(stress_tension, angle, c=0.0):
    """f = M p - K + J alpha cos(acos(beta sin 3 theta) / 3), tension positive (as the C code), with
    the closed-form constants; for a plastic potential c does not matter."""
    s = np.sin(np.radians(angle))
    M = 2 * np.sqrt(3) * s / (3 - s)
    K = 2 * np.sqrt(3) * c * np.cos(np.radians(angle)) / (3 - s)
    alpha = 2 * np.sqrt(3 + s * s) / (3 - s)
    beta = s * (9 - s * s) / (3 + s * s) ** 1.5
    p = stress_tension[:3].mean()
    dev = tensor(stress_tension) - p * np.eye(3)
    j2 = 0.5 * np.sum(dev * dev)
    if j2 <= 0.0:
        return M * p - K
    xi = np.clip(1.5 * np.sqrt(3) * np.linalg.det(dev) / j2 ** 1.5, -1, 1)
    return M * p - K + np.sqrt(j2) * alpha * np.cos(np.arccos(beta * xi) / 3)


def plane_strain_failure(phi, psi, c, sigma3):
    """Stationary failure state in plane strain with sigma_3 given: on the yield surface, with
    d eps_2^p = 0, i.e. dg/dsigma_2 = 0 (solved by Newton iteration on sigma_1 and sigma_2)."""
    a = attraction(c, phi)

    def residual(x):
        s1, s2 = x
        st = -np.array([s1, s2, sigma3, 0, 0, 0])
        h = 1e-6 * (abs(s1) + 1)
        dg2 = (lagioia(st + [0, h, 0, 0, 0, 0], psi) - lagioia(st - [0, h, 0, 0, 0, 0], psi)) / (2 * h)
        return np.array([mn_ratio([s1 + a, s2 + a, sigma3 + a]) / mn_k(phi) - 1.0, dg2])

    x = np.array([3.5 * sigma3, 2.0 * sigma3])
    for _ in range(50):
        r = residual(x)
        jac = np.zeros((2, 2))
        for j in range(2):
            dx = np.zeros(2)
            dx[j] = 1e-6 * x[j]
            jac[:, j] = (residual(x + dx) - r) / dx[j]
        x = x - np.linalg.solve(jac, r)
        if np.max(np.abs(r)) < 1e-14:
            break
    return x


# --------------------------------------------------------------------------- #
# Element test paths                                                            #
# --------------------------------------------------------------------------- #
def mixed_step(stress, strain, dstrain, groups, targets, props, tol=1e-10):
    """Strain increment in which each group of components (sharing one strain increment) is
    adjusted by Newton iteration with DDSDDE so that the stress of its first component equals the
    target (as in tests/test_matsuoka_nakai.py)."""
    d = fu.mixed_strain_increment(lambda d: step(stress, d, props)[:2], dstrain, groups, targets, tol)
    s, _, state = step(stress, d, props)
    return s, state, strain + d


def drained_triaxial(props, cell, axial_strain, n_steps, extension=False):
    """Drained triaxial compression along x (sigma_y = sigma_z = cell) or extension along z
    (sigma_x = sigma_y = cell)."""
    stress, strain = np.array([cell, cell, cell, 0, 0, 0], float), np.zeros(6)
    if extension:
        d, groups = [0, 0, -axial_strain / n_steps, 0, 0, 0], [[0, 1]]
    else:
        d, groups = [axial_strain / n_steps, 0, 0, 0, 0, 0], [[1, 2]]
    h = History(strain, stress, 0.0)
    for _ in range(n_steps):
        stress, state, strain = mixed_step(stress, strain, d, groups, np.array([cell]), props)
        h.add(strain, stress, state)
    return h.arrays()


def plane_strain(props, cell, axial_strain, n_steps):
    """eps_y = 0, sigma_z = cell, compression along x."""
    stress, strain = np.array([cell, cell, cell, 0, 0, 0], float), np.zeros(6)
    h = History(strain, stress, 0.0)
    for _ in range(n_steps):
        stress, state, strain = mixed_step(stress, strain, [axial_strain / n_steps, 0, 0, 0, 0, 0], [[2]],
                                           np.array([cell]), props)
        h.add(strain, stress, state)
    return h.arrays()


def undrained_triaxial(props, p0, axial_strain, n_steps):
    stress, strain = np.array([p0, p0, p0, 0, 0, 0], float), np.zeros(6)
    d = axial_strain / n_steps
    dstrain = np.array([d, -d / 2, -d / 2, 0, 0, 0])
    h = History(strain, stress, 0.0)
    for _ in range(n_steps):
        stress, _, state = step(stress, dstrain, props)
        strain = strain + dstrain
        h.add(strain, stress, state)
    return h.arrays()


def lode_direction(lode):
    """Deviatoric principal direction: lode = -pi/6 is triaxial compression along x, +pi/6
    triaxial extension along z (compression positive)."""
    return np.array([np.sin(lode + 2 * np.pi / 3), np.sin(lode), np.sin(lode - 2 * np.pi / 3)])


def lode_sweep(props, p0=100.0, magnitude=0.02, n=49, umat=None):
    """Single large deviatoric strain increments in all directions of the deviatoric plane."""
    out = []
    for lode in np.linspace(-np.pi, np.pi, n):  # the full plane, including the other sextants
        s, _, state = step([p0, p0, p0, 0, 0, 0], magnitude * np.r_[lode_direction(lode), 0, 0, 0], props, umat)
        out.append((s, state))
    return out


# --------------------------------------------------------------------------- #
# Figures                                                                       #
# --------------------------------------------------------------------------- #
def fig_yield_surface():
    """Deviatoric section and the strength as a function of b, from returns of single large
    deviatoric increments in all directions."""
    c = BASE["c"]

    # (a) deviatoric plane at phi = 30, stresses shifted by a and normalised by p + a
    phi = BASE["phi"]
    a = attraction(c, phi)
    sweep = lode_sweep(props_of())

    angles = np.linspace(0, 2 * np.pi, 721)
    mn_curve = []
    for w in angles:
        e = np.array([np.cos(w), np.cos(w - 2 * np.pi / 3), np.cos(w + 2 * np.pi / 3)]) * np.sqrt(2 / 3)
        rho = bisect(lambda r: mn_ratio(1.0 + r * e) - mn_k(phi), 1e-9, (1 - 1e-12) / -e.min())
        mn_curve.append(1.0 + rho * e)
    mx, my = fu.pi_plane(np.array(mn_curve))
    cx, cy = fu.mohr_coulomb_section(phi)

    pts, worst = [], 0.0
    for s, state in sweep:
        sig = np.diag(tensor(s)) + a  # principal axes are x, y, z for these paths
        pts.append(sig / sig.mean())
        worst = max(worst, abs(mn_criterion(s, c, phi)))
    px, py = fu.pi_plane(np.array(pts))
    summary(f"\n[yield surface] phi = {phi}: {len(sweep)} directions, states {sorted(set(st for _, st in sweep))},"
            f" max |I1 I2 / (I3 k) - 1| = {worst:.1e}")

    fig, ax = plt.subplots()
    fu.deviatoric_axes(ax, 1.2 * np.hypot(cx, cy).max())
    ax.plot(cx, cy, **ENVELOPE, label="Mohr-Coulomb")
    ax.plot(mx, my, **REF, label=r"Matsuoka-Nakai, $I_1 I_2 / I_3 = 9 + 8\tan^2\varphi$")
    ax.plot(px, py, **marker(SERIES[0], 4.5), label="UMAT, returned stress")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=1, fontsize=7)
    save(fig, "deviatoric_section")

    # (b) mobilised friction angle at failure against b, three friction angles
    fig, ax = plt.subplots()
    b_line = np.linspace(0, 1, 101)
    summary("[strength vs b] phi, max |phi_m(UMAT) - phi_m(MN)| [deg], phi_m at b = 0.5 (MN), plane strain"
            " (associated) phi_ps")
    for phi_i, color in zip([20.0, 30.0, 40.0], SERIES):
        a_i = attraction(c, phi_i)
        mn_line = [np.degrees(np.arcsin((r - 1) / (r + 1))) for r in (mn_strength_ratio(phi_i, b) for b in b_line)]
        ax.plot(b_line, mn_line, **REF)
        ax.axhline(phi_i, **ENVELOPE)
        b_pts, phi_pts, dev = [], [], 0.0
        for lode in np.linspace(-np.pi / 6, np.pi / 6, 13):  # one sextant: b from 0 to 1
            s, _, _ = step([100.0, 100.0, 100.0, 0, 0, 0], 0.02 * np.r_[lode_direction(lode), 0, 0, 0],
                           props_of(phi=phi_i, psi=phi_i / 3))
            s1, s2, s3 = principal(s)
            b = (s2 - s3) / (s1 - s3)
            phi_m = np.degrees(np.arcsin(sin_mobilised(s1, s3, a_i)))
            r_mn = mn_strength_ratio(phi_i, b)
            dev = max(dev, abs(phi_m - np.degrees(np.arcsin((r_mn - 1) / (r_mn + 1)))))
            b_pts.append(b)
            phi_pts.append(phi_m)
        ax.plot(b_pts, phi_pts, **marker(color, 4.5), label=rf"UMAT, $\varphi$ = {phi_i:.0f}$^\circ$")
        k = mn_k(phi_i)
        t = 0.5 * ((np.sqrt(k) - 1) + np.sqrt((np.sqrt(k) - 1) ** 2 - 4))
        r_mid = mn_strength_ratio(phi_i, 0.5)
        summary(f"  {phi_i:4.0f} {dev:.1e} {np.degrees(np.arcsin((r_mid - 1) / (r_mid + 1))):.2f}"
                f" {np.degrees(np.arcsin((t * t - 1) / (t * t + 1))):.2f} (b = {1 / (t + 1):.3f})")
    ax.plot([], [], **REF, label="Matsuoka-Nakai")
    ax.plot([], [], **ENVELOPE, label="Mohr-Coulomb")
    ax.set_xlabel(r"$b = (\sigma_2 - \sigma_3)/(\sigma_1 - \sigma_3)$ [-]")
    ax.set_ylabel(r"$\varphi_m$ at failure [$^\circ$]")
    ax.set_xlim(0, 1)
    ax.set_ylim(15, 57)
    ax.legend(loc="upper center", ncol=3, columnspacing=0.8, fontsize=6.8, handlelength=1.8)
    save(fig, "strength_b")


def fig_triaxial():
    """Drained triaxial compression and extension at sigma_3 = 100 kPa."""
    c, phi = BASE["c"], BASE["phi"]
    a = attraction(c, phi)
    sin_phi = np.sin(np.radians(phi))
    n_f = (1 + sin_phi) / (1 - sin_phi)
    cell = 100.0
    runs = {(ext, psi): drained_triaxial(props_of(psi=psi), cell, 0.03, 60, ext)
            for ext in (False, True) for psi in (0.0, 10.0, 20.0)}

    summary("\n[triaxial] path, psi, (s1+a)/(s3+a) / N_phi - 1, -dev/dgamma UMAT, 4 sin(psi)/(3 -/+ sin(psi))")
    for (ext, psi), h in runs.items():
        s1, s2, s3 = principal(h.sig[-1])
        ev = h.eps[:, :3].sum(axis=1)
        gamma = h.eps[:, 0] - (h.eps[:, 2] if ext else h.eps[:, 1])
        ratio = -(ev[-1] - ev[-11]) / (gamma[-1] - gamma[-11])
        sp = np.sin(np.radians(psi))
        summary(f"  {'extension  ' if ext else 'compression'} {psi:4.1f} {(s1 + a) / (s3 + a) / n_f - 1:+.1e}"
                f" {ratio:.8f} {4 * sp / ((3 + sp) if ext else (3 - sp)):.8f}")

    # (a) deviator stress against shear strain (psi does not change it)
    fig, ax = plt.subplots()
    for ext, color, ls, name in [(False, SERIES[0], "-", "compression"), (True, SERIES[1], (0, (6, 2)), "extension")]:
        h = runs[(ext, 10.0)]
        gamma = h.eps[:, 0] - (h.eps[:, 2] if ext else h.eps[:, 1])
        ax.plot(100 * gamma, [principal(s)[0] - principal(s)[2] for s in h.sig], color=color, ls=ls,
                label=f"UMAT, triaxial {name}")
    q_c = (cell + a) * n_f - a - cell
    q_e = cell - ((cell + a) / n_f - a)
    ax.axhline(q_c, **ENVELOPE)
    ax.axhline(q_e, **ENVELOPE)
    ax.plot([], [], **ENVELOPE, label=r"$q_f$, $(\sigma_1 + a)/(\sigma_3 + a) = N_\varphi$")
    ax.set_xlabel(r"shear strain $\varepsilon_1 - \varepsilon_3$ [%]")
    ax.set_ylabel(r"$q = \sigma_1 - \sigma_3$ [kPa]")
    ax.set_xlim(0, 4.5)
    ax.set_ylim(0, 260)
    ax.legend(loc="lower right")
    save(fig, "triaxial_q")
    summary(f"  q_f compression {q_c:.2f} kPa, extension {q_e:.2f} kPa")

    # (b) volumetric strain in compression for three dilation angles
    fig, ax = plt.subplots()
    for psi, color, ls in zip((0.0, 10.0, 20.0), SERIES, LINESTYLES):
        h = runs[(False, psi)]
        gamma = h.eps[:, 0] - h.eps[:, 1]
        ax.plot(100 * gamma, 100 * h.eps[:, :3].sum(axis=1), color=color, ls=ls,
                label=rf"$\psi$ = {psi:.0f}$^\circ$")
    ax.axhline(0.0, color=AXIS, lw=0.8)
    ax.set_xlabel(r"shear strain $\varepsilon_1 - \varepsilon_3$ [%]")
    ax.set_ylabel(r"volumetric strain $\varepsilon_v$ [%] (compression +)")
    ax.set_xlim(0, 4.5)
    ax.legend(loc="lower left")
    save(fig, "triaxial_ev")

    # (c) dilatancy ratio at failure against psi
    psis = np.arange(-20.0, 31.0, 5.0)
    fig, ax = plt.subplots()
    psi_line = np.linspace(-20, 30, 100)
    sp_line = np.sin(np.radians(psi_line))
    ax.plot(psi_line, 4 * sp_line / (3 - sp_line), **REF, label=r"$4\sin\psi/(3 - \sin\psi)$")
    ax.plot(psi_line, 4 * sp_line / (3 + sp_line), **ENVELOPE, label=r"$4\sin\psi/(3 + \sin\psi)$")
    worst = 0.0
    for ext, color, name in [(False, SERIES[0], "compression"), (True, SERIES[1], "extension")]:
        ratios = []
        for psi in psis:
            h = drained_triaxial(props_of(psi=psi), cell, 0.03, 30, ext)
            ev = h.eps[:, :3].sum(axis=1)
            gamma = h.eps[:, 0] - (h.eps[:, 2] if ext else h.eps[:, 1])
            ratios.append(-(ev[-1] - ev[-6]) / (gamma[-1] - gamma[-6]))
            sp = np.sin(np.radians(psi))
            worst = max(worst, abs(ratios[-1] - 4 * sp / ((3 + sp) if ext else (3 - sp))))
        ax.plot(psis, ratios, **marker(color, 5), label=f"UMAT, {name}")
    summary(f"[dilatancy] max deviation of -dev/dgamma from the closed form, psi = -20..30: {worst:.1e}")
    ax.axhline(0.0, color=AXIS, lw=0.8)
    ax.axvline(0.0, color=AXIS, lw=0.8)
    ax.set_xlabel(r"dilation angle $\psi$ [$^\circ$]")
    ax.set_ylabel(r"$-\Delta\varepsilon_v / \Delta(\varepsilon_1 - \varepsilon_3)$ at failure [-]")
    ax.set_xlim(-21, 31)
    ax.set_ylim(-0.6, 0.85)
    ax.legend(loc="upper left")
    save(fig, "triaxial_dilatancy")


def fig_plane_strain():
    """Plane strain: the intermediate principal stress at failure follows from d eps_2^p = 0."""
    c, phi = BASE["c"], BASE["phi"]
    a = attraction(c, phi)
    cell = 100.0
    fig_b, ax_b = plt.subplots()
    fig_m, ax_m = plt.subplots()
    summary("\n[plane strain] psi, b end (UMAT), b stationary, phi_m end (UMAT), phi_m stationary,"
            " (s2+a)/sqrt((s1+a)(s3+a)) - 1, s2/((s1+s3)/2) - 1")
    for psi, color, ls in zip((30.0, 10.0, 0.0), SERIES, LINESTYLES):
        h = plane_strain(props_of(psi=psi), cell, 0.30, 120)
        s1, s2, s3 = h.sig[:, 0], h.sig[:, 1], h.sig[:, 2]
        with np.errstate(invalid="ignore", divide="ignore"):
            b = (s2 - s3) / (s1 - s3)
        phi_m = np.degrees(np.arcsin(sin_mobilised(s1, s3, a)))
        x1, x2 = plane_strain_failure(phi, psi, c, cell)
        b_stat = (x2 - cell) / (x1 - cell)
        phi_stat = np.degrees(np.arcsin(sin_mobilised(x1, cell, a)))
        label = (r"$\psi = \varphi$ = " if psi == phi else r"$\psi$ = ") + f"{psi:.0f}" + r"$^\circ$"
        ax_b.plot(100 * h.eps[1:, 0], b[1:], color=color, ls=ls, label=label)
        ax_b.axhline(b_stat, xmin=0.75, **REF)
        ax_m.plot(100 * h.eps[:, 0], phi_m, color=color, ls=ls, label=label)
        ax_m.axhline(phi_stat, xmin=0.75, **REF)
        summary(f"  {psi:4.1f} {b[-1]:.6f} {b_stat:.6f} {phi_m[-1]:.4f} {phi_stat:.4f}"
                f" {(s2[-1] + a) / np.sqrt((s1[-1] + a) * (s3[-1] + a)) - 1:+.1e} {s2[-1] / (0.5 * (s1[-1] + s3[-1])) - 1:+.1e}")
    ax_b.plot([], [], **REF, label=r"stationary, $\partial g / \partial \sigma_2 = 0$")
    ax_b.set_xlabel(r"axial strain $\varepsilon_1$ [%]")
    ax_b.set_ylabel(r"$b = (\sigma_2 - \sigma_3)/(\sigma_1 - \sigma_3)$ [-]")
    ax_b.set_xlim(0, 30)
    ax_b.set_ylim(0, 0.7)
    ax_b.legend(loc="lower right")
    save(fig_b, "plane_strain_b")
    ax_m.axhline(phi, **ENVELOPE)
    ax_m.annotate(r"triaxial compression, $\varphi$ = %.0f$^\circ$" % phi, (29.5, phi), xytext=(0, 3),
                  textcoords="offset points", ha="right", va="bottom", color=INK2, fontsize=7)
    ax_m.plot([], [], **REF, label="stationary")
    ax_m.set_xlabel(r"axial strain $\varepsilon_1$ [%]")
    ax_m.set_ylabel(r"mobilised friction angle $\varphi_m$ [$^\circ$]")
    ax_m.set_xlim(0, 30)
    ax_m.set_ylim(29, 35)
    ax_m.legend(loc="lower right", bbox_to_anchor=(1.0, 0.25))
    save(fig_m, "plane_strain_phi")


def fig_undrained():
    """Undrained (isochoric) triaxial compression from p = 100 kPa."""
    c, phi = BASE["c"], BASE["phi"]
    a = attraction(c, phi)
    sin_phi = np.sin(np.radians(phi))
    m_q = 6 * sin_phi / (3 - sin_phi)
    fig_pq, ax_pq = plt.subplots()
    fig_q, ax_q = plt.subplots()
    summary("\n[undrained] psi, p' range (psi = 0) or p' at the end, q/(M (p+a)) - 1 on the cone,"
            " eps_1 at which the apex is reached")
    for psi, color, ls in zip((-5.0, 0.0, 10.0), SERIES, LINESTYLES):
        h = undrained_triaxial(props_of(psi=psi), 100.0, 0.05, 100)
        p, q = p_q(h.sig)
        label = rf"$\psi = {psi:.0f}^\circ$"
        ax_pq.plot(p, q, color=color, ls=ls, label=label)
        ax_q.plot(100 * h.eps[:, 0], q, color=color, ls=ls, label=label)
        extra = f"{p.min():.10f} .. {p.max():.10f}" if psi == 0 else f"{p[-1]:.2f}"
        cone = h.sv == 1
        dev = np.max(np.abs(q[cone] / (m_q * (p[cone] + a)) - 1)) if cone.any() else np.nan
        apex = h.eps[np.argmax(h.sv == 2), 0] if (h.sv == 2).any() else np.nan
        summary(f"  {psi:4.1f} {extra} {dev:.1e} {apex:.4f}")
    pl = np.array([-a, 450.0])
    ax_pq.plot(pl, m_q * (pl + a), **ENVELOPE, label=rf"$q = {m_q:.1f}\,(p' + a)$")
    ax_pq.axvline(0.0, color=AXIS, lw=0.8)
    ax_pq.set_xlabel(r"mean effective stress $p'$ [kPa]")
    ax_pq.set_ylabel(r"deviator stress $q$ [kPa]")
    ax_pq.set_xlim(-20, 450)
    ax_pq.set_ylim(0, 600)
    ax_pq.legend(loc="upper left")
    save(fig_pq, "undrained_pq")
    ax_q.set_xlabel(r"axial strain $\varepsilon_1$ [%]")
    ax_q.set_ylabel(r"deviator stress $q$ [kPa]")
    ax_q.set_xlim(0, 5)
    ax_q.set_ylim(0, 600)
    ax_q.legend(loc="upper left")
    save(fig_q, "undrained_q")


def global_newton(umat, tangent="consistent", max_iter=15):
    """A mixed-control strain increment from a state at failure, solved as in a finite element
    code: eps_x and eps_y prescribed, sigma_z = 100 kPa and sigma_xy = 5 kPa, sigma_yz = sigma_xz = 0
    as targets. Returns the residual norm per iteration."""
    props = props_of()
    stress = np.array([100.0, 100.0, 100.0, 0, 0, 0])
    for _ in range(8):  # plane strain loading up to failure
        stress, _, _ = step(stress, [2.5e-3, 0, 1e-3, 0, 0, 0], props, umat)
    free = [2, 3, 4, 5]
    target = np.array([100.0, 5.0, 0.0, 0.0])
    de = np.array([2.5e-3, 0, 0, 0, 0, 0])
    residuals = []
    D_el = fu.elastic_stiffness(BASE["E"], BASE["nu"])
    for _ in range(max_iter):
        s, ddsdde, _ = step(stress, de, props, umat)
        r = s[free] - target
        residuals.append(np.linalg.norm(r))
        if residuals[-1] < 1e-12 * 100.0:
            break
        D = ddsdde if tangent == "consistent" else D_el
        de[free] -= np.linalg.solve(D[np.ix_(free, free)], r)
    return np.array(residuals)


def fig_numerics():
    # (a) global Newton iteration
    fig, ax = plt.subplots()
    summary("\n[global Newton] tangent: residual per iteration [kPa]")
    for tangent, color, ls in [("consistent", SERIES[0], "-"), ("elastic", SERIES[1], (0, (6, 2)))]:
        res = global_newton(UMAT, tangent)
        summary(f"  {tangent:10s} " + " ".join(f"{v:.1e}" for v in res))
        ax.semilogy(np.arange(1, len(res) + 1), res, color=color, ls=ls, marker="o", ms=4, mec="white", mew=0.8,
                    label=f"DDSDDE" if tangent == "consistent" else "elastic stiffness")
    ax.axhline(1e-10, **ENVELOPE)
    ax.set_xlabel("global iteration [-]")
    ax.set_ylabel(r"residual $\|\Delta\sigma\|$ [kPa]")
    ax.set_xlim(0.5, 15.5)
    ax.set_ylim(1e-14, 1e3)
    ax.legend(loc="center right")
    save(fig, "global_newton")

    # (b) step size, plane strain with psi = 10
    props = props_of()
    fine = plane_strain(props, 100.0, 0.10, 1000)
    b_fine = (fine.sig[1:, 1] - fine.sig[1:, 2]) / (fine.sig[1:, 0] - fine.sig[1:, 2])
    fig, ax = plt.subplots()
    ax.plot(100 * fine.eps[1:, 0], b_fine, color=INK2, lw=0.9, label="1000 steps (reference)")
    summary("[step size] plane strain psi = 10, eps_1 = 10 %: n_steps, max |b - b_ref|, max |q - q_ref| / q_f")
    q_f = fine.sig[-1, 0] - fine.sig[-1, 2]
    for n, color, ls in zip([5, 20, 100], SERIES, LINESTYLES):
        h = plane_strain(props, 100.0, 0.10, n)
        b = (h.sig[1:, 1] - h.sig[1:, 2]) / (h.sig[1:, 0] - h.sig[1:, 2])
        b_ref = np.interp(h.eps[1:, 0], fine.eps[1:, 0], b_fine)
        q_ref = np.interp(h.eps[1:, 0], fine.eps[1:, 0], fine.sig[1:, 0] - fine.sig[1:, 2])
        summary(f"  {n:4d} {np.max(np.abs(b - b_ref)):.1e} {np.max(np.abs(h.sig[1:, 0] - h.sig[1:, 2] - q_ref)) / q_f:.1e}")
        ax.plot(100 * h.eps[1:, 0], b, color=color, ls=ls, marker="o" if n <= 20 else None, ms=4, mec="white",
                mew=0.8, label=rf"{n} steps ($\Delta\varepsilon_1$ = {0.10 / n:.0e})")
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [%]")
    ax.set_ylabel(r"$b = (\sigma_2 - \sigma_3)/(\sigma_1 - \sigma_3)$ [-]")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 0.6)
    ax.legend(loc="lower right")
    save(fig, "step_size")


# --------------------------------------------------------------------------- #
# Verification summary (Table 1 of the report)                                  #
# --------------------------------------------------------------------------- #
def verification_summary(umat):
    """Main verification checks of the report."""
    props = props_of()
    c, phi = BASE["c"], BASE["phi"]
    out = {}

    # returned stress on the surface, all Lode angles
    out["lode sweep: max |I1 I2/(I3 k) - 1|"] = max(abs(mn_criterion(s, c, phi))
                                                    for s, _ in lode_sweep(props, umat=umat))

    # one large general increment: surface and tangent
    s0 = np.array([120.0, 90.0, 100.0, 10.0, -5.0, 8.0])
    d0 = np.array([1.2e-2, -6e-3, 3e-3, 9e-3, -3e-3, 6e-3])
    s, ddsdde, _ = step(s0, d0, props, umat)
    out["general increment: |I1 I2/(I3 k) - 1|"] = abs(mn_criterion(s, c, phi))
    h = 1e-7
    fd = np.array([(step(s0, d0 + h * e, props, umat)[0] - step(s0, d0 - h * e, props, umat)[0]) / (2 * h)
                   for e in np.eye(6)]).T
    out["general increment: |D - D_FD| / |D_FD|"] = np.linalg.norm(ddsdde - fd) / np.linalg.norm(fd)

    # frame invariance: the same increment in a rotated frame
    a_, b_ = 0.4, -0.7
    R = np.array([[np.cos(a_), -np.sin(a_), 0], [np.sin(a_), np.cos(a_), 0], [0, 0, 1]]) @ \
        np.array([[1, 0, 0], [0, np.cos(b_), -np.sin(b_)], [0, np.sin(b_), np.cos(b_)]])
    rot = lambda v: voigt(R @ tensor(v) @ R.T)
    rot_strain = lambda e: voigt(R @ tensor(np.r_[e[:3], e[3:] / 2]) @ R.T) * np.r_[1, 1, 1, 2, 2, 2]
    s_rot, _, _ = step(rot(s0), rot_strain(d0), props, umat)
    out["rotated frame: |R^T s_rot R - s| / |s|"] = np.linalg.norm(voigt(R.T @ tensor(s_rot) @ R) - s) / np.linalg.norm(s)

    # tension beyond the apex
    s_t, _, state = step([50.0, 50.0, 50.0, 0, 0, 0], [-1e-2, -1e-2, -1e-2, 0, 0, 0], props, umat)
    out["isotropic tension: stress (apex: -a = %.2f)" % -attraction(c, phi)] = \
        np.array2string(s_t[:3], precision=3) + f", state {state:.0f}"

    # phi = 0
    s_0, _, _ = step([100.0, 100.0, 100.0, 0, 0, 0], [0.01, -0.005, -0.005, 0, 0, 0], props_of(c=20.0, phi=0.0, psi=0.0),
                     umat)
    out["phi = 0, c = 20: q after triaxial compression (2c = 40)"] = s_0[0] - s_0[2]

    # global Newton iterations with the tangent of the library
    res = global_newton(umat, "consistent", max_iter=40)
    out["global Newton: iterations to 1e-10 kPa"] = len(res) if res[-1] < 1e-10 else f">{len(res)}"

    # elastic strain energy of an elastic increment from zero stress
    _, _, _, sse, _ = umat(np.zeros(6), [0.0], np.zeros(6), -np.array([1e-4, 0, 0, 0, 0, 0]), props)
    sigma = fu.elastic_stiffness(BASE["E"], BASE["nu"]) @ np.array([1e-4, 0, 0, 0, 0, 0])
    out["SSE / (sigma : eps / 2), elastic from zero stress"] = sse / (0.5 * sigma @ np.array([1e-4, 0, 0, 0, 0, 0]))
    return out


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    fu.apply_style()
    print(f"UMAT: {DLL}")
    summary("\n[verification summary]")
    for key, value in verification_summary(UMAT).items():
        summary(f"  {key:58s} {value:.1e}" if isinstance(value, float) else f"  {key:58s} {value}")
    fig_yield_surface()
    fig_triaxial()
    fig_plane_strain()
    fig_undrained()
    fig_numerics()


if __name__ == "__main__":
    main()
