"""
Figures for the Hardening Soil test report (docs/hs_test_report/hs_test_report.tex).

Drives the Hardening Soil UMAT (build_C/lib/hardening_soil.dll or .so) along element test paths:
the analytic checks of tests/test_hardening_soil.py and the Hostun sand calibration and
verification of Schanz, Vermeer & Bonnier (1999), "The hardening soil model: Formulation and
verification". Equation and figure numbers refer to that paper.

    python docs/hs_test_report/make_figures.py

Stresses and strains are compression positive, as in the paper and in the tests; the UMAT interface
is tension positive and `step` converts. The figures are written to docs/hs_test_report/figures and
the values quoted in the report are printed.
"""

import logging
import os
import sys

import cffi
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIG_DIR = os.path.join(HERE, "figures")
DLL = os.path.join(ROOT, "build_C", "lib", "hardening_soil." + ("dll" if sys.platform == "win32" else "so"))

# --------------------------------------------------------------------------- #
# UMAT                                                                          #
# --------------------------------------------------------------------------- #
_ffi = cffi.FFI()
_ffi.cdef("""
void umat(double* STRESS, double* STATEV, double* DDSDDE,
          double* SSE, double* SPD, double* SCD, double* RPL,
          double* DDSDDT, double* DRPLDE, double* DRPLDT,
          double* STRAN, double* DSTRAN, double* TIME, double* DTIME,
          double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
          char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV,
          double* PROPS, int* NPROPS, double* COORDS, double* DROT,
          double* PNEWDT, double* CELENT, double* DFGRD0, double* DFGRD1,
          int* NOEL, int* NPT, int* LAYER, int* KSPT, int* KSTEP, int* KINC);
""")
_lib = None


def _umat(stress, statev, strain, dstrain, props):
    """One UMAT call in the tension-positive convention of the interface (same arguments as
    tests/utils.py, Utils.run_c_umat, but with the library loaded once)."""
    global _lib
    if _lib is None:
        if not os.path.exists(DLL):
            sys.exit(f"{DLL} not found: build the C models first (see README.md)")
        _lib = _ffi.dlopen(DLL)
    ints = [_ffi.new("int*", v) for v in (3, 3, 6, len(statev), len(props), 1, 1, 1, 1)]
    c_stress = _ffi.new("double[]", list(stress))
    c_statev = _ffi.new("double[]", list(statev))
    c_ddsdde = _ffi.new("double[]", 36)
    scalars = [_ffi.new("double*", 0.0) for _ in range(3)]
    _lib.umat(c_stress, c_statev, c_ddsdde, *scalars, _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL,
              _ffi.new("double[]", list(strain)), _ffi.new("double[]", list(dstrain)),
              _ffi.new("double[]", [0.0, 0.0]), _ffi.new("double*", 1.0),
              _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.new("char[]", b"HS".ljust(80)),
              ints[0], ints[1], ints[2], ints[3], _ffi.new("double[]", list(props)), ints[4],
              _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL,
              ints[5], ints[6], _ffi.NULL, _ffi.NULL, ints[7], ints[8])
    return np.array(list(c_stress)), np.array(list(c_statev))


def step(stress, statev, strain, dstrain, props):
    """UMAT call in the compression-positive convention of the paper."""
    s, sv = _umat(-np.asarray(stress, float), statev, -np.asarray(strain, float),
                  -np.asarray(dstrain, float), props)
    return -s, sv


# --------------------------------------------------------------------------- #
# Parameters                                                                    #
# --------------------------------------------------------------------------- #
NAMES = ["E50_ref", "Eur_ref", "m", "c", "phi", "psi", "p_ref", "Rf", "nu", "M_cap", "K_ratio", "e0", "e_cv"]

# unit-test parameter set (tests/test_hardening_soil.py)
TEST = dict(E50_ref=30000.0, Eur_ref=90000.0, m=0.55, c=0.0, phi=42.0, psi=16.0, p_ref=100.0, Rf=0.9,
            nu=0.25, M_cap=1.0, K_ratio=1.84)

# loose Hostun sand, Table 1 of the paper (Rf = 0.9 is the default of Sec. 2, K0_NC = 1 - sin(phi) of
# Sec. 5.2); M_cap and K_ratio are calibrated below to K0_NC and Eoed_ref = 0.8 E50_ref
HOSTUN = dict(E50_ref=20000.0, Eur_ref=60000.0, m=0.65, c=0.0, phi=34.0, psi=0.0, p_ref=100.0, Rf=0.9,
              nu=0.2, M_cap=None, K_ratio=None)
HOSTUN_EOED_REF = 0.8 * HOSTUN["E50_ref"]
HOSTUN_K0_NC = 1.0 - np.sin(np.radians(HOSTUN["phi"]))
# dilatant variant (Table 1 gives psi = 0)
HOSTUN_PSI_VARIANT = 14


def props_of(base, **overrides):
    values = dict(base, **overrides)
    return [values[n] for n in NAMES if n in values]


def k_f(phi_deg):
    """qf = k_f (sigma_3 + a), Eq. 2 written in sigma_3."""
    s = np.sin(np.radians(phi_deg))
    return 2.0 * s / (1.0 - s)


def sin_psi_mobilised(sin_phi_m, sin_phi, sin_psi):
    """Flow rule of the UMAT: Rowe (Eq. 11) bounded so that the shear surface does not contract."""
    sin_phi_m = np.asarray(sin_phi_m, float)
    sin_cv = (sin_phi - sin_psi) / (1.0 - sin_phi * sin_psi)
    s = np.minimum(sin_phi_m, sin_phi)
    above = sin_psi if sin_psi <= 0.0 else np.maximum((s - sin_cv) / (1.0 - s * sin_cv), 0.0)
    return np.where(sin_phi_m < 0.75 * sin_phi, 0.0, above)


# --------------------------------------------------------------------------- #
# Element test paths                                                            #
# --------------------------------------------------------------------------- #
class History:
    def __init__(self, strain, stress, statev):
        self.eps, self.sig, self.sv = [strain.copy()], [stress.copy()], [statev.copy()]

    def add(self, strain, stress, statev):
        self.eps.append(strain.copy())
        self.sig.append(stress.copy())
        self.sv.append(statev.copy())

    def arrays(self):
        self.eps, self.sig, self.sv = np.array(self.eps), np.array(self.sig), np.array(self.sv)
        return self


def mixed_step(stress, statev, strain, dstrain, groups, targets, props, tol=1e-9):
    """Strain increment in which each group of components (sharing one strain increment) is
    adjusted by Newton iteration so that the stress of its first component equals the target
    (as in tests/test_hardening_soil.py)."""
    d = np.array(dstrain, dtype=float)
    u = np.array([d[g[0]] for g in groups])
    first = [g[0] for g in groups]

    def apply(values):
        out = d.copy()
        for g, v in zip(groups, values):
            out[list(g)] = v
        return out

    for _ in range(50):
        s, _ = step(stress, statev, strain, apply(u), props)
        r = s[first] - targets
        if np.max(np.abs(r)) < tol * (1.0 + np.max(np.abs(targets))):
            break
        jac = np.zeros((len(groups), len(groups)))
        h = 1e-8
        for a in range(len(groups)):
            u_h = u.copy()
            u_h[a] += h
            s_h, _ = step(stress, statev, strain, apply(u_h), props)
            jac[:, a] = (s_h[first] - s[first]) / h
        u -= np.linalg.solve(jac, r)

    d = apply(u)
    s, sv = step(stress, statev, strain, d, props)
    return s, sv, strain + d


def drained_triaxial(props, cell, axial_strain, n_steps, unload_steps=0, d_unload=2e-4):
    """Drained triaxial compression at constant cell pressure (axial = x), optionally followed by
    axial unloading until q = 0."""
    stress, statev, strain = np.full(6, cell) * [1, 1, 1, 0, 0, 0], np.zeros(3), np.zeros(6)
    h = History(strain, stress, statev)
    for _ in range(n_steps):
        stress, statev, strain = mixed_step(stress, statev, strain, [axial_strain / n_steps, 0, 0, 0, 0, 0],
                                            [[1, 2]], np.array([cell]), props)
        h.add(strain, stress, statev)
    for _ in range(unload_steps):
        stress, statev, strain = mixed_step(stress, statev, strain, [-d_unload, 0, 0, 0, 0, 0],
                                            [[1, 2]], np.array([cell]), props)
        h.add(strain, stress, statev)
        if stress[0] <= cell:
            break
    return h.arrays()


def strain_path(props, stress0, increments, stop=None):
    """Strain-controlled path; `stop(stress)` ends a segment early."""
    stress, statev, strain = np.array(stress0, float), np.zeros(3), np.zeros(6)
    h = History(strain, stress, statev)
    for d in increments:
        stress, statev = step(stress, statev, strain, d, props)
        strain = strain + d
        h.add(strain, stress, statev)
        if stop is not None and stop(stress):
            break
    return h.arrays()


def stress_targets(props, stress0, direction, targets, d_eps, component=0):
    """Strain-controlled loading along `direction`, reversing at each target value of the stress
    component (load / unload / reload sequences of the oedometer and isotropic tests)."""
    stress, statev, strain = np.array(stress0, float), np.zeros(3), np.zeros(6)
    h = History(strain, stress, statev)
    sign = 1.0
    for target in targets:
        sign = 1.0 if target > stress[component] else -1.0
        for _ in range(100000):
            d = sign * d_eps * np.asarray(direction, float)
            stress, statev = step(stress, statev, strain, d, props)
            strain = strain + d
            h.add(strain, stress, statev)
            if sign * (stress[component] - target) >= 0.0:
                break
    return h.arrays()


def undrained_triaxial(props, cell, axial_strain, n_steps):
    """Undrained (isochoric) triaxial compression from an isotropic state; the total lateral
    stress stays equal to the cell pressure, so the excess pore pressure is cell - sigma_3'."""
    d = axial_strain / n_steps
    return strain_path(props, [cell, cell, cell, 0, 0, 0], [np.array([d, -d / 2, -d / 2, 0, 0, 0])] * n_steps)


def oedometer(props, sigma_v0, k0, targets, d_eps=2e-5):
    return stress_targets(props, [sigma_v0, k0 * sigma_v0, k0 * sigma_v0, 0, 0, 0], [1, 0, 0, 0, 0, 0],
                          targets, d_eps)


def p_q(sig):
    p = sig[:, :3].mean(axis=1)
    q = np.sqrt(0.5 * ((sig[:, 0] - sig[:, 1]) ** 2 + (sig[:, 1] - sig[:, 2]) ** 2
                       + (sig[:, 2] - sig[:, 0]) ** 2) + 3 * (sig[:, 3:] ** 2).sum(axis=1))
    return p, q


# --------------------------------------------------------------------------- #
# Calibration of the cap parameters (Sec. 4-5)                                  #
# --------------------------------------------------------------------------- #
def oedometer_response(props, sigma_v=None):
    """Tangent K0 = d sigma_h / d sigma_v and tangent Eoed = d sigma_v / d eps_v of normally
    consolidated oedometer loading at sigma_v (default p_ref)."""
    p_ref = props[NAMES.index("p_ref")]
    sigma_v = p_ref if sigma_v is None else sigma_v
    h = oedometer(props, 0.2 * sigma_v, HOSTUN_K0_NC, [1.6 * sigma_v], d_eps=1e-5)
    sv, sh, e1 = h.sig[:, 0], h.sig[:, 1], h.eps[:, 0]
    mid = 0.5 * (sv[1:] + sv[:-1])
    k0 = np.interp(sigma_v, mid, np.diff(sh) / np.diff(sv))
    eoed = np.interp(sigma_v, mid, np.diff(sv) / np.diff(e1))
    return k0, eoed


def calibrate_cap(base, k0_target, eoed_target, x0=(1.5, 2.0)):
    """Newton iteration on (M_cap, K_ratio) such that NC oedometer loading has K0 = k0_target and
    Eoed(sigma_v = p_ref) = eoed_target, from the input parameters K0_NC and Eoed_ref (Sec. 4)."""
    x = np.array(x0, float)

    def residual(x):
        k0, eoed = oedometer_response(props_of(base, M_cap=x[0], K_ratio=x[1]))
        return np.array([k0 - k0_target, eoed / eoed_target - 1.0])

    for _ in range(30):
        r = residual(x)
        if np.max(np.abs(r)) < 1e-6:
            break
        jac = np.zeros((2, 2))
        for j in range(2):
            dx = np.zeros(2)
            dx[j] = 1e-4 * x[j]
            jac[:, j] = (residual(x + dx) - r) / dx[j]
        dx = -np.linalg.solve(jac, r)
        t = 1.0
        while np.any(x + t * dx <= [0.0, 1.0]):  # M_cap > 0, K_ratio > 1
            t *= 0.5
        x = x + t * dx
    return x, residual(x)


# --------------------------------------------------------------------------- #
# Plot style (static PDF figures; palette validated for colour-vision         #
# deficiency, see the report)                                                   #
# --------------------------------------------------------------------------- #
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
REF = dict(color=INK, lw=0.9, ls=(0, (4, 2)))       # analytic solutions
ENVELOPE = dict(color=INK2, lw=0.9, ls=(0, (1, 1.5)))  # failure envelopes and limits

plt.rcParams.update({
    "figure.figsize": (3.9, 2.9), "font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 8.5,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.5, "lines.solid_capstyle": "round", "legend.frameon": False,
    "legend.fontsize": 7.5, "legend.handlelength": 2.4, "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})
logging.getLogger("fontTools").setLevel(logging.ERROR)  # font timestamp warnings when embedding


def log_axes(ax, xticks, yticks):
    """Log-log axes with plain tick labels."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    for axis, ticks in [(ax.xaxis, xticks), (ax.yaxis, yticks)]:
        axis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())


def save(fig, name):
    fig.savefig(os.path.join(FIG_DIR, name + ".pdf"))
    plt.close(fig)


def summary(text):
    print(text)


# --------------------------------------------------------------------------- #
# Figures                                                                       #
# --------------------------------------------------------------------------- #
def fig_hyperbola():
    """Drained triaxial without dilatancy (cut-off active) and without cap: Eq. 1 and Eqs. 3-4."""
    props = props_of(TEST, M_cap=0.0, e0=0.8, e_cv=0.5)
    ei_ref = 2.0 * TEST["E50_ref"] / (2.0 - TEST["Rf"])

    fig, ax = plt.subplots()
    summary("\n[hyperbola] cell, q_end/qf, max rel. deviation from Eq. 1")
    for cell, color in zip([100.0, 200.0, 300.0], SERIES):
        h = drained_triaxial(props, cell, 0.05, 100)
        e1, q = h.eps[:, 0], h.sig[:, 0] - h.sig[:, 2]
        factor = (cell / TEST["p_ref"]) ** TEST["m"]
        qf = k_f(TEST["phi"]) * cell
        q_line = np.linspace(0.0, q[-1], 200)
        e_line = q_line / (ei_ref * factor * (1.0 - q_line * TEST["Rf"] / qf))
        dev = np.max(np.abs(e1[1:] - q[1:] / (ei_ref * factor * (1.0 - q[1:] * TEST["Rf"] / qf))) / e1[1:])
        summary(f"  {cell:5.0f} {q[-1] / qf:.4f} {dev:.2e}")
        ax.plot(100 * e1, q, color=color, label=rf"UMAT, $\sigma_3$ = {cell:.0f} kPa")
        ax.plot(100 * e_line, q_line, **REF)
        ax.hlines(qf, 3.0, 5.0, **ENVELOPE)  # kept clear of the legend
    ax.plot([], [], **REF, label=r"hyperbola, $E_i = 2E_{50}/(2-R_f)$")
    ax.plot([], [], **ENVELOPE, label=r"$q_f$ (Eq. 2)")
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [%]")
    ax.set_ylabel(r"deviator stress $q = \sigma_1 - \sigma_3$ [kPa]")
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 1500)
    ax.legend(loc="upper left")
    save(fig, "triaxial_hyperbola")

    # secant stiffness at q = qf / 2 as a function of the confining pressure (Eq. 3)
    cells = np.array([25.0, 50.0, 100.0, 200.0, 400.0, 800.0])
    e50 = []
    for cell in cells:
        qf = k_f(TEST["phi"]) * cell
        e_half = 0.5 * qf / (TEST["E50_ref"] * (cell / TEST["p_ref"]) ** TEST["m"])
        h = drained_triaxial(props, cell, 2.5 * e_half, 100)
        q = h.sig[:, 0] - h.sig[:, 2]
        e50.append(0.5 * qf / np.interp(0.5 * qf, q, h.eps[:, 0]))
    e50 = np.array(e50)
    eq3 = TEST["E50_ref"] * (cells / TEST["p_ref"]) ** TEST["m"]
    summary("[E50] sigma3, E50 UMAT, E50 Eq. 3, rel. diff")
    for c, a, b in zip(cells, e50, eq3):
        summary(f"  {c:5.0f} {a:10.1f} {b:10.1f} {a / b - 1:+.2e}")

    fig, ax = plt.subplots()
    s_line = np.geomspace(20, 1000, 50)
    ax.plot(s_line, TEST["E50_ref"] * (s_line / TEST["p_ref"]) ** TEST["m"] / 1000, **REF,
            label=rf"Eq. 3, $m$ = {TEST['m']}")
    ax.plot(cells, e50 / 1000, "o", color=SERIES[0], ms=6, mec="white", mew=1.0,
            label=r"UMAT, secant at $q_f/2$")
    log_axes(ax, [25, 50, 100, 200, 400, 800], [10, 20, 50, 100])
    ax.set_xlabel(r"confining pressure $\sigma_3$ [kPa]")
    ax.set_ylabel(r"$E_{50}$ [MPa]")
    ax.legend(loc="upper left")
    save(fig, "triaxial_e50")


def fig_dilatancy():
    """Rowe stress dilatancy (Eqs. 9-15) and the dilatancy cut-off (Eqs. 38-39), with m = 0 so the
    plastic strains follow from the total strains with a constant elastic stiffness."""
    base = props_of(TEST, m=0.0, M_cap=0.0)
    e0, e_cv = 0.60, 0.65
    h = drained_triaxial(base, 100.0, 0.12, 120)
    hc = drained_triaxial(props_of(TEST, m=0.0, M_cap=0.0, e0=e0, e_cv=e_cv), 100.0, 0.12, 120)

    fig, ax = plt.subplots()
    ev = 100 * h.eps[:, :3].sum(axis=1)
    evc = 100 * hc.eps[:, :3].sum(axis=1)
    ev_cv = -100 * np.log((1 + e_cv) / (1 + e0))  # Eq. 39
    ax.plot(100 * h.eps[:, 0], ev, color=SERIES[0], label="cut-off off")
    ax.plot(100 * hc.eps[:, 0], evc, color=SERIES[1], ls=(0, (6, 2)),
            label=rf"cut-off on ($e_0$ = {e0}, $e_{{cv}}$ = {e_cv})")
    ax.axhline(ev_cv, **ENVELOPE, label=r"$e = e_{cv}$ (Eq. 39)")
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [%]")
    ax.set_ylabel(r"volumetric strain $\varepsilon_v$ [%] (compression +)")
    ax.legend(loc="lower left")
    save(fig, "dilatancy_ev")
    summary(f"\n[dilatancy] eps_v at e_cv: {ev_cv:.3f} %, final eps_v with cut-off: {evc[-1]:.3f} %,"
            f" without: {ev[-1]:.3f} %")

    # plastic strains, mobilised friction and dilatancy angles
    E, nu = TEST["Eur_ref"], TEST["nu"]
    G, lam = E / (2 * (1 + nu)), E * nu / ((1 + nu) * (1 - 2 * nu))
    De = np.zeros((6, 6))
    De[:3, :3] = lam
    De[np.arange(3), np.arange(3)] += 2 * G
    De[3:, 3:] = G * np.eye(3)
    eps_p = h.eps - np.linalg.solve(De, (h.sig - h.sig[0]).T).T
    dev_p, dgamma = np.diff(eps_p[:, :3].sum(axis=1)), np.diff(h.sv[:, 0])
    s1, s3 = h.sig[:, 0], h.sig[:, 2]
    sin_phi_m = (s1 - s3) / (s1 + s3)
    sin_phi_m_mid = 0.5 * (sin_phi_m[1:] + sin_phi_m[:-1])
    psi_m_umat = np.degrees(np.arcsin(-dev_p / dgamma))

    sin_phi, sin_psi = np.sin(np.radians(TEST["phi"])), np.sin(np.radians(TEST["psi"]))
    sin_cv = (sin_phi - sin_psi) / (1 - sin_phi * sin_psi)
    sp = np.linspace(0.0, sin_phi, 400)
    psi_rowe = np.degrees(np.arcsin((sp - sin_cv) / (1 - sp * sin_cv)))
    psi_bounded = np.degrees(np.arcsin(sin_psi_mobilised(sp, sin_phi, sin_psi)))

    gamma_check = np.max(np.abs(h.sv[1:, 0] - (eps_p[1:, 0] - eps_p[1:, 1] - eps_p[1:, 2])) / h.sv[1:, 0])
    # psi_m is evaluated at the start of each sub-step, so the ratio of a load step lies between the
    # values of the flow rule at its start and at its end
    ratio_umat = -dev_p / dgamma
    ends = np.sort([sin_psi_mobilised(sin_phi_m[:-1], sin_phi, sin_psi),
                    sin_psi_mobilised(sin_phi_m[1:], sin_phi, sin_psi)], axis=0)
    outside = np.max(np.maximum(0.0, np.maximum(ends[0] - ratio_umat, ratio_umat - ends[1])))
    below_cv = ends[1] == 0.0
    sin_phi_start = sin_phi_m[1:][~below_cv][0]
    summary(f"[rowe] phi_cv = {np.degrees(np.arcsin(sin_cv)):.2f} deg, 3/4 sin(phi): phi_m ="
            f" {np.degrees(np.arcsin(0.75 * sin_phi)):.2f} deg, max rel. error Eq. 9: {gamma_check:.1e},"
            f" max |sin(psi_m)| below the bound: {np.max(np.abs(ratio_umat[below_cv])):.1e} ({below_cv.sum()} steps),"
            f" max distance from the rule: {outside:.1e}, dilatant from phi_m ="
            f" {np.degrees(np.arcsin(sin_phi_start)):.2f} deg, psi_m at failure: {psi_m_umat[-1]:.3f} deg")

    fig, ax = plt.subplots()
    ax.plot(np.degrees(np.arcsin(sp)), psi_rowe, **ENVELOPE, label="Eq. 11 (Rowe)")
    ax.plot(np.degrees(np.arcsin(sp)), psi_bounded, **REF, label="Eq. 11 bounded")
    ax.plot(np.degrees(np.arcsin(sin_phi_m_mid)), psi_m_umat, "o", color=SERIES[0], ms=4.5, mec="white",
            mew=0.8, label=r"UMAT, $-\Delta\varepsilon_v^p / \Delta\gamma^p$")
    ax.axvline(np.degrees(np.arcsin(sin_cv)), color=AXIS, lw=0.8)
    ax.axhline(0.0, color=AXIS, lw=0.8)
    ax.annotate(r"$\varphi_{cv}$", (np.degrees(np.arcsin(sin_cv)), -27), xytext=(4, 0),
                textcoords="offset points", color=INK2, fontsize=7.5)
    ax.set_xlabel(r"mobilised friction angle $\varphi_m$ [$^\circ$]")
    ax.set_ylabel(r"mobilised dilatancy angle $\psi_m$ [$^\circ$]")
    ax.legend(loc="upper left")
    save(fig, "dilatancy_rowe")


def fig_cap():
    """Isotropic compression with an unloading-reloading loop (Eqs. 31-36)."""
    props = props_of(TEST)
    m, p_ref = TEST["m"], TEST["p_ref"]
    ks_ref = TEST["Eur_ref"] / (3 * (1 - 2 * TEST["nu"]))
    kc_ref = ks_ref / TEST["K_ratio"]
    h = stress_targets(props, [100, 100, 100, 0, 0, 0], [1, 1, 1, 0, 0, 0], [400, 150, 800], 2e-5)
    p, ev = h.sig[:, 0], h.eps[:, :3].sum(axis=1)

    def branch(p0, e0, p1, k_ref, n=200):
        """p(eps_v) for dp/deps_v = k_ref (p / p_ref)^m from (e0, p0) to p1."""
        pp = np.linspace(p0, p1, n)
        return e0 + (pp ** (1 - m) - p0 ** (1 - m)) / ((1 - m) * k_ref * p_ref ** (-m)), pp

    # analytic solution through the reversal points of the UMAT run (the last increment of each
    # branch overshoots the target stress slightly)
    i_top = np.argmax(p >= 400)
    i_bottom = i_top + np.argmax(p[i_top:] <= 150)
    p_top, p_bottom = p[i_top], p[i_bottom]
    e_exact, p_exact, ends = [], [], []
    e_a, p_a = 0.0, 100.0
    for p_b, k in [(p_top, kc_ref), (p_bottom, ks_ref), (p_top, ks_ref), (p[-1], kc_ref)]:
        e_line, p_line = branch(p_a, e_a, p_b, k)
        e_exact.append(e_line)
        p_exact.append(p_line)
        e_a, p_a = e_line[-1], p_line[-1]
        ends.append(e_a)
    e_exact, p_exact = np.concatenate(e_exact), np.concatenate(p_exact)
    summary("\n[cap] p, eps_v UMAT, eps_v analytic, rel. diff")
    for i, e_ref in zip([i_top, i_bottom, len(p) - 1], [ends[0], ends[1], ends[3]]):
        summary(f"  {p[i]:7.2f} {ev[i]:.6f} {e_ref:.6f} {ev[i] / e_ref - 1:+.2e}")
    summary(f"  p_c at the end: {h.sv[-1, 1]:.2f}, gamma_p: {h.sv[-1, 0]:.1e}")

    fig, ax = plt.subplots()
    ax.plot(100 * ev, p, color=SERIES[0], label="UMAT")
    ax.plot(100 * e_exact, p_exact, **REF, label="Eqs. 31-36, $K_s/K_c$ = %.2f" % TEST["K_ratio"])
    ax.set_xlabel(r"volumetric strain $\varepsilon_v$ [%]")
    ax.set_ylabel(r"mean stress $p$ [kPa]")
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.legend(loc="upper left")
    save(fig, "cap_isotropic")

    # tangent bulk modulus; primary loading when p_c grows from a state on the cap, elastic when p_c
    # is unchanged (the increment in which reloading reaches the cap is partly elastic and is left out)
    k_t = np.diff(p) / np.diff(ev)
    p_mid = 0.5 * (p[1:] + p[:-1])
    p_c = h.sv[:, 1]
    virgin = (np.diff(p_c) > 0) & (p[:-1] >= p_c[:-1] * (1 - 1e-9))
    elastic = np.diff(p_c) == 0
    fig, ax = plt.subplots()
    pl = np.geomspace(80, 1000, 50)
    ax.plot(pl, ks_ref * (pl / p_ref) ** m / 1000, **REF, label=r"$K_s = E_{ur}/3(1-2\nu)$")
    ax.plot(pl, kc_ref * (pl / p_ref) ** m / 1000, **ENVELOPE, label=r"$K_c = K_s / (K_s/K_c)$")
    sel = np.arange(len(k_t))
    for mask, color, label, marker in [(virgin, SERIES[0], "UMAT, primary loading", "o"),
                                       (elastic, SERIES[1], "UMAT, unloading / reloading", "s")]:
        idx = sel[mask][::25]
        ax.plot(p_mid[idx], k_t[idx] / 1000, marker, color=color, ms=4.5, mec="white", mew=0.8, ls="none",
                label=label)
    log_axes(ax, [100, 200, 400, 800], [20, 50, 100, 200, 500])
    ax.set_ylim(25, 600)
    ax.set_xlabel(r"mean stress $p$ [kPa]")
    ax.set_ylabel(r"tangent bulk modulus $dp/d\varepsilon_v$ [MPa]")
    ax.legend(loc="upper left")
    save(fig, "cap_bulk_modulus")
    summary("[cap] max rel. deviation of the tangent modulus: primary %.2e, un/reloading %.2e"
            " (%d of %d increments are transitions)" % (
                np.max(np.abs(k_t[virgin] / (kc_ref * (p_mid[virgin] / p_ref) ** m) - 1)),
                np.max(np.abs(k_t[elastic] / (ks_ref * (p_mid[elastic] / p_ref) ** m) - 1)),
                np.sum(~virgin & ~elastic), len(k_t)))


def fig_hostun(m_cap, k_ratio):
    hostun = dict(HOSTUN, M_cap=m_cap, K_ratio=k_ratio)
    props = props_of(hostun)
    phi = hostun["phi"]
    sin_phi = np.sin(np.radians(phi))

    # --- oedometer test with an unloading-reloading loop (Fig. 5) ---
    h = oedometer(props, 2.0, HOSTUN_K0_NC, [200.0, 10.0, 300.0], d_eps=1e-5)
    sv, sh, e1 = h.sig[:, 0], h.sig[:, 1], h.eps[:, 0]
    fig, ax = plt.subplots()
    ax.plot(100 * e1, sv, color=SERIES[0], label="UMAT")
    sl = np.linspace(2.0, 300.0, 200)
    m, p_ref = hostun["m"], hostun["p_ref"]
    e_eq37 = (sl ** (1 - m) - 2.0 ** (1 - m)) / ((1 - m) * HOSTUN_EOED_REF * p_ref ** (-m))
    ax.plot(100 * e_eq37, sl, **REF, label=r"Eq. 37, $E_{oed}^{ref}$ = %.0f MPa" % (HOSTUN_EOED_REF / 1000))
    ax.set_xlabel(r"vertical strain $\varepsilon_1$ [%]")
    ax.set_ylabel(r"vertical stress $\sigma_1$ [kPa]")
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.legend(loc="upper left")
    save(fig, "hostun_oedometer")

    eoed = np.diff(sv) / np.diff(e1)
    mid = 0.5 * (sv[1:] + sv[:-1])
    virgin = (np.diff(sv) > 0) & (np.maximum.accumulate(sv[:-1]) <= sv[:-1] + 1e-9) & (mid > 20)
    k0_t = np.diff(sh) / np.diff(sv)
    summary("\n[hostun oedometer] sigma_v, Eoed UMAT, Eoed Eq. 37, K0 tangent, K0 secant")
    for s in [50.0, 100.0, 200.0, 280.0]:
        i = np.argmin(np.abs(np.where(virgin, mid, np.inf) - s))
        summary(f"  {mid[i]:6.1f} {eoed[i]:9.0f} {HOSTUN_EOED_REF * (mid[i] / p_ref) ** m:9.0f}"
                f" {k0_t[i]:.4f} {sh[i] / sv[i]:.4f}")
    i_top = np.argmax(sv >= 200.0)
    i_bottom = i_top + np.argmin(sv[i_top:])
    summary(f"  eps_1 at 200 kPa: {e1[i_top]:.5f}, after unloading to {sv[i_bottom]:.1f} kPa: {e1[i_bottom]:.5f}"
            f" (Eq. 37 at 200 kPa: {np.interp(200.0, sl, e_eq37):.5f}); K0 at end of unloading:"
            f" {sh[i_bottom] / sv[i_bottom]:.3f}")

    fig, ax = plt.subplots()
    ax.plot(sv, sh, color=SERIES[0], label="UMAT")
    ax.plot([0, 300], [0, 300 * HOSTUN_K0_NC], **REF, label=r"$K_0^{NC} = 1 - \sin\varphi$ = %.3f" % HOSTUN_K0_NC)
    ax.set_xlabel(r"vertical stress $\sigma_1$ [kPa]")
    ax.set_ylabel(r"horizontal stress $\sigma_3$ [kPa]")
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.legend(loc="upper left")
    save(fig, "hostun_oedometer_path")

    # --- drained triaxial tests, sigma_3 = 300 kPa (Fig. 6) and 600 kPa ---
    # psi = 0 of Table 1 and the dilatant variant
    props_psi = props_of(hostun, psi=HOSTUN_PSI_VARIANT)
    psi_label = rf"$\psi = {HOSTUN_PSI_VARIANT:g}^\circ$"
    variant_ls = (0, (5, 2))
    ratio_f = (1 + sin_phi) / (1 - sin_phi)
    cells = [300.0, 600.0]
    drained = {(psi, cell): drained_triaxial(pv, cell, 0.18, 360, unload_steps=200, d_unload=1e-4)
               for psi, pv in [(0.0, props), (HOSTUN_PSI_VARIANT, props_psi)] for cell in cells}
    e50 = hostun["E50_ref"] * (300.0 / p_ref) ** m
    qf = k_f(phi) * 300.0
    summary("\n[hostun drained] psi, sigma_3, max sigma1/sigma3 (MC: %.4f), eps_1 at 0.999 of max, max eps_v"
            " at eps_1, eps_v at end of loading" % ratio_f)
    for (psi, cell), hv in drained.items():
        r, ev = hv.sig[:, 0] / hv.sig[:, 2], hv.eps[:, :3].sum(axis=1)
        summary(f"  {psi:4.1f} {cell:.0f} {r.max():.4f} {hv.eps[np.argmax(r >= 0.999 * r.max()), 0]:.4f}"
                f" {ev.max():.5f} {hv.eps[np.argmax(ev), 0]:.4f} {ev[360]:.5f}")

    # contribution of the cap to the secant stiffness (with psi = 0 the shear surface has no plastic
    # volume change below failure)
    summary(f"[hostun E50] variant, E50 secant at qf/2 (Eq. 3: {e50:.0f}), eps_v at qf/2")
    for name, variant in [("no cap", props_of(hostun, M_cap=0.0)), ("cap (full model)", props),
                          ("cap, psi variant", props_psi)]:
        hv = drained_triaxial(variant, 300.0, 0.03, 300)
        qv = hv.sig[:, 0] - hv.sig[:, 2]
        summary(f"  {name:22s} {0.5 * qf / np.interp(0.5 * qf, qv, hv.eps[:, 0]):8.0f}"
                f" {np.interp(0.5 * qf, qv, hv.eps[:, :3].sum(axis=1)):.5f}")

    def drained_legend(ax, loc, **kwargs):
        for cell, color in zip(cells, SERIES):
            ax.plot([], [], color=color, label=rf"$\sigma_3$ = {cell:.0f} kPa")
        ax.plot([], [], color=INK2, label=r"$\psi = 0$ (Table 1)")
        ax.plot([], [], color=INK2, ls=variant_ls, label=psi_label)
        ax.legend(loc=loc, **kwargs)

    def drained_lines(ax, y_of):
        for (psi, cell), hv in drained.items():
            ax.plot(hv.eps[:, 0], y_of(hv), color=SERIES[cells.index(cell)], ls="-" if psi == 0.0 else variant_ls)

    fig, ax = plt.subplots()
    drained_lines(ax, lambda hv: hv.sig[:, 0] / hv.sig[:, 2])
    ax.axhline(ratio_f, **ENVELOPE, label=r"failure, $(1+\sin\varphi)/(1-\sin\varphi)$")
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [-]")
    ax.set_ylabel(r"stress ratio $\sigma_1/\sigma_3$ [-]")
    ax.set_xlim(0, 0.2)
    ax.set_ylim(1, 4)
    drained_legend(ax, "lower left", bbox_to_anchor=(0.12, 0.02))  # clear of the unloading branch
    save(fig, "hostun_drained_ratio")

    fig, ax = plt.subplots()
    drained_lines(ax, lambda hv: hv.eps[:, :3].sum(axis=1))
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [-]")
    ax.set_ylabel(r"volumetric strain $\varepsilon_v$ [-] (compression +)")
    ax.set_xlim(0, 0.2)
    ax.axhline(0.0, color=AXIS, lw=0.8)
    drained_legend(ax, "lower left")
    save(fig, "hostun_drained_ev")

    # --- undrained triaxial tests, sigma_c = 300 and 600 kPa (Figs. 7-8) ---
    runs = {(psi, cell): undrained_triaxial(pv, cell, 0.20, 800)
            for psi, pv in [(0.0, props), (HOSTUN_PSI_VARIANT, props_psi)] for cell in cells}

    def undrained_legend(ax, loc, ncol=1):
        for cell, color in zip(cells, SERIES):
            ax.plot([], [], color=color, label=rf"$\sigma_c$ = {cell:.0f} kPa")
        ax.plot([], [], color=INK2, label=r"$\psi = 0$ (Table 1)")
        ax.plot([], [], color=INK2, ls=variant_ls, label=psi_label)
        ax.legend(loc=loc, ncol=ncol, columnspacing=1.0)

    def undrained_lines(ax, x_of, y_of):
        for (psi, cell), hv in runs.items():
            ax.plot(x_of(hv, cell), y_of(hv, cell), color=SERIES[cells.index(cell)],
                    ls="-" if psi == 0.0 else variant_ls)

    fig, ax = plt.subplots()
    undrained_lines(ax, lambda hv, cell: p_q(hv.sig)[0], lambda hv, cell: p_q(hv.sig)[1])
    m_mc = 6 * sin_phi / (3 - sin_phi)
    pl = np.array([10.0, 2e4 / m_mc])
    ax.plot(pl, m_mc * pl, **ENVELOPE, label=rf"MC, $q = {m_mc:.2f}\,p'$")
    ax.set_xlabel(r"mean effective stress $p'$ [kPa]")
    ax.set_ylabel(r"deviator stress $q$ [kPa]")
    log_axes(ax, [100, 200, 500, 1000, 2000, 5000, 10000], [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000])
    ax.set_xlim(100, 10000)
    ax.set_ylim(10, 10000)
    undrained_legend(ax, "lower right")
    save(fig, "hostun_undrained_pq")

    fig, ax = plt.subplots()
    undrained_lines(ax, lambda hv, cell: hv.eps[:, 0], lambda hv, cell: hv.sig[:, 0] / hv.sig[:, 2])
    ax.axhline(ratio_f, **ENVELOPE, label="failure")
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [-]")
    ax.set_ylabel(r"stress ratio $\sigma_1'/\sigma_3'$ [-]")
    ax.set_xlim(0, 0.2)
    ax.set_ylim(1, 4)
    undrained_legend(ax, "lower right")
    save(fig, "hostun_undrained_ratio")

    fig, ax = plt.subplots()
    undrained_lines(ax, lambda hv, cell: hv.eps[:, 0], lambda hv, cell: (cell - hv.sig[:, 2]) / cell)
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [-]")
    ax.set_ylabel(r"excess pore pressure $\Delta u / \sigma_c$ [-]")
    ax.set_xlim(0, 0.2)
    ax.axhline(0.0, color=AXIS, lw=0.8)
    undrained_legend(ax, "lower left", ncol=2)
    save(fig, "hostun_undrained_du")

    summary("[hostun undrained] psi, cell | on failure line: eps_1, p', q, du/sc | du/sc max at eps_1 |"
            " eps_1 = 0.14: p', du/sc | end: p', q, du/sc, phi'_mob")
    for (psi, cell), hv in runs.items():
        p, qq = p_q(hv.sig)
        du = (cell - hv.sig[:, 2]) / cell
        e1 = hv.eps[:, 0]
        i_f = np.argmax(hv.sv[:, 2] > 0.5)
        i_14 = np.argmin(np.abs(e1 - 0.14))
        s1, s3 = hv.sig[-1, 0], hv.sig[-1, 2]
        summary(f"  {psi:3.1f} {cell:.0f} | {e1[i_f]:.4f} {p[i_f]:6.1f} {qq[i_f]:6.1f} {du[i_f]:.3f} |"
                f" {du.max():.3f} {e1[np.argmax(du)]:.4f} | {p[i_14]:6.1f} {du[i_14]:.3f} |"
                f" {p[-1]:6.1f} {qq[-1]:6.1f} {du[-1]:.3f} {np.degrees(np.arcsin((s1 - s3) / (s1 + s3))):.3f}")
    for cell in cells:
        p_nocap = p_q(undrained_triaxial(props_of(hostun, M_cap=0.0), cell, 0.20, 800).sig)[0]
        summary(f"  without cap, psi = 0, {cell:.0f} kPa: p' between {p_nocap.min():.4f} and {p_nocap.max():.4f}")
    # sigma_c scaling: with c = 0 and power-law stiffness the model has no stress scale, so the response
    # normalised by sigma_c at 600 kPa is the one at 300 kPa with the strains multiplied by k = 2^(1 - m)
    k_model = 2.0 ** (1.0 - m)
    h3 = runs[(HOSTUN_PSI_VARIANT, 300.0)]
    h6 = undrained_triaxial(props_psi, 600.0, 0.20 * k_model, 800)
    dev = max(np.max(np.abs((600.0 - h6.sig[:, 2]) / 600.0 - (300.0 - h3.sig[:, 2]) / 300.0)),
              np.max(np.abs(h6.sig[:, 0] / h6.sig[:, 2] - h3.sig[:, 0] / h3.sig[:, 2])))
    summary(f"[sigma_c scaling] k = 2^(1-m) = {k_model:.4f}; UMAT (psi variant), 600 kPa at k eps_1 against"
            f" 300 kPa at eps_1: max difference in du/sc and sigma1'/sigma3' {dev:.1e}")

    summary("[hostun undrained step size] n_steps, p' and du/sc at the end (sigma_c = 300 kPa), psi = 0 and variant")
    for n in [200, 800, 3200]:
        out = []
        for pv in [props, props_psi]:
            hn = undrained_triaxial(pv, 300.0, 0.20, n)
            out.append(f"{p_q(hn.sig)[0][-1]:8.3f} {(300.0 - hn.sig[-1, 2]) / 300.0:.5f}")
        summary(f"  {n:5d} " + "  ".join(out))

    # --- step size (Hostun drained triaxial) ---
    fig, ax = plt.subplots()
    fine = drained_triaxial(props, 300.0, 0.18, 720)
    summary("[step size] n_steps, max |q_n - q_fine| / qf, max |ev_n - ev_fine|")
    for n, color, ls in zip([10, 40, 360], SERIES, ["-", (0, (6, 2)), (0, (1, 1.2))]):
        hn = drained_triaxial(props, 300.0, 0.18, n)
        ax.plot(hn.eps[:, 0], hn.sig[:, 0] / hn.sig[:, 2], color=color, ls=ls, marker="o" if n <= 10 else None,
                ms=4, mec="white", mew=0.8,
                label=rf"{n} steps ($\Delta\varepsilon_1$ = {0.18 / n:.1e})")
        qn = hn.sig[:, 0] - hn.sig[:, 2]
        q_fine = np.interp(hn.eps[:, 0], fine.eps[:, 0], fine.sig[:, 0] - fine.sig[:, 2])
        ev_fine = np.interp(hn.eps[:, 0], fine.eps[:, 0], fine.eps[:, :3].sum(axis=1))
        summary(f"  {n:4d} {np.max(np.abs(qn - q_fine)) / qf:.2e} {np.max(np.abs(hn.eps[:, :3].sum(axis=1) - ev_fine)):.2e}")
    ax.set_xlabel(r"axial strain $\varepsilon_1$ [-]")
    ax.set_ylabel(r"stress ratio $\sigma_1/\sigma_3$ [-]")
    ax.set_xlim(0, 0.18)
    ax.set_ylim(1, 4)
    ax.legend(loc="lower right")
    save(fig, "step_size")


def fig_general_stress_states():
    """Drained compression, plane strain and extension from 100 kPa: all reach Mohr-Coulomb."""
    props_cap = props_of(TEST)
    props_nocap = props_of(TEST, M_cap=0.0)
    cell = 100.0
    paths = {}
    # triaxial compression: sigma_y = sigma_z = cell
    paths["triaxial compression"] = drained_triaxial(props_cap, cell, 0.12, 120)
    # plane strain: eps_y = 0, sigma_z = cell (cap off, as in the unit test)
    stress, statev, strain = np.array([cell, cell, cell, 0, 0, 0]), np.zeros(3), np.zeros(6)
    h = History(strain, stress, statev)
    for _ in range(120):
        stress, statev, strain = mixed_step(stress, statev, strain, [1e-3, 0, 0, 0, 0, 0], [[2]],
                                            np.array([cell]), props_nocap)
        h.add(strain, stress, statev)
    paths["plane strain"] = h.arrays()
    # triaxial extension: sigma_x = sigma_y = cell, eps_z decreases
    stress, statev, strain = np.array([cell, cell, cell, 0, 0, 0]), np.zeros(3), np.zeros(6)
    h = History(strain, stress, statev)
    for _ in range(120):
        stress, statev, strain = mixed_step(stress, statev, strain, [0, 0, -1e-3, 0, 0, 0], [[0, 1]],
                                            np.array([cell]), props_cap)
        h.add(strain, stress, statev)
    paths["triaxial extension"] = h.arrays()

    sin_phi = np.sin(np.radians(TEST["phi"]))
    fig, ax = plt.subplots()
    summary("\n[stress states] path, (s1-s3)/(sin(phi)(s1+s3)) end, b end, at_failure flag")
    for (name, h), color, ls in zip(paths.items(), SERIES, ["-", (0, (6, 2)), (0, (1, 1.2))]):
        s = np.sort(h.sig[:, :3], axis=1)[:, ::-1]
        gamma = h.eps[:, :3].max(axis=1) - h.eps[:, :3].min(axis=1)
        mob = (s[:, 0] - s[:, 2]) / (sin_phi * (s[:, 0] + s[:, 2]))
        b = (s[-1, 1] - s[-1, 2]) / (s[-1, 0] - s[-1, 2])
        summary(f"  {name:22s} {mob[-1]:.8f} {b:.3f} {h.sv[-1, 2]:.0f}")
        ax.plot(100 * gamma, mob, color=color, ls=ls, label=f"{name} ($b$ = {b:.2f})")
    ax.axhline(1.0, **ENVELOPE)
    ax.set_xlabel(r"shear strain $\varepsilon_1 - \varepsilon_3$ [%]")
    ax.set_ylabel(r"$\sin\varphi_m / \sin\varphi$ [-]")
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 1.1)
    ax.legend(loc="lower right")
    save(fig, "stress_states_mobilisation")

    # deviatoric plane, stresses normalised by p (c = 0: the Mohr-Coulomb section is a fixed hexagon)
    def pi_plane(sig):
        p = sig[:, :3].mean(axis=1)
        x = (sig[:, 1] - sig[:, 2]) / np.sqrt(2) / p
        y = (2 * sig[:, 0] - sig[:, 1] - sig[:, 2]) / np.sqrt(6) / p
        return x, y

    # corners of the Mohr-Coulomb section at p = 1: triaxial compression (sigma_1 > sigma_2 = sigma_3)
    # and extension (sigma_1 = sigma_2 > sigma_3), each in the three principal directions
    tc, te = 2 * sin_phi / (3 - sin_phi), 2 * sin_phi / (3 + sin_phi)
    corners = []
    for i in range(3):
        compression, extension = np.full(3, 1 - tc), np.full(3, 1 + te)
        compression[i], extension[i] = 1 + 2 * tc, 1 - 2 * te
        corners += [compression, extension]
    corners = np.array(corners)
    cx, cy = pi_plane(np.c_[corners, np.zeros((len(corners), 3))])
    order = np.argsort(np.arctan2(cy, cx))
    cx, cy = np.r_[cx[order], cx[order][0]], np.r_[cy[order], cy[order][0]]

    fig, ax = plt.subplots()
    r = 1.15 * np.hypot(cx, cy).max()
    # projections of the principal axes: sigma_x up, sigma_y lower right, sigma_z lower left
    for ang, lab in [(90, r"$\sigma_x$"), (330, r"$\sigma_y$"), (210, r"$\sigma_z$")]:
        ax.plot([0, r * np.cos(np.radians(ang))], [0, r * np.sin(np.radians(ang))], color=AXIS, lw=0.8)
        ax.annotate(lab, (r * np.cos(np.radians(ang)), r * np.sin(np.radians(ang))), color=INK2,
                    ha="center", va="center", xytext=(0, 7 if ang == 90 else -7), textcoords="offset points")
    ax.plot(cx, cy, **ENVELOPE, label="Mohr-Coulomb\n($q = q_f$)")
    for (name, h), color, ls in zip(paths.items(), SERIES, ["-", (0, (6, 2)), (0, (1, 1.2))]):
        x, y = pi_plane(h.sig)
        ax.plot(x, y, color=color, ls=ls, label=name)
        ax.plot(x[-1], y[-1], "o", color=color, ms=5, mec="white", mew=0.8)
    ax.set_aspect("equal")
    ax.set_xlim(-r, r)
    ax.set_ylim(-0.75 * r, 1.12 * r)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7)
    save(fig, "stress_states_pi_plane")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    print(f"UMAT: {DLL}")
    fig_hyperbola()
    fig_dilatancy()
    fig_cap()
    fig_general_stress_states()

    (m_cap, k_ratio), res = calibrate_cap(HOSTUN, HOSTUN_K0_NC, HOSTUN_EOED_REF)
    k0, eoed = oedometer_response(props_of(HOSTUN, M_cap=m_cap, K_ratio=k_ratio))
    print(f"\n[calibration] M_cap = {m_cap:.4f}, K_ratio = {k_ratio:.4f} -> K0 = {k0:.4f}"
          f" (target {HOSTUN_K0_NC:.4f}), Eoed(p_ref) = {eoed:.1f} (target {HOSTUN_EOED_REF:.1f})")
    fig_hostun(m_cap, k_ratio)


if __name__ == "__main__":
    main()
