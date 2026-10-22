"""
Shared parts of the figure scripts of the test reports (docs/*_test_report/make_figures.py): calling
a UMAT shared library, element test bookkeeping, stress measures and the plot style of the reports.

A script imports it with

    sys.path.insert(0, os.path.dirname(HERE))  # docs
    import figure_utils as fu

Stresses and strains are compression positive in the helpers that take a `step` or `response`
function; `Umat` itself works in the tension-positive convention of the interface.
"""

import logging
import os
import sys

import cffi
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def library_path(name):
    """Shared library of a C model in build_C/lib (.dll on Windows, .so otherwise)."""
    return os.path.join(ROOT, "build_C", "lib", name + "." + ("dll" if sys.platform == "win32" else "so"))


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


class Umat:
    """The UMAT of one shared library, loaded once at the first call (same arguments as
    tests/utils.py, Utils.run_c_umat)."""

    def __init__(self, path, cmname="UMAT"):
        self.path = path
        self.cmname = cmname.encode().ljust(80)
        self.lib = None

    def __call__(self, stress, statev, strain, dstrain, props):
        """One call in the tension-positive convention of the interface. Returns the stress, the
        state variables, DDSDDE and the increments of SSE and SPD."""
        if self.lib is None:
            if not os.path.exists(self.path):
                sys.exit(f"{self.path} not found: build the C models first (see README.md)")
            self.lib = _ffi.dlopen(self.path)
        ints = [_ffi.new("int*", v) for v in (3, 3, 6, len(statev), len(props), 1, 1, 1, 1)]
        c_stress = _ffi.new("double[]", [float(v) for v in stress])
        c_statev = _ffi.new("double[]", [float(v) for v in statev])
        c_ddsdde = _ffi.new("double[]", 36)
        sse, spd, scd = [_ffi.new("double*", 0.0) for _ in range(3)]
        self.lib.umat(c_stress, c_statev, c_ddsdde, sse, spd, scd, _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL,
                      _ffi.new("double[]", [float(v) for v in strain]),
                      _ffi.new("double[]", [float(v) for v in dstrain]),
                      _ffi.new("double[]", [0.0, 0.0]), _ffi.new("double*", 1.0),
                      _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.new("char[]", self.cmname),
                      ints[0], ints[1], ints[2], ints[3], _ffi.new("double[]", [float(v) for v in props]), ints[4],
                      _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL, _ffi.NULL,
                      ints[5], ints[6], _ffi.NULL, _ffi.NULL, ints[7], ints[8])
        return (np.array(list(c_stress)), np.array(list(c_statev)), np.array(list(c_ddsdde)).reshape(6, 6),
                sse[0], spd[0])


# --------------------------------------------------------------------------- #
# Element test paths                                                            #
# --------------------------------------------------------------------------- #
class History:
    """Strain, stress and state variables along a path (lists while adding, arrays after
    `arrays`)."""

    def __init__(self, strain, stress, statev):
        self.eps, self.sig, self.sv = [], [], []
        self.add(strain, stress, statev)

    def add(self, strain, stress, statev):
        self.eps.append(np.array(strain, float))
        self.sig.append(np.array(stress, float))
        self.sv.append(np.array(statev, float))

    def arrays(self):
        self.eps, self.sig, self.sv = np.array(self.eps), np.array(self.sig), np.array(self.sv)
        return self


def mixed_strain_increment(response, dstrain, groups, targets, tol, max_iter=50, h=1e-8):
    """Strain increment in which each group of components (sharing one strain increment) is
    adjusted by Newton iteration so that the stress of its first component equals the target.

    `response(dstrain)` returns the stress and DDSDDE after the increment. The Jacobian is taken
    from DDSDDE, or by forward differences with step `h` when `response` returns None for it.
    Returns the strain increment."""
    d = np.array(dstrain, dtype=float)
    u = np.array([d[g[0]] for g in groups])
    first = [g[0] for g in groups]

    def apply(values):
        out = d.copy()
        for g, v in zip(groups, values):
            out[list(g)] = v
        return out

    for _ in range(max_iter):
        s, ddsdde = response(apply(u))
        r = s[first] - targets
        if np.max(np.abs(r)) < tol * (1.0 + np.max(np.abs(targets))):
            break
        if ddsdde is not None:
            jac = np.array([[ddsdde[i, list(g)].sum() for g in groups] for i in first])
        else:
            jac = np.zeros((len(groups), len(groups)))
            for a in range(len(groups)):
                u_h = u.copy()
                u_h[a] += h
                jac[:, a] = (response(apply(u_h))[0][first] - s[first]) / h
        u -= np.linalg.solve(jac, r)
    return apply(u)


# --------------------------------------------------------------------------- #
# Stress measures and elasticity                                                #
# --------------------------------------------------------------------------- #
def p_q(sig):
    """Mean stress and deviator stress of an (n, 6) array of stresses."""
    p = sig[:, :3].mean(axis=1)
    q = np.sqrt(0.5 * ((sig[:, 0] - sig[:, 1]) ** 2 + (sig[:, 1] - sig[:, 2]) ** 2
                       + (sig[:, 2] - sig[:, 0]) ** 2) + 3 * (sig[:, 3:] ** 2).sum(axis=1))
    return p, q


def elastic_stiffness(E, nu):
    """Isotropic elastic stiffness in Voigt notation with engineering shear strains."""
    G, lam = E / (2 * (1 + nu)), E * nu / ((1 + nu) * (1 - 2 * nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] += 2 * G
    D[3:, 3:] = G * np.eye(3)
    return D


def pi_plane(sig):
    """Projection of the normal stresses (sigma_x, sigma_y, sigma_z) of an (n, >= 3) array on the
    deviatoric plane, sigma_x up."""
    sig = np.atleast_2d(sig)
    return (sig[:, 1] - sig[:, 2]) / np.sqrt(2), (2 * sig[:, 0] - sig[:, 1] - sig[:, 2]) / np.sqrt(6)


def mohr_coulomb_section(phi):
    """Closed Mohr-Coulomb hexagon in the deviatoric plane at p + a = 1 (`pi_plane` coordinates):
    corners at triaxial compression (sigma_1 > sigma_2 = sigma_3) and extension
    (sigma_1 = sigma_2 > sigma_3), each in the three principal directions."""
    sin_phi = np.sin(np.radians(phi))
    tc, te = 2 * sin_phi / (3 - sin_phi), 2 * sin_phi / (3 + sin_phi)
    corners = []
    for i in range(3):
        compression, extension = np.full(3, 1 - tc), np.full(3, 1 + te)
        compression[i], extension[i] = 1 + 2 * tc, 1 - 2 * te
        corners += [compression, extension]
    cx, cy = pi_plane(np.array(corners))
    order = np.argsort(np.arctan2(cy, cx))
    return np.r_[cx[order], cx[order][0]], np.r_[cy[order], cy[order][0]]


# --------------------------------------------------------------------------- #
# Plot style (static PDF figures; palette validated for colour-vision          #
# deficiency, see the HS report)                                               #
# --------------------------------------------------------------------------- #
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
LINESTYLES = ["-", (0, (6, 2)), (0, (1, 1.2))]
REF = dict(color=INK, lw=0.9, ls=(0, (4, 2)))          # analytic solutions
ENVELOPE = dict(color=INK2, lw=0.9, ls=(0, (1, 1.5)))  # failure envelopes and limits


def apply_style():
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


def figure_saver(directory):
    """save(fig, name): writes <directory>/<name>.pdf and closes the figure."""
    os.makedirs(directory, exist_ok=True)

    def save(fig, name):
        fig.savefig(os.path.join(directory, name + ".pdf"))
        plt.close(fig)

    return save


def summary(text):
    """Values quoted in a report are printed."""
    print(text)


def marker(color, size=5.0):
    """Markers without line, for UMAT results at discrete points."""
    return dict(color=color, marker="o", ls="none", ms=size, mec="white", mew=0.8)


def log_axes(ax, xticks, yticks):
    """Log-log axes with plain tick labels."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    for axis, ticks in [(ax.xaxis, xticks), (ax.yaxis, yticks)]:
        axis.set_major_locator(matplotlib.ticker.FixedLocator(ticks))
        axis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())


def deviatoric_axes(ax, r):
    """Deviatoric-plane plot of radius r: projections of the principal axes (sigma_x up, sigma_y
    lower right, sigma_z lower left), equal aspect, no ticks, grid or spines."""
    for ang, lab in [(90, r"$\sigma_x$"), (330, r"$\sigma_y$"), (210, r"$\sigma_z$")]:
        ax.plot([0, r * np.cos(np.radians(ang))], [0, r * np.sin(np.radians(ang))], color=AXIS, lw=0.8)
        ax.annotate(lab, (r * np.cos(np.radians(ang)), r * np.sin(np.radians(ang))), color=INK2,
                    ha="center", va="center", xytext=(0, 7 if ang == 90 else -7), textcoords="offset points")
    ax.set_aspect("equal")
    ax.set_xlim(-r, r)
    ax.set_ylim(-0.75 * r, 1.12 * r)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
