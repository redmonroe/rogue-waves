"""Published closed-form 3D MHD equilibria, as SymPy expressions in Cartesian x, y, z.

Each function returns (B, psi, p, in_domain) where B is a 3x1 Matrix, psi a flux label,
p the pressure, and in_domain(point, params) says whether a numeric point is usable.
Formulas are transcribed from Landreman, arXiv:2609.26742 (v2), and checked against the
SymPy notebooks in his supplemental repository (github.com/landreman/analytic_3d_equilibria).
"""
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
X = (x, y, z)


def _re(e):
    return (e + e.subs(sp.I, -sp.I)) / 2


def _im(e):
    return (e - e.subs(sp.I, -sp.I)) / (2 * sp.I)


def landreman_iota2(a, b):
    """Eqs. (2.1)-(2.4). The paper's family is a = sqrt(1+eps), b = sqrt(1-eps); a = b = 1 is axisymmetric."""
    s = x**2 / a**2 + y**2 / b**2
    G = 1 - (1 - s) ** 2 - 4 * z**2
    F = sp.sqrt(G)
    B = sp.Matrix([(2 * z * x - (a / b) * F * y) / s, (2 * z * y + (b / a) * F * x) / s, 1 - s])
    B2 = (B.T * B)[0]
    psi = (x**2 + y**2 + 4 * z**2 + B2 - 2) / 4  # constant offset dropped
    p = -2 * psi
    G_fn = sp.lambdify(X, G, "mpmath")
    s_fn = sp.lambdify(X, s, "mpmath")
    return B, psi, p, lambda pt: G_fn(*pt) > 0.05 and s_fn(*pt) > 0.05


def landreman_sheared(eps, lam, S):
    """Eqs. (3.1)-(3.3). eps = 0 is axisymmetric."""
    w, wb = x + sp.I * y, x - sp.I * y
    K = wb * sp.sqrt(1 + eps / wb**2)
    Xi = w * K + sp.pi / 2 - S
    h = sp.I * sp.exp(-sp.I * lam * z) * sp.sin(Xi) / (2 * K)
    c = _re(sp.exp(-sp.I * lam * z) * sp.cos(Xi))
    B = sp.Matrix([_re(h), _im(h), c / lam])
    psi = (sp.sin(lam * z) ** 2 + c**2) / 2
    p = -psi / lam**2
    rad_fn = sp.lambdify(X, 1 + eps / wb**2, "mpmath")
    r2_fn = sp.lambdify(X, x**2 + y**2, "mpmath")

    def in_domain(pt):
        rad = rad_fn(*pt)
        return r2_fn(*pt) > 0.05 and rad.real > 0.05  # stay off the square-root branch cut

    return B, psi, p, in_domain


def issan_iota2(a, b, tau0):
    """Issan, Citrin, Landreman, Tracey, arXiv:2610.07304, eqs. (f1_s_F)-(f1_pressure_psi):
    Landreman's iota=2 family with a stellarator-symmetry-breaking parameter tau0.
    tau0 = 0 recovers landreman_iota2; a = b is axisymmetric."""
    s = x**2 / a**2 + y**2 / b**2
    G = 1 - (1 - s) ** 2 - 4 * z**2 - 4 * tau0 * s * z
    F = sp.sqrt(G)
    B = sp.Matrix([((2 * z + tau0 * s) * x - (a / b) * F * y) / s,
                   ((2 * z + tau0 * s) * y + (b / a) * F * x) / s,
                   1 - s - 2 * tau0 * z])
    Om2 = 1 - tau0**2
    Q = Om2 / 2 * (x**2 + y**2 + 4 * z**2) + 2 * tau0 * z
    B2 = (B.T * B)[0]
    psi = Om2 / 2 * (Q + B2 / 2)  # constant offset dropped
    p = -2 * psi / Om2
    G_fn = sp.lambdify(X, G, "mpmath")
    s_fn = sp.lambdify(X, s, "mpmath")
    return B, psi, p, lambda pt: G_fn(*pt) > 0.05 and s_fn(*pt) > 0.05


def issan_sheared(eps, lam, S, tau1):
    """Issan et al., arXiv:2610.07304, eqs. (f2_w_K_W_Xi_def)-(f2_Q_pressure_psi):
    Landreman's sheared family with a symmetry-breaking phase tau1. tau1 = 0 recovers
    landreman_sheared; eps = 0 is axisymmetric."""
    w, wb = x + sp.I * y, x - sp.I * y
    K = wb * sp.sqrt(1 + eps / wb**2)
    Xi = w * K + sp.pi / 2 - S + sp.I * tau1
    h = sp.exp(-sp.I * lam * z) * sp.I * sp.sin(Xi) / (2 * K)
    Bz = _re(sp.exp(-sp.I * lam * z) * sp.cos(Xi)) / lam
    B = sp.Matrix([_re(h), _im(h), Bz])
    psi = (sp.sin(lam * z) ** 2 + lam**2 * Bz**2) / 2
    p = -psi / lam**2
    rad_fn = sp.lambdify(X, 1 + eps / wb**2, "mpmath")
    r2_fn = sp.lambdify(X, x**2 + y**2, "mpmath")

    def in_domain(pt):
        return r2_fn(*pt) > 0.05 and rad_fn(*pt).real > 0.05

    return B, psi, p, in_domain
