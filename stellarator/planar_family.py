"""Exact quasisymmetric equilibria with no rigid-motion symmetry, in slab / tube geometry.

Found while designing the "symmetry by construction" recipe. The construction:

  1. Take a 2D field B0 in the (x, z) plane (no y component, divergence-free) whose field-line
     force B0 . grad B0 is the gradient of a function Pi0 that is unchanged by some group of
     motions of the plane (shifts in x, or rotations about the origin).
  2. In the plane at height y, use B0 moved by a motion g(y) from that group. Any g(y).

Then B . grad B = grad Pi0 everywhere, so B is in force balance with p = Pi0 - |B|^2 / 2, and
every plane carries the same field-strength pattern, moved. The vector field u that generates
"go up in y and apply the motion" is divergence-free, satisfies B x u = grad psi and
u . grad |B|^2 = 0: weak quasisymmetry, exactly. If g(y) is not a one-parameter group (shift
not linear in y, angle not linear in y), u is not a rigid motion and B has no rigid-motion
symmetry at all.

Two families:
  shift : B = exp(kappa z / 2) (cosh(mu X), 0, -(2 mu / kappa) sinh(mu X)),  X = x - m(y)
          open field lines, sheets z = c - (2 / kappa) log cosh(mu X)
  twist : B = R(y) M R(y)^T (x, z),  M = [[0, -2 beta], [2 alpha, 0]],  R(y) = rotation by theta(y)
          closed elliptical field lines around the y axis; nested elliptical tubes whose
          orientation turns with y at an arbitrary rate. A pinch with no field along its axis.
          With tan(theta / 2) = y the field is a rational function, which the Lean generator
          can certify.

What this is not: a toroidal example. The surfaces are sheets or infinite tubes (tori only if y
is made periodic, which is not a torus in ordinary space), and every field line lies in a
plane y = const, so there is no rotational transform. The open problem is untouched. What it
does show: force balance plus quasisymmetry do not force a rigid-motion symmetry near a point.

This script checks every equation symbolically (exact zero) and tests for a rigid-motion
symmetry at 40 digits. Usage:  python planar_family.py
"""
import random

import mpmath as mp
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
r = (x, y, z)


def grad(f):
    return sp.Matrix([sp.diff(f, v) for v in r])


def div(V):
    return sum(sp.diff(V[i], r[i]) for i in range(3))


def curl(V):
    return sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], z),
                      sp.diff(V[0], z) - sp.diff(V[2], x),
                      sp.diff(V[1], x) - sp.diff(V[0], y)])


def shift_family(m, kappa, mu):
    X = x - m
    e = sp.exp(kappa * z / 2)
    B = sp.Matrix([e * sp.cosh(mu * X), 0, -(2 * mu / kappa) * e * sp.sinh(mu * X)])
    psi = (2 / kappa) * e * sp.cosh(mu * X)
    p = -(mu ** 2 / 2 + kappa ** 2 / 8) * psi ** 2
    u = sp.Matrix([sp.diff(m, y), 1, 0])
    return B, u, psi, p


def twist_family(c, s, dtheta, alpha, beta):
    """c, s = cos and sin of the twist angle theta(y); dtheta = theta'(y)."""
    xp = c * x + s * z            # coordinates turned back by theta
    zp = -s * x + c * z
    b0 = (-2 * beta * zp, 2 * alpha * xp)
    B = sp.Matrix([c * b0[0] - s * b0[1], 0, s * b0[0] + c * b0[1]])
    psi = -(alpha * xp ** 2 + beta * zp ** 2)   # sign so that B x u = +grad psi
    p = 2 * (alpha + beta) * psi
    u = sp.Matrix([-dtheta * z, 1, dtheta * x])
    return B, u, psi, p


def equations(B, u, psi, p):
    return {
        "div B": div(B),
        "force balance": curl(B).cross(B) - grad(p),
        "div u": div(u),
        "B x u - grad psi": B.cross(u) - grad(psi),
        "u . grad |B|^2": u.dot(grad(B.dot(B))),
    }


def rigid_symmetry(B, npts=40, seed=0, box=0.5):
    """Smallest singular value of (a, w) -> L_K B at random points, K = a + w x r, relative to
    the largest. Zero iff some rigid motion leaves B unchanged."""
    a = sp.symbols("a1:4")
    w = sp.symbols("w1:4")
    K = sp.Matrix(a) + sp.Matrix(w).cross(sp.Matrix(r))
    L = sp.Matrix([sum(K[j] * sp.diff(B[i], r[j]) - B[j] * sp.diff(K[i], r[j]) for j in range(3))
                   for i in range(3)])
    f = sp.lambdify(r, L.jacobian(list(a) + list(w)), "mpmath")
    rnd = random.Random(seed)
    A = mp.matrix(3 * npts, 6)
    for n in range(npts):
        Mi = f(*[mp.mpf(rnd.uniform(-box, box)) for _ in range(3)])
        for i in range(3):
            for j in range(6):
                A[3 * n + i, j] = Mi[i, j]
    S = mp.svd_r(A, compute_uv=False)
    return S[5] / S[0]


def report(name, fields):
    B, u, psi, p = fields
    zero = all(sp.simplify(v) == (0 if not isinstance(v, sp.MatrixBase) else sp.zeros(*v.shape))
               for v in equations(B, u, psi, p).values())
    print(f"{name:46s} all equations exactly zero: {zero}   rigid-symmetry measure: "
          f"{mp.nstr(rigid_symmetry(B), 3)}")


if __name__ == "__main__":
    mp.mp.dps = 40
    q = sp.Rational
    print("A rigid-symmetry measure of order 1e-40 means symmetric; of order 0.1, not.\n")
    report("shift, m = 0.3 y (control, symmetric)", shift_family(q(3, 10) * y, q(7, 5), q(9, 10)))
    report("shift, m = y^2", shift_family(y ** 2, q(7, 5), q(9, 10)))
    report("shift, m = sin(2 y) / 2", shift_family(sp.sin(2 * y) / 2, q(7, 5), q(9, 10)))
    th = q(4, 5) * y
    report("twist, theta = 0.8 y (control, helical)", twist_family(sp.cos(th), sp.sin(th), q(4, 5), 1, q(1, 3)))
    th = sp.sin(y)
    report("twist, theta = sin y (periodic in y)", twist_family(sp.cos(th), sp.sin(th), sp.cos(y), 1, q(1, 3)))
    c, s = (1 - y ** 2) / (1 + y ** 2), 2 * y / (1 + y ** 2)
    report("twist, tan(theta/2) = y (rational field)", twist_family(c, s, 2 / (1 + y ** 2), 1, q(1, 3)))
    report("twist, alpha = beta (control, circular)", twist_family(c, s, 2 / (1 + y ** 2), 1, 1))
