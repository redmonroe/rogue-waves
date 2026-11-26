"""Quasisymmetry residuals for two published families that are not given as explicit
functions of (x, y, z).

1. Lopez, arXiv:2610.29921: a vacuum field B = Im grad log zeta, where zeta(x) is a selected
   root of a cubic whose coefficients depend on position. Derivatives are taken with the
   implicit-function rule grad zeta = N/D from the paper (eq. main-harmonic-chain).
2. Issan, Jose, Landreman, Citrin, Tracey, arXiv:2611.07672: equilibria with N field periods,
   given through a coordinate map x(rho, zeta, eta) with B = d x / d zeta. Derivatives are
   taken through the numerically inverted Jacobian.

Columns as in qs_residuals.py:
  check  : how well the transcription satisfies the paper's own equations (should be ~0)
  triple : grad(psi) x grad|B|^2 . grad(B . grad|B|^2), normalised, where a flux label is given
  frob   : N . curl N with N = B x (grad|B|^2 x grad(B . grad|B|^2)), normalised; must vanish
           for a weakly quasisymmetric field whatever flux label is used
Neither family has an axisymmetric member to use as a control; the `check` column and the
controls in qs_residuals.py are what validate the method.
"""
import random

import mpmath as mp
import sympy as sp

mp.mp.dps = 40
I = sp.I


def nrm(v):
    return mp.sqrt(sum(abs(c) ** 2 for c in v))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def dot(a, b):
    return sum(p * q for p, q in zip(a, b))


# ------------------------------------------------------------------ 1. Lopez vacuum field

def lopez(n_points=12, seed=0):
    x, y, z = sp.symbols("x y z", real=True)
    zt, zb = sp.symbols("zeta zetabar")  # zeta and its complex conjugate, treated as independent
    X = (x, y, z)
    w, wb = x + I * y, x - I * y
    N = [1 - zt**2, I * (1 + zt**2), 2 * zt]
    D = 6 * I * zt**2 + 2 * (wb + 2) * zt - 2 * z
    Nb = [1 - zb**2, -I * (1 + zb**2), 2 * zb]
    Db = -6 * I * zb**2 + 2 * (w + 2) * zb - 2 * z

    def d(e, i):  # total derivative with respect to X[i], using grad zeta = N / D
        return sp.diff(e, X[i]) + sp.diff(e, zt) * N[i] / D + sp.diff(e, zb) * Nb[i] / Db

    grad = lambda e: [d(e, i) for i in range(3)]
    curl = lambda V: [d(V[2], 1) - d(V[1], 2), d(V[0], 2) - d(V[2], 0), d(V[1], 0) - d(V[0], 1)]
    B = [(N[j] / (zt * D) - Nb[j] / (zb * Db)) / (2 * I) for j in range(3)]
    divB = sum(d(B[j], j) for j in range(3))
    B2 = dot(B, B)
    g = grad(B2)
    T = dot(B, g)
    Nq = cross(B, cross(g, grad(T)))
    args = (x, y, z, zt, zb)
    f = sp.lambdify(args, [[divB] + curl(B), Nq, curl(Nq), B], "mpmath", cse=True)
    Gpoly = sp.lambdify(args, 2 * I * zt**3 + (wb + 2) * zt**2 - 2 * z * zt - w, "mpmath")

    rng = random.Random(seed)
    worst = dict(check=mp.mpf(0), frob=mp.mpf(0))
    lam0, beta0 = mp.mpf(3) / 20, mp.mpf(1) / 50
    for _ in range(n_points):
        r, th, phi = mp.sqrt(rng.random()) * 0.95, rng.uniform(0, 2 * mp.pi), rng.uniform(0, 2 * mp.pi)
        lam, beta = lam0 + r * mp.cos(th) / 30, beta0 + r * mp.sin(th) / 10  # inside the ellipse Q < 1
        z0 = mp.sqrt(lam) * mp.expj(phi)
        u, v = z0.real, z0.imag
        Nv = [1 - z0**2, 1j * (1 + z0**2), 2 * z0]
        A, C = [c.real for c in Nv], [c.imag for c in Nv]
        H = 2 * z0**2 + 2j * z0**3
        n = [2 * u / (1 + lam), 2 * v / (1 + lam), (lam - 1) / (1 + lam)]
        kappa = 4 * lam * u * (1 - v) / (1 + lam) ** 2
        pt = [(H.real * A[j] + H.imag * C[j]) / (1 + lam) ** 2 + (beta - kappa) * n[j] for j in range(3)]
        a = (*pt, z0, mp.conj(z0))
        chk, Nn, cN, Bv = f(*a)
        worst["check"] = max(worst["check"], abs(Gpoly(*a)), max(abs(c) for c in chk) / nrm(Bv))
        Nn, cN = [mp.re(c) for c in Nn], [mp.re(c) for c in cN]
        if nrm(Nn) > 0 and nrm(cN) > 0:
            worst["frob"] = max(worst["frob"], abs(dot(Nn, cN)) / (nrm(Nn) * nrm(cN)))
    return worst


# ------------------------------------------------------------------ 2. N field periods

def issan_nfp(Nfp, k, rho_a, sgn=1, n_points=10, seed=0, spread=0.03):
    rho, zeta, eta = sp.symbols("rho zeta eta", real=True)
    q = (rho, zeta, eta)
    rho0 = 1
    m1, m2, m3 = 1 + (k - 1) * Nfp, 1 + k * Nfp, 1 + (k + 1) * Nfp
    C1 = sp.sign(m1) * rho0 / sp.Integer(m1) ** 2 * (rho / rho0) ** m1
    C3 = rho0 / sp.Integer(m3) ** 2 * (rho0 / rho) ** m3
    C5 = 2 * rho0 / (Nfp * sp.sqrt(abs(m1) * m3)) * (rho0 / rho) ** Nfp
    Dp = m1 * (m1 + m2) * C1 + m3 * (m3 + m2) * C3
    Dm = m1 * (m1 + m2) * C1 - m3 * (m3 + m2) * C3
    at = lambda e: e.subs(rho, rho_a)
    rad = (m1 * at(C1) + m3 * at(C3)) * (m1**2 * at(C1) - m3**2 * at(C3)) * at(Dm) / (m2**3 * at(Dp))
    if not (rad > 0):
        return None  # C2(rho_a) is not real for these parameters
    C2 = sgn * sp.sqrt(rad) * (rho_a / rho) ** Nfp * at(Dm) / Dm
    C4 = Nfp**2 * C5 / Dp
    xy = C1 * sp.exp(I * m1 * zeta) + (C2 + I * eta * C4) * sp.exp(I * m2 * zeta) - C3 * sp.exp(I * m3 * zeta)
    re = lambda e: sp.simplify((e + e.subs(I, -I)) / 2)
    im = lambda e: sp.simplify((e - e.subs(I, -I)) / (2 * I))
    pos = sp.Matrix([re(xy), im(xy), eta + C5 * sp.sin(Nfp * zeta)])
    B = pos.diff(zeta)
    Hf = (m1 * C1 + m3 * C3) ** 2 + m2**2 * (C2**2 + eta**2 * C4**2)
    pr = -Hf  # p = p0 - 2 psi, psi = (H - H_axis)/2
    B2 = (B.T * B)[0]
    T = sp.diff(B2, zeta)  # B . grad = d/d zeta, because B = d x / d zeta
    dq = lambda e: [sp.diff(e, v) for v in q]
    f = sp.lambdify(q, [pos.jacobian(q).tolist(), list(B), [dq(B[i]) for i in range(3)],
                        dq(pr), dq(B2), dq(T), dq(Hf)], "mpmath", cse=True)

    def pieces(pt):
        M, Bv, dB, dp, dB2, dT, dH = f(*pt)
        Minv = mp.matrix(M) ** -1
        gradq = lambda dv: [sum(dv[kk] * Minv[kk, j] for kk in range(3)) for j in range(3)]
        return Bv, [gradq(dB[i]) for i in range(3)], gradq(dp), gradq(dB2), gradq(dT), gradq(dH), Minv

    def Nvec(pt):
        Bv, _, _, g, gT, _, _ = pieces(pt)
        return cross(Bv, cross(g, gT))

    rng = random.Random(seed)
    worst = dict(check=mp.mpf(0), triple=mp.mpf(0), frob=mp.mpf(0))
    for _ in range(n_points):
        pt = [mp.mpf(rho_a) * (1 + rng.uniform(-spread, spread)), mp.mpf(rng.uniform(0, 2 * mp.pi)),
              mp.mpf(rng.uniform(-spread, spread))]
        Bv, J, gp, g, gT, gH, Minv = pieces(pt)  # J[i][j] = d B_i / d x_j
        divB = J[0][0] + J[1][1] + J[2][2]
        curlB = [J[2][1] - J[1][2], J[0][2] - J[2][0], J[1][0] - J[0][1]]
        force = cross(curlB, Bv)
        scale = nrm(Bv) * nrm(curlB) + nrm(gp)
        worst["check"] = max(worst["check"], abs(divB) / nrm(Bv), nrm([force[i] - gp[i] for i in range(3)]) / scale)
        sc = nrm(gH) * nrm(g) * nrm(gT)
        if sc > 0:
            worst["triple"] = max(worst["triple"], abs(dot(gH, cross(g, gT))) / sc)
        Nn = cross(Bv, cross(g, gT))
        dN = [[mp.diff(lambda t, i=i, kk=kk: Nvec([t if m == kk else pt[m] for m in range(3)])[i], pt[kk])
               for kk in range(3)] for i in range(3)]
        JN = [[sum(dN[i][kk] * Minv[kk, j] for kk in range(3)) for j in range(3)] for i in range(3)]
        cN = [JN[2][1] - JN[1][2], JN[0][2] - JN[2][0], JN[1][0] - JN[0][1]]
        if nrm(Nn) > 0 and nrm(cN) > 0:
            worst["frob"] = max(worst["frob"], abs(dot(Nn, cN)) / (nrm(Nn) * nrm(cN)))
    return worst


if __name__ == "__main__":
    import time
    print(f"{'field':58s} {'check':>10s} {'triple':>10s} {'frob':>10s}")
    t0 = time.time()
    w = lopez()
    print(f"{'Lopez vacuum field (12 points in the torus)':58s} {float(w['check']):10.1e} {'n/a':>10s} {float(w['frob']):10.1e}   [{time.time() - t0:.0f}s]", flush=True)
    R = sp.Rational
    for Nfp, k, rho_a in [(1, 1, R(3, 2)), (2, 1, R(3, 2)), (3, 1, R(13, 10)), (2, 0, 1), (3, 0, R(1, 2)), (3, -1, R(4, 5)), (5, -1, R(4, 5))]:
        t0 = time.time()
        w = issan_nfp(Nfp, k, rho_a)
        name = f"Issan N={Nfp}, k={k}, rho_a={rho_a}"
        if w is None:
            print(f"{name:58s} parameters outside the family's valid range", flush=True)
        else:
            print(f"{name:58s} {float(w['check']):10.1e} {float(w['triple']):10.1e} {float(w['frob']):10.1e}   [{time.time() - t0:.0f}s]", flush=True)
