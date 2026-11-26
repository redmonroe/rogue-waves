"""Local (Taylor-series) test: can a field be in force balance and quasisymmetric near a point
without having an ordinary Euclidean symmetry?

Unknowns: Taylor coefficients at the origin of B, u (vector fields), psi, p (scalars).
Equations, imposed on every Taylor coefficient up to degree K-1:

    div B = 0                    (curl B) x B = grad p          (equilibrium)
    div u = 0     B x u = grad psi     u . grad|B|^2 = 0        (weak quasisymmetry)

Normalisation (uses up rotations and the three scalings, and keeps away from the degenerate
solutions with B . grad|B| = 0, where quasisymmetry says nothing):

    B(0) = (0, 0, 1),   u(0) = (1, 0, c),   dB_z/dz (0) = 1,   psi(0) = p(0) = 0.

A degree-K solution is found by damped Gauss-Newton from a random start. For each solution two
gauge-independent numbers are reported:

  sym_B   : smallest singular value of the map (a, w) -> low-order Taylor coefficients of the Lie
            derivative of B along the rigid motion a + w x X. Zero iff B is invariant, to that
            order, under some translation / rotation / screw motion.
  kill_u  : the part of the symmetrised gradient of u at the origin that cannot be removed by
            rescaling u -> h(psi) u. Zero iff u is, to first order, a function of psi times a
            Killing field.

A truncated solution constrains its low-order coefficients only through the equations up to
degree K-1, so what matters is the trend as K grows. If sym_B stays of order one for random
solutions at every K, there is no local obstruction through that order. If it collapses to zero
beyond some K, quasisymmetry plus force balance is forcing a Euclidean symmetry locally.
Neither outcome is a proof: this is a numerical experiment on truncated series.

Usage:  python local_jets.py [Kmax] [starts]
"""
import sys
import time

import numpy as np


class Jets:
    """Truncated power series in (x, y, z), stored as coefficient vectors over the monomials of
    total degree <= K+1. Products and derivatives use precomputed index tables."""

    def __init__(self, K):
        self.K = K
        D = K + 1
        mons = [(i, j, d - i - j) for d in range(D + 1) for i in range(d, -1, -1) for j in range(d - i, -1, -1)]
        self.mons = mons
        self.N = len(mons)
        pos = {m: n for n, m in enumerate(mons)}
        self.pos = pos
        self.deg = np.array([sum(m) for m in mons])
        # product table: out[O] += a[I] * b[J] for every pair of monomials with degree sum <= D
        I, J, O = [], [], []
        for n1, m1 in enumerate(mons):
            for n2, m2 in enumerate(mons):
                if self.deg[n1] + self.deg[n2] <= D:
                    I.append(n1); J.append(n2); O.append(pos[tuple(p + q for p, q in zip(m1, m2))])
        self.I, self.J, self.O = np.array(I), np.array(J), np.array(O)
        # derivative tables: (d a / d x_ax)[to] = factor * a[frm]
        self.dtab = []
        for ax in range(3):
            frm, to, fac = [], [], []
            for n, m in enumerate(mons):
                if m[ax] > 0:
                    lower = list(m); lower[ax] -= 1
                    frm.append(n); to.append(pos[tuple(lower)]); fac.append(m[ax])
            self.dtab.append((np.array(frm), np.array(to), np.array(fac, dtype=float)))
        # unknown layout: B and u up to degree K, psi and p in degrees 1 .. K+1
        self.masks = [self.deg <= K] * 6 + [(self.deg >= 1) & (self.deg <= K + 1)] * 2
        self.sizes = [int(m.sum()) for m in self.masks]
        self.offsets = np.concatenate([[0], np.cumsum(self.sizes)])
        self.nvar = int(self.offsets[-1])
        self.eq_mask = self.deg <= K - 1
        self.i0 = pos[(0, 0, 0)]
        self.lin = [pos[tuple(1 if t == ax else 0 for t in range(3))] for ax in range(3)]

    # -- truncated power-series arithmetic
    def mul(self, a, b):
        w = a[self.I] * b[self.J]
        if np.iscomplexobj(w):  # real and imaginary parts are kept separate for the complex-step Jacobian
            return np.bincount(self.O, w.real, self.N) + 1j * np.bincount(self.O, w.imag, self.N)
        return np.bincount(self.O, w, self.N)

    def d(self, a, ax):
        frm, to, fac = self.dtab[ax]
        out = np.zeros_like(a)
        out[to] = a[frm] * fac
        return out

    def unpack(self, x):
        fields = []
        for m, o, s in zip(self.masks, self.offsets, self.sizes):
            a = np.zeros(self.N, dtype=x.dtype)
            a[m] = x[o:o + s]
            fields.append(a)
        return fields[0:3], fields[3:6], fields[6], fields[7]

    def cross(self, a, b):
        m = self.mul
        return [m(a[1], b[2]) - m(a[2], b[1]), m(a[2], b[0]) - m(a[0], b[2]), m(a[0], b[1]) - m(a[1], b[0])]

    def curl(self, V):
        d = self.d
        return [d(V[2], 1) - d(V[1], 2), d(V[0], 2) - d(V[2], 0), d(V[1], 0) - d(V[0], 1)]

    def residual(self, x):
        B, u, psi, p = self.unpack(x)
        d, m = self.d, self.mul
        res = [sum(d(B[i], i) for i in range(3)), sum(d(u[i], i) for i in range(3))]
        f = self.cross(self.curl(B), B)
        res += [f[i] - d(p, i) for i in range(3)]
        g = self.cross(B, u)
        res += [g[i] - d(psi, i) for i in range(3)]
        B2 = sum(m(B[i], B[i]) for i in range(3))
        res.append(sum(m(u[i], d(B2, i)) for i in range(3)))
        out = [r[self.eq_mask] for r in res]
        i0 = self.i0
        pins = [B[0][i0], B[1][i0], B[2][i0] - 1, u[0][i0] - 1, u[1][i0], B[2][self.lin[2]] - 1]
        return np.concatenate(out + [np.array(pins, dtype=x.dtype)])

    def jacobian(self, x):
        h = 1e-30  # complex step: exact for polynomial residuals
        xc = x.astype(complex)
        cols = []
        for j in range(self.nvar):
            xc[j] += 1j * h
            cols.append(self.residual(xc).imag / h)
            xc[j] = x[j]
        return np.array(cols).T

    # -- solver
    def solve(self, rng, scale=0.5, tol=1e-12, max_iter=60):
        x = np.concatenate([rng.normal(size=s) * (scale ** self.deg[m]) for m, s in zip(self.masks, self.sizes)])
        r = self.residual(x)
        for it in range(max_iter):
            nr = np.linalg.norm(r)
            if nr < tol:
                return x, nr, it
            J = self.jacobian(x)
            step = np.linalg.lstsq(J, r, rcond=1e-12)[0]
            t = 1.0
            while t > 1e-6:
                xn = x - t * step
                rn = self.residual(xn)
                if np.linalg.norm(rn) < (1 - 0.3 * t) * nr:
                    break
                t /= 2
            x, r = xn, rn
        return x, np.linalg.norm(r), max_iter

    # -- diagnostics
    def sym_B(self, x, dmax):
        """Smallest singular value of (a, w) -> Taylor coefficients of L_K B up to degree dmax,
        K = a + w x X. Zero iff B is invariant, to that order, under some rigid motion."""
        B, _, _, _ = self.unpack(x)
        X = []
        for ax in range(3):
            a = np.zeros(self.N)
            a[self.lin[ax]] = 1.0
            X.append(a)
        one = np.zeros(self.N)
        one[self.i0] = 1.0
        mask = self.deg <= dmax
        cols = []
        for q in range(6):
            a, w = np.zeros(3), np.zeros(3)
            (a if q < 3 else w)[q % 3] = 1.0
            Kf = [a[i] * one + (w[(i + 1) % 3] * X[(i + 2) % 3] - w[(i + 2) % 3] * X[(i + 1) % 3]) for i in range(3)]
            wxB = [w[(i + 1) % 3] * B[(i + 2) % 3] - w[(i + 2) % 3] * B[(i + 1) % 3] for i in range(3)]
            L = [sum(self.mul(Kf[j], self.d(B[i], j)) for j in range(3)) - wxB[i] for i in range(3)]
            cols.append(np.concatenate([l[mask] for l in L]))
        return float(np.linalg.svd(np.array(cols).T, compute_uv=False)[-1])

    def kill_u(self, x):
        _, u, _, _ = self.unpack(x)
        G = np.array([[u[j][self.lin[i]] for j in range(3)] for i in range(3)])
        S = G + G.T  # S_ij = d_i u_j + d_j u_i at the origin
        c = u[2][self.i0]
        M = np.zeros((3, 3))  # grad psi (x) u + u (x) grad psi at the origin, grad psi = e_y, u = (1, 0, c)
        M[0, 1] = M[1, 0] = 1.0
        M[1, 2] = M[2, 1] = c
        kappa = np.sum(S * M) / np.sum(M * M)
        return float(np.linalg.norm(S - kappa * M))


def check_exact_solution():
    """Residuals must vanish on the screw pinch B = (-y, x, 1), p = -(x^2 + y^2), u = e_z,
    psi = (x^2 + y^2)/2 (a degenerate exact solution, used only to test the arithmetic)."""
    J = Jets(4)
    B, u, psi, p = J.unpack(np.zeros(J.nvar))
    P = J.pos
    B[0][P[(0, 1, 0)]], B[1][P[(1, 0, 0)]], B[2][P[(0, 0, 0)]] = -1.0, 1.0, 1.0
    u[2][P[(0, 0, 0)]] = 1.0
    psi[P[(2, 0, 0)]] = psi[P[(0, 2, 0)]] = 0.5
    p[P[(2, 0, 0)]] = p[P[(0, 2, 0)]] = -1.0
    x = np.concatenate([f[m] for f, m in zip(B + u + [psi, p], J.masks)])
    r = J.residual(x)[:-6]  # ignore the normalisation pins
    assert np.abs(r).max() < 1e-14, np.abs(r).max()
    assert J.sym_B(x, 2) < 1e-12  # this field has translation and rotation symmetry


if __name__ == "__main__":
    check_exact_solution()
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    starts = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    print(f"{'K':>2s} {'unknowns':>8s} {'equations':>9s} {'solved':>7s} {'free dims':>9s} "
          f"{'sym_B deg<=1 (min / median)':>28s} {'sym_B deg<=2 (min / median)':>28s} {'kill_u (min / median)':>22s} {'time':>6s}")
    for K in range(2, Kmax + 1):
        J = Jets(K)
        rng = np.random.default_rng(K)
        t0 = time.time()
        s1, s2, ku, null = [], [], [], []
        for _ in range(starts):
            x, nr, it = J.solve(rng)
            if nr > 1e-10:
                continue
            s1.append(J.sym_B(x, 1))
            s2.append(J.sym_B(x, min(2, K - 1)))
            ku.append(J.kill_u(x))
            sv = np.linalg.svd(J.jacobian(x), compute_uv=False)
            null.append(J.nvar - int((sv > 1e-9 * sv[0]).sum()))
        neq = len(J.residual(np.zeros(J.nvar)))
        if s1:
            fmt = lambda v: f"{min(v):12.2e} / {np.median(v):9.2e}"
            print(f"{K:2d} {J.nvar:8d} {neq:9d} {len(s1):3d}/{starts:<3d} {int(np.median(null)):9d} "
                  f"{fmt(s1):>28s} {fmt(s2):>28s} {fmt(ku):>22s} {time.time() - t0:5.0f}s", flush=True)
        else:
            print(f"{K:2d} {J.nvar:8d} {neq:9d}   0/{starts:<3d}  no solutions found {time.time() - t0:5.0f}s", flush=True)
