"""Bottom-right cell, cleanest case: vacuum fields with strong quasisymmetry, near a point.

In a vacuum B = grad Phi with Phi harmonic, and the strong form of quasisymmetry reduces (by
hand, see LOG.md "Reduction of the strong form") to two equations for two functions and one
constant C:

    Laplace(Phi) = 0
    grad psi . grad Phi = 0                                         (V1)
    (C grad Phi + grad psi x grad Phi) . grad |grad Phi|^2 = 0      (V2)

The symmetry field is u = (C grad Phi + grad psi x grad Phi) / |grad Phi|^2; then div u = 0,
B x u = grad psi and u . B = C hold identically. Unknowns here: Taylor coefficients of Phi
(degrees 1..K+1) and psi (degrees 1..K), and C: about 500 numbers at order 9, against 1700 in
local_jets.py. So solves from random starts directly at high order are cheap.

Normalisation: B(0) = (0, 0, 1), grad psi(0) = (1, 0, 0), dB_z/dz(0) = 1 (so |B| varies along B).

For each converged solution: sym_B (zero iff some rigid motion leaves B unchanged, to second
order) and strain (norm of the symmetrised gradient of u at the origin; zero iff u is a rigid
motion to first order).

Usage:  python local_vacuum.py [--kmin 5] [--kmax 10] [--seeds 0:8] [--steps 200] [--check] [--ladder]
        --check   runs the equations on an exact translation-symmetric solution first.
        --ladder  solves at order kmin from random starts, then carries each solution up to kmax.
        Set OMP_NUM_THREADS=1: the seeds already run in parallel.
"""
import argparse
import multiprocessing
import time

import numpy as np

from local_jets import Jets


class Vac:
    def __init__(self, K):
        self.K = K
        J = self.J = Jets(K)
        self.mPhi = (J.deg >= 1) & (J.deg <= K + 1)
        self.mpsi = (J.deg >= 1) & (J.deg <= K)
        self.nPhi, self.npsi = int(self.mPhi.sum()), int(self.mpsi.sum())
        self.nvar = self.nPhi + self.npsi + 1
        self.eq = J.deg <= K - 1
        self.var_deg = np.concatenate([J.deg[self.mPhi], J.deg[self.mpsi], [0]])

    def unpack(self, x):
        J = self.J
        Phi = np.zeros(J.N, dtype=x.dtype)
        psi = np.zeros(J.N, dtype=x.dtype)
        Phi[self.mPhi] = x[:self.nPhi]
        psi[self.mpsi] = x[self.nPhi:self.nPhi + self.npsi]
        return Phi, psi, x[-1]

    def fields(self, x):
        """B, grad psi, and n^2 u = C B + grad psi x B."""
        J = self.J
        Phi, psi, C = self.unpack(x)
        B = [J.d(Phi, i) for i in range(3)]
        g = [J.d(psi, i) for i in range(3)]
        gxB = J.cross(g, B)
        return B, g, [C * B[i] + gxB[i] for i in range(3)]

    def residual(self, x, pins=True):
        J = self.J
        m, d = J.mul, J.d
        Phi, psi, C = self.unpack(x)
        B, g, w = self.fields(x)
        lap = sum(d(B[i], i) for i in range(3))
        v1 = sum(m(g[i], B[i]) for i in range(3))
        N = sum(m(B[i], B[i]) for i in range(3))
        v2 = sum(m(w[i], d(N, i)) for i in range(3))
        out = [lap[self.eq], v1[self.eq], v2[self.eq]]
        if pins:
            i0, lin = J.i0, J.lin
            zz = J.pos[(0, 0, 2)]
            out.append(np.array([B[0][i0], B[1][i0], B[2][i0] - 1, g[0][i0] - 1, g[1][i0], 2 * Phi[zz] - 1],
                                dtype=x.dtype))
        return np.concatenate(out)

    def jacobian(self, x):
        xc = x.astype(complex)
        cols = []
        for j in range(self.nvar):
            xc[j] += 1e-30j
            cols.append(self.residual(xc).imag / 1e-30)
            xc[j] = x[j]
        return np.array(cols).T

    def solve(self, x, max_steps=200, log=None):
        """Minimum-norm Gauss-Newton with backtracking (as local_extend.py)."""
        r = self.residual(x)
        nr = np.linalg.norm(r)
        for it in range(1, max_steps + 1):
            U, s, Vt = np.linalg.svd(self.jacobian(x), full_matrices=False)
            Utr = U.T @ r
            found = False
            for cut in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4):
                keep = s > cut * s[0]
                full = Vt[keep].T @ (Utr[keep] / s[keep])
                t = 1.0
                while t > 1e-3:
                    xn = x - t * full
                    rn = self.residual(xn)
                    nn = np.linalg.norm(rn)
                    if nn < (1 - 0.3 * t) * nr:
                        found = True
                        break
                    t /= 2
                if found:
                    break
            size = max(1.0, np.abs(x).max())
            if not found:
                return x, nr, it, nr < 1e-11 * size
            x, r, nr = xn, rn, nn
            if log:
                log(f"   step {it:3d} residual {nr:9.3e} cut {cut:5.0e} length {t:4.2f}")
            if nr < 1e-12 * max(1.0, np.abs(x).max()):
                return x, nr, it, True
        return x, nr, max_steps, False

    # -- diagnostics
    def sym_B(self, x, dmax=2):
        J = self.J
        B, _, _ = self.fields(x)
        X = []
        for ax in range(3):
            a = np.zeros(J.N)
            a[J.lin[ax]] = 1.0
            X.append(a)
        one = np.zeros(J.N)
        one[J.i0] = 1.0
        mask = J.deg <= dmax
        cols = []
        for q in range(6):
            a, w = np.zeros(3), np.zeros(3)
            (a if q < 3 else w)[q % 3] = 1.0
            Kf = [a[i] * one + (w[(i + 1) % 3] * X[(i + 2) % 3] - w[(i + 2) % 3] * X[(i + 1) % 3]) for i in range(3)]
            wxB = [w[(i + 1) % 3] * B[(i + 2) % 3] - w[(i + 2) % 3] * B[(i + 1) % 3] for i in range(3)]
            L = [sum(J.mul(Kf[j], J.d(B[i], j)) for j in range(3)) - wxB[i] for i in range(3)]
            cols.append(np.concatenate([l[mask] for l in L]))
        return float(np.linalg.svd(np.array(cols).T, compute_uv=False)[-1])

    def strain(self, x):
        """Norm of the symmetrised gradient of u at the origin."""
        J = self.J
        B, _, w = self.fields(x)
        N = sum(J.mul(B[i], B[i]) for i in range(3))
        one = np.zeros(J.N)
        one[J.i0] = 1.0
        inv = one / N[J.i0]
        for _ in range(3):  # reciprocal of a power series, good to first order after one pass
            inv = J.mul(inv, 2 * one - J.mul(N, inv))
        u = [J.mul(w[i], inv) for i in range(3)]
        G = np.array([[u[j][J.lin[i]] for j in range(3)] for i in range(3)])
        return float(np.linalg.norm(G + G.T) / 2)


def exact_check(K=8):
    """Translation-symmetric vacuum field: Phi = c y + chi(x, z), chi harmonic; u = e_y; C = c."""
    V = Vac(K)
    J = V.J
    c, eps = 0.7, 0.3
    Phi = np.zeros(J.N)
    psi = np.zeros(J.N)
    for mon, v in {(0, 1, 0): c, (0, 0, 1): 1.0, (0, 0, 2): 0.5, (2, 0, 0): -0.5,
                   (0, 0, 3): eps, (2, 0, 1): -3 * eps}.items():
        Phi[J.pos[mon]] = v
    for mon, v in {(1, 0, 0): -1.0, (1, 0, 1): -1.0, (1, 0, 2): -3 * eps, (3, 0, 0): eps}.items():
        psi[J.pos[mon]] = v
    x = np.concatenate([Phi[V.mPhi], psi[V.mpsi], [c]])
    print(f"exact symmetric solution, order {K}: residual of the equations "
          f"{np.linalg.norm(V.residual(x, pins=False)):.1e}, sym_B {V.sym_B(x):.1e}, strain {V.strain(x):.1e}")
    y = x.copy()
    y[J.pos[(1, 1, 1)] - 1] += 0.1          # spoil it: the equations must notice
    print(f"same with one coefficient changed: residual {np.linalg.norm(V.residual(y, pins=False)):.1e}")


def run(job):
    K, seed, steps, scale = job
    V = Vac(K)
    rng = np.random.default_rng(1000 * K + seed)
    x0 = rng.normal(size=V.nvar) * scale ** V.var_deg
    t0 = time.time()
    x, nr, it, ok = V.solve(x0, steps)
    if not ok:
        return f"K={K:2d} seed {seed:3d}  no convergence (residual {nr:7.1e} after {it} steps)  [{time.time() - t0:.0f}s]"
    return (f"K={K:2d} seed {seed:3d}  sym_B {V.sym_B(x):8.2e}  strain {V.strain(x):8.2e}  C {x[-1]:+8.3f}"
            f"  max|x| {np.abs(x).max():7.1e}  residual {nr:7.1e}  ({it} steps)  [{time.time() - t0:.0f}s]")


def pad(lower, x, upper):
    """The same polynomial solution as unknowns of a higher order (new coefficients zero)."""
    Phi, psi, C = lower.unpack(x)
    P, Q = np.zeros(upper.J.N), np.zeros(upper.J.N)
    P[:lower.J.N], Q[:lower.J.N] = Phi, psi
    return np.concatenate([P[upper.mPhi], Q[upper.mpsi], [C]])


def ladder(job):
    """Solve at order K0 from a random start, then carry the solution up order by order."""
    K0, seed, steps, scale, kmax = job
    V = Vac(K0)
    rng = np.random.default_rng(1000 * K0 + seed)
    x, nr, it, ok = V.solve(rng.normal(size=V.nvar) * scale ** V.var_deg, steps)
    if not ok:
        return f"seed {seed:3d}  K={K0}: no convergence from the random start"
    out = [f"seed {seed:3d}  K={K0}: sym_B {V.sym_B(x):8.2e} strain {V.strain(x):8.2e} C {x[-1]:+7.3f}"]
    for K in range(K0 + 1, kmax + 1):
        W = Vac(K)
        t0 = time.time()
        x, nr, it, ok = W.solve(pad(V, x, W), steps)
        V = W
        if not ok:
            out.append(f"          K={K}: NOT CONVERGED (residual {nr:7.1e}, {it} steps): sym_B {V.sym_B(x):8.2e} "
                       f"strain {V.strain(x):8.2e} max|x| {np.abs(x).max():7.1e}")
            break
        out.append(f"          K={K}: sym_B {V.sym_B(x):8.2e} strain {V.strain(x):8.2e} C {x[-1]:+7.3f} "
                   f"max|x| {np.abs(x).max():7.1e} ({it} steps) [{time.time() - t0:.0f}s]")
    return "\n".join(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmin", type=int, default=5)
    ap.add_argument("--kmax", type=int, default=10)
    ap.add_argument("--seeds", default="0:8")
    ap.add_argument("--steps", type=int, default=200)
    ap.add_argument("--scale", type=float, default=0.4)
    ap.add_argument("--workers", type=int, default=multiprocessing.cpu_count())
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--ladder", action="store_true",
                    help="solve at --kmin from random starts, then carry each solution up to --kmax")
    a = ap.parse_args()
    if a.check:
        exact_check()
    lo, hi = (int(v) for v in a.seeds.split(":"))
    if a.ladder:
        jobs, fn = [(a.kmin, s, a.steps, a.scale, a.kmax) for s in range(lo, hi)], ladder
    else:
        jobs, fn = [(K, s, a.steps, a.scale) for K in range(a.kmin, a.kmax + 1) for s in range(lo, hi)], run
    with multiprocessing.get_context("fork").Pool(a.workers) as pool:
        for line in pool.imap(fn, jobs):
            print(line, flush=True)
