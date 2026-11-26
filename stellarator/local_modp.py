"""First-order rigidity at an explicit symmetric equilibrium, in exact arithmetic.

Question: take an axisymmetric equilibrium (which is quasisymmetric, with u the rotation). To
first order, in how many ways can the Taylor polynomial of B of degree <= k be changed so that
the field stays a solution of

    sym    : force balance + quasisymmetry with u a rigid motion      (stays symmetric)
    strong : force balance + strong quasisymmetry
    weak   : force balance + weak quasisymmetry

If `strong` equals `sym`, no non-symmetric strongly quasisymmetric family can branch off this
equilibrium, to first order.

Everything is done with integers modulo the prime p = 2^31 - 1: the equilibrium's Taylor
coefficients, the residuals (which must vanish mod p: checked), the Jacobian, and the ranks by
Gaussian elimination. There is no rounding and no singular-value cutoff, which is what limited
local_jets_project.py and local_rigidity.py. A count computed mod p equals the count over the
rationals except for an unlucky prime or unlucky parameters; `--prime` and `--seed` let the run
be repeated to rule that out. (It is exact arithmetic on a random instance, not yet a proof.)

The equilibrium is a Solov'ev solution: psi = A R^4/8 + C0 z^2/2 + h1 R^2 + h2 (R^4 - 4 R^2 z^2)
+ h3 R^2 z + h4 z,  B = grad psi x grad phi + F grad phi,  F^2 = F0^2 + 2 s C0 psi,
p = s' A psi, u = (y, -x, 0), expanded about a random point with random parameters. No
normalisation is imposed, so the counts include moving, turning and rescaling the field; these
are symmetric changes and appear in all three rows alike.

Special cases of the equilibrium, for contrast: --case vacuum (no pressure, no current),
force-free (no pressure), F-constant (pressure, but toroidal field function constant).

Usage:  python local_modp.py [--kmin 5] [--kmax 10] [--seed 0] [--prime 2147483647] [--case general]
"""
import argparse
import time

import numpy as np

from local_jets_count import Variant


class Mod(Variant):
    """Variant's equations on int64 arrays modulo P."""

    def __init__(self, K, mode, P):
        super().__init__(K, mode)
        self.P = P
        order = np.argsort(self.O, kind="stable")
        self.Is, self.Js = self.I[order], self.J[order]
        Os = self.O[order]
        self.starts = np.flatnonzero(np.r_[True, Os[1:] != Os[:-1]])
        self.outs = Os[self.starts]
        self.dint = [(frm, to, fac.astype(np.int64)) for frm, to, fac in self.dtab]
        base = 9 * int(self.eq_mask.sum())
        self.pin_rows = np.arange(base, base + 6)

    def mul(self, a, b):
        P = self.P
        w = (a % P)[self.Is] * (b % P)[self.Js] % P
        out = np.zeros(self.N, dtype=np.int64)
        out[self.outs] = np.add.reduceat(w, self.starts) % P
        return out

    def d(self, a, ax):
        frm, to, fac = self.dint[ax]
        out = np.zeros(self.N, dtype=np.int64)
        out[to] = (a % self.P)[frm] * fac % self.P
        return out

    def res(self, x):
        """Residual without the six normalisation rows."""
        return np.delete(self.residual(x), self.pin_rows) % self.P

    # -- series helpers
    def const(self, c):
        a = np.zeros(self.N, dtype=np.int64)
        a[self.i0] = c % self.P
        return a

    def inv(self, a):
        P = self.P
        r = self.const(pow(int(a[self.i0]), P - 2, P))
        for _ in range(6):
            r = self.mul(r, (2 * self.const(1) - self.mul(a, r)) % P)
        return r

    def sqrt(self, a, root0):
        """Square root of a series whose constant term is root0^2."""
        P = self.P
        s = self.const(root0)
        half = pow(2, P - 2, P)
        for _ in range(6):
            s = (s + self.mul(a, self.inv(s))) % P * half % P
        return s


def solovev(J, rng, sign_F, sign_p, case="general"):
    """Taylor series mod P of the Solov'ev equilibrium about a random point."""
    P = J.P
    m = J.mul
    r = lambda: int(rng.integers(2, P - 1))
    x0, A, C0, h1, h2, h3, h4, f0 = [r() for _ in range(3)], r(), r(), r(), r(), r(), r(), r()
    if case in ("vacuum", "force-free"):   # no pressure gradient
        A = 0
    if case in ("vacuum", "F-constant"):   # toroidal field function F constant
        C0 = 0
    X = []
    for i in range(3):
        a = J.const(x0[i])
        a[J.lin[i]] = 1
        X.append(a)
    x, y, z = X
    R2 = (m(x, x) + m(y, y)) % P
    R4 = m(R2, R2)
    z2 = m(z, z)
    inv8, inv2 = pow(8, P - 2, P), pow(2, P - 2, P)
    psi = (A * inv8 % P * R4 + C0 * inv2 % P * z2 + h1 * R2 + h2 * ((R4 - 4 * m(R2, z2)) % P)
           + h3 * m(R2, z) + h4 * z) % P
    iR2 = J.inv(R2)
    gphi = [(-m(y, iR2)) % P, m(x, iR2), np.zeros(J.N, dtype=np.int64)]
    gpsi = [J.d(psi, i) for i in range(3)]
    # F^2 = f0^2 + 2 sign_F C0 (psi - psi(x0)), so that F(x0) = f0
    F2 = (J.const(f0 * f0) + 2 * sign_F * C0 * ((psi - J.const(psi[J.i0])) % P)) % P
    F = J.sqrt(F2, f0)
    c = J.cross(gpsi, gphi)
    B = [(c[i] + m(F, gphi[i])) % P for i in range(3)]
    u = [y % P, (-x) % P, np.zeros(J.N, dtype=np.int64)]
    p = sign_p * A * psi % P
    fields = B + u + [psi, p]
    return np.concatenate([f[mk] for f, mk in zip(fields, J.masks)]) % P


def jacobian(J, x):
    """Exact Jacobian mod P. The residual has degree <= 3 in the unknowns, so
    R'(x) v = [8 (R(x+v) - R(x-v)) - (R(x+2v) - R(x-2v))] / 12."""
    P = J.P
    inv12 = pow(12, P - 2, P)
    cols = np.empty((len(J.res(x)), J.nvar), dtype=np.int64)
    xv = x.copy()
    for j in range(J.nvar):
        vals = []
        for t in (1, -1, 2, -2):
            xv[j] = (x[j] + t) % P
            vals.append(J.res(xv))
        xv[j] = x[j]
        cols[:, j] = (8 * (vals[0] - vals[1]) - (vals[2] - vals[3])) % P * inv12 % P
    return cols


def pivots(M, P):
    """Column-ordered Gaussian elimination mod P. Returns the list of pivot columns."""
    M = M.copy() % P
    rows, cols = M.shape
    piv, r = [], 0
    for c in range(cols):
        if r == rows:
            break
        nz = np.flatnonzero(M[r:, c])
        if nz.size == 0:
            continue
        i = r + nz[0]
        if i != r:
            M[[r, i]] = M[[i, r]]
        M[r, c:] = M[r, c:] * pow(int(M[r, c]), P - 2, P) % P
        f = M[r + 1:, c].copy()
        nzr = np.flatnonzero(f)
        if nzr.size:
            M[r + 1 + nzr, c:] = (M[r + 1 + nzr, c:] - f[nzr, None] * M[r, c:][None, :]) % P
        piv.append(c)
        r += 1
    return piv


def counts(J, x, kmax):
    """Free parameters, to first order, in the Taylor polynomial of B of degree <= k."""
    Jac = jacobian(J, x)
    var_deg = np.concatenate([J.deg[m] for m in J.masks])
    is_B = np.zeros(J.nvar, dtype=bool)
    is_B[: J.offsets[3]] = True
    # eliminate everything else first, then B from high degree to low: a B coefficient of
    # degree <= k is free exactly when it is not a pivot
    key = np.where(is_B, -var_deg, -10 ** 6)
    order = np.argsort(key, kind="stable")
    piv = set(pivots(Jac[:, order], J.P))
    is_piv = np.zeros(J.nvar, dtype=bool)
    is_piv[order[sorted(piv)]] = True
    return [int((is_B & (var_deg <= k) & ~is_piv).sum()) for k in range(1, kmax + 1)]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmin", type=int, default=5)
    ap.add_argument("--kmax", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--prime", type=int, default=2147483647)
    ap.add_argument("--case", default="general", choices=["general", "vacuum", "force-free", "F-constant"])
    a = ap.parse_args()
    P = a.prime
    # fix the two sign conventions by which one makes the residual vanish
    signs = None
    for sF in (1, -1):
        for sp in (1, -1):
            J = Mod(4, "killing", P)
            if not np.any(J.res(solovev(J, np.random.default_rng(a.seed), sF, sp, "general"))):
                signs = (sF, sp)
    print(f"case {a.case}, prime {P}, seed {a.seed}; sign conventions {signs} (residual vanishes mod p)", flush=True)
    for K in range(a.kmin, a.kmax + 1):
        t0 = time.time()
        rows = {}
        for name, mode in (("sym", "killing"), ("strong", "strong"), ("weak", "qs")):
            J = Mod(K, mode, P)
            x = solovev(J, np.random.default_rng(a.seed), *signs, a.case)
            assert not np.any(J.res(x)), f"{mode}: the equilibrium is not a solution mod p"
            rows[name] = counts(J, x, min(K - 1, 8))
        print(f"K={K}   [{time.time() - t0:.0f}s]", flush=True)
        for name, v in rows.items():
            print(f"   {name:6s} " + " ".join(f"{c:4d}" for c in v), flush=True)
