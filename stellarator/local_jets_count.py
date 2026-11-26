"""How much freedom do local quasisymmetric equilibria have, compared with symmetric ones?

For a degree-K truncated solution found by local_jets.py, the tangent space of the solution set
is the kernel of the Jacobian. Its projection onto the Taylor coefficients of B of degree <= k
has some dimension d(k, K). For fixed k this can only fall as K grows, because higher-degree
equations add compatibility conditions; once it stops falling it is the number of free
parameters in the degree-k Taylor polynomial of B among genuine (formal) solutions.

Three systems are compared, with the same normalisation:
  mhs      : equilibrium only (div B = 0, (curl B) x B = grad p)
  qs       : equilibrium + weak quasisymmetry (the system in local_jets.py)
  killing  : qs, plus u required to be a Killing field (rigid motion) to every order
  iso      : qs, plus |B| constant on each flux surface (isodynamic). Such fields are
             quasisymmetric near a point for a trivial reason: any u in the surface will do.
If `qs` and `killing` agree for large K, quasisymmetry is forcing a symmetry locally.
If `qs` stays well above `killing`, non-symmetric local solutions exist, at least formally.

Usage:  python local_jets_count.py [Kmax]
"""
import sys
import time

import numpy as np

from local_jets import Jets


class Variant(Jets):
    def __init__(self, K, mode):
        super().__init__(K)
        self.mode = mode

    def residual(self, x):
        r = super().residual(x)
        if self.mode == "qs":
            return r
        B, u, psi, p = self.unpack(x)
        if self.mode == "killing":
            d = self.d
            extra = [(d(u[j], i) + d(u[i], j))[self.eq_mask] for i in range(3) for j in range(i, 3)]
            return np.concatenate([r] + extra)
        if self.mode == "iso":  # field strength constant on flux surfaces: grad|B|^2 x grad psi = 0
            d, m = self.d, self.mul
            B2 = sum(m(B[i], B[i]) for i in range(3))
            g = self.cross([d(B2, i) for i in range(3)], [d(psi, i) for i in range(3)])
            # the qs normalisation d_z B_3(0) = 1 says |B| varies along the field line, which
            # excludes this case; normalise the variation across the surface instead
            r = r.copy()
            r[-1] = B[2][self.lin[1]] - 1
            return np.concatenate([r] + [c[self.eq_mask] for c in g])
        if self.mode == "mhs":  # switch the quasisymmetry equations off by pinning u and psi to zero
            n_eq = int(self.eq_mask.sum())
            keep = np.ones(len(r), dtype=bool)
            keep[n_eq:2 * n_eq] = False            # div u
            keep[5 * n_eq:9 * n_eq] = False        # B x u = grad psi, u . grad|B|^2
            keep[-3:-1] = False                    # pins on u(0)
            pin = np.concatenate([f for f in u] + [psi])
            return np.concatenate([r[keep], pin])
        raise ValueError(self.mode)


def gap_rank(s, floor=1e-13):
    """Numerical rank from the largest gap in the singular values. Returns (rank, gap), where
    gap = (smallest value kept) / (largest value dropped). A gap below about 1e3 means the count
    is not trustworthy at this order."""
    if len(s) < 2:
        return len(s), np.inf
    s = np.maximum(s, floor * s[0])
    ratios = s[:-1] / s[1:]
    i = int(np.argmax(ratios))
    if ratios[i] < 10:  # no gap at all: full rank
        return len(s), 1.0
    return i + 1, float(ratios[i])


def tangent_dims(J, x, kmax, rho=0.4):
    """Dimension of the projection of the solution set's tangent space onto the Taylor
    coefficients of B of degree <= k, for k = 1 .. kmax, with the reliability gap of each count.

    The Jacobian is rescaled as if lengths were measured in units of 1/rho (coefficients of
    degree d scale by rho^d). Without this, free directions that need large high-degree
    coefficients have singular values that shrink tenfold per order and get miscounted."""
    Jac = J.jacobian(x)
    col = np.concatenate([rho ** J.deg[m] for m in J.masks])
    n_eq = int(J.eq_mask.sum())
    blocks = (len(Jac) - 6) // n_eq if J.mode != "mhs" else None
    if blocks is None:
        row = np.ones(len(Jac))
    else:
        row = np.concatenate([rho ** (-(J.deg[J.eq_mask] + 1.0))] * 9 + [np.ones(6)]
                             + [rho ** (-(J.deg[J.eq_mask] + 1.0))] * (blocks - 9))
    _, s, Vt = np.linalg.svd(Jac * col * row[:, None])
    if rho == 1.0:  # plain thresholds; adequate for the symmetric reference, whose counts settle at once
        rank, gap0 = int((s > 1e-9 * s[0]).sum()), np.inf
    else:
        rank, gap0 = gap_rank(s)
        if rank == len(s) and Jac.shape[1] > len(s):
            gap0 = np.inf
    ker = Vt[rank:].T
    dims, gaps = [], [gap0]
    for k in range(1, kmax + 1):
        rows = np.concatenate([J.offsets[f] + np.cumsum(J.masks[f])[J.masks[f] & (J.deg <= k)] - 1 for f in range(3)])
        sv = np.linalg.svd(ker[rows], compute_uv=False)
        r, g = (int((sv > 1e-7).sum()), np.inf) if rho == 1.0 else gap_rank(sv)
        dims.append(r)
        gaps.append(g)
    return dims, J.nvar - rank, gaps


def ladder(mode, Kmax, kmax, seed=0, Kstart=3, out=print):
    """Find one solution at order Kstart, then carry it up one order at a time."""
    rng = np.random.default_rng(seed)
    J = Variant(Kstart, mode)
    for _ in range(20):
        x, nr, _ = J.solve(rng)
        if nr < 1e-10:
            break
    else:
        out(f"{mode}: no solution found at order {Kstart}")
        return
    K = Kstart
    while True:
        t0 = time.time()
        dims, free, gaps = tangent_dims(J, x, min(kmax, K), rho=0.4 if mode == "qs" else 1.0)
        extra = ""
        if mode != "mhs":
            extra = f"  sym_B(deg<=2) {J.sym_B(x, 2):8.2e}  kill_u {J.kill_u(x):8.2e}"
        out(f"{mode:8s} K={K:2d}  free dims {free:5d}  B-polynomial freedom k=1.. : "
            + " ".join(f"{v:4d}" for v in dims) + extra
            + "  log10 gaps: " + " ".join("inf" if not np.isfinite(g) else f"{np.log10(g):.0f}" for g in gaps)
            + f"   [{time.time() - t0:.0f}s]")
        if K == Kmax:
            return
        Jn = Variant(K + 1, mode)
        t0 = time.time()
        xn, nr, it = Jn.solve(rng, x0=Jn.extend_from(J, x, rng, noise=0.0))
        if nr > 1e-10:  # retry with a little randomness in the new coefficients
            xn, nr, it = Jn.solve(rng, x0=Jn.extend_from(J, x, rng, noise=1e-2))
        if nr > 1e-10:
            out(f"{mode:8s} K={K + 1:2d}  could not be extended (residual {nr:.1e} after {it} iterations, {time.time() - t0:.0f}s)")
            return
        J, x, K = Jn, xn, K + 1


if __name__ == "__main__":
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    modes = sys.argv[2].split(",") if len(sys.argv) > 2 else ["qs", "killing", "mhs"]
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    say = lambda s: print(s, flush=True)
    for mode in modes:
        ladder(mode, Kmax, kmax, seed=seed, out=say)
