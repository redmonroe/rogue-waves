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
        if self.mode == "mhs":  # switch the quasisymmetry equations off by pinning u and psi to zero
            n_eq = int(self.eq_mask.sum())
            keep = np.ones(len(r), dtype=bool)
            keep[n_eq:2 * n_eq] = False            # div u
            keep[5 * n_eq:9 * n_eq] = False        # B x u = grad psi, u . grad|B|^2
            keep[-3:-1] = False                    # pins on u(0)
            pin = np.concatenate([f for f in u] + [psi])
            return np.concatenate([r[keep], pin])
        raise ValueError(self.mode)


def count(K, mode, seed=0, kmax=4):
    J = Variant(K, mode)
    rng = np.random.default_rng(seed)
    for _ in range(8):
        x, nr, _ = J.solve(rng)
        if nr < 1e-10:
            break
    else:
        return None
    Jac = J.jacobian(x)
    _, s, Vt = np.linalg.svd(Jac)
    rank = int((s > 1e-9 * s[0]).sum())
    ker = Vt[rank:].T  # columns span the tangent space of the solution set
    dims = []
    for k in range(1, kmax + 1):
        rows = np.concatenate([np.where(J.masks[f] & (J.deg <= k))[0] * 0 + J.offsets[f]
                               + np.cumsum(J.masks[f])[J.masks[f] & (J.deg <= k)] - 1 for f in range(3)])
        sv = np.linalg.svd(ker[rows], compute_uv=False)
        dims.append(int((sv > 1e-7).sum()))
    return dims


if __name__ == "__main__":
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    modes = sys.argv[2].split(",") if len(sys.argv) > 2 else ["mhs", "qs", "killing"]
    kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    Kmin = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    ks = ",".join(str(k) for k in range(1, kmax + 1))
    print("free parameters in the Taylor polynomial of B of degree <= k, among degree-K truncated solutions")
    print(f"{'K':>2s}  " + "  ".join(f"{m + '  k=' + ks:>{5 * kmax + 4}s}" for m in modes))
    for K in range(Kmin, Kmax + 1):
        t0 = time.time()
        row = [count(K, m, kmax=kmax) for m in modes]
        cells = [" ".join(f"{v:4d}" for v in r) if r else "no solution found" for r in row]
        print(f"{K:2d}  " + "  ".join(f"{c:>{5 * kmax + 4}s}" for c in cells) + f"   [{time.time() - t0:.0f}s]", flush=True)
