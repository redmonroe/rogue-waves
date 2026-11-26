"""Can a symmetric local solution be bent into a non-symmetric quasisymmetric one?

Take a truncated solution that has a rigid-motion symmetry (mode `killing`), and count the
directions in which it can be changed, to first order, while staying a solution of

    sym : force balance + quasisymmetry + "u is a rigid motion"   (the symmetric neighbours)
    qs  : force balance + quasisymmetry only

Both counts are for the Taylor polynomial of B of degree <= k, by the stable elimination of
local_jets_project.py. If `qs` equals `sym` once the order K is high enough, every first-order
quasisymmetric change of this solution is itself symmetric: nothing non-symmetric branches off
it. If `qs` stays larger, the extra directions are first-order non-symmetric quasisymmetric
perturbations, the place to look for a non-symmetric family (whether they survive at second
order is a separate question).

Unlike the carried non-symmetric solution of local_jets_project.py, a symmetric solution
extends to the next order without trouble, so this reaches high order cheaply.

Usage:  python local_rigidity.py [Kmax] [seed]
"""
import os
import sys
import time

import numpy as np

from local_jets_count import Variant
from local_jets_project import projected_counts, trust_region


def run(Kmax, seed=0, kmax=9, out=print, ckpt_dir="runs"):
    os.makedirs(ckpt_dir, exist_ok=True)
    path = lambda K: os.path.join(ckpt_dir, f"local_killing_seed{seed}_K{K}.npy")
    rng = np.random.default_rng(seed)
    saved = [K for K in range(3, Kmax + 1) if os.path.exists(path(K))]
    if saved:
        K = max(saved)
        S = Variant(K, "killing")
        x = np.load(path(K))
        out(f"resuming from saved order {K}")
    else:
        K = 3
        S = Variant(K, "killing")
        for _ in range(20):
            x, nr, _ = S.solve(rng)
            if nr < 1e-10:
                break
        else:
            out("no symmetric solution found at order 3")
            return
        np.save(path(K), x)
    while True:
        t0 = time.time()
        Q = Variant(K, "qs")
        k = min(kmax, K - 1)
        sym, g1 = projected_counts(S, x, k)
        qs, g2 = projected_counts(Q, x, k)
        out(f"K={K:2d}  sym " + " ".join(f"{v:4d}" for v in sym) + f"   gap {np.log10(g1):4.1f}\n"
            f"      qs  " + " ".join(f"{v:4d}" for v in qs) + f"   gap {np.log10(g2):4.1f}"
            f"   sym_B {S.sym_B(x, 2):7.1e}  residual {np.linalg.norm(S.residual(x)):7.1e}  [{time.time() - t0:.0f}s]")
        if K == Kmax:
            return
        Sn = Variant(K + 1, "killing")
        xn, nr, it = Sn.solve(rng, x0=Sn.extend_from(S, x, rng, noise=0.0))
        if nr > 1e-10:
            xn, nr, it = trust_region(Sn, Sn.extend_from(S, x, rng, noise=0.0))
        if nr > 1e-10:
            out(f"K={K + 1:2d}  could not be extended (residual {nr:.1e})")
            return
        S, x, K = Sn, xn, K + 1
        np.save(path(K), x)


if __name__ == "__main__":
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    run(Kmax, seed, out=lambda s: print(s, flush=True))
