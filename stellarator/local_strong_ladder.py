"""Carry one strong-form solution up the orders with the logged step of local_extend.py,
printing the freedom count, the symmetry measures and the coefficient size at each order.

Start from a saved solution runs/local_strong_seed{seed}_K{K0}.npy (seed 1 was picked at order 4
as the best-scaled of 39 random non-symmetric solutions).

Usage:  python local_strong_ladder.py K0 Kmax [seed] [mode]
"""
import os
import sys

import numpy as np

from local_extend import extend
from local_jets_count import Variant
from local_jets_project import projected_counts

K0, Kmax = int(sys.argv[1]), int(sys.argv[2])
seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
mode = sys.argv[4] if len(sys.argv) > 4 else "strong"
for K in range(K0, Kmax + 1):
    path = f"runs/local_{mode}_seed{seed}_K{K}.npy"
    if not os.path.exists(path):
        steps = []
        ok = extend(K - 1, seed, 300, mode=mode, out=steps.append)
        print(f"   {K - 1} -> {K}: {len(steps) - 1} steps; last: {steps[-2] if ok else steps[-1]}", flush=True)
        if not ok:
            print("   " + steps[-2], flush=True)
            break
    J = Variant(K, mode)
    x = np.load(path)
    dims, gap = projected_counts(J, x, min(9, K - 1))
    print(f"{mode} K={K:2d}  freedom in B, k=1.. : " + " ".join(f"{v:4d}" for v in dims)
          + f"   sym_B {J.sym_B(x, 2):8.2e}  kill_u {J.kill_u(x):8.2e}  max|x| {np.abs(x).max():7.1e}"
          + f"  residual {np.linalg.norm(J.residual(x)):7.1e}  min log10 gap {np.log10(gap):4.1f}", flush=True)
