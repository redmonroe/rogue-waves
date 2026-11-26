"""Strong-form local test over several independent solutions.

For each seed: find non-symmetric solutions at order 4, keep the best scaled, carry it upward
with the step of local_extend.py (capped), and print at each order the symmetry measure of B,
the stretch in u, the largest coefficient and whether the step converged. The question: does
every solution end on the symmetric family, and at which order?

Usage:  python local_strong_seeds.py [first_seed] [n_seeds] [Kmax] [max_steps]
"""
import sys
import time

import numpy as np

from local_extend import extend
from local_jets_count import Variant
from local_jets_project import trust_region

first = int(sys.argv[1]) if len(sys.argv) > 1 else 10
n = int(sys.argv[2]) if len(sys.argv) > 2 else 6
Kmax = int(sys.argv[3]) if len(sys.argv) > 3 else 9
max_steps = int(sys.argv[4]) if len(sys.argv) > 4 else 120
MODE = "strong"

for seed in range(first, first + n):
    K = 4
    J = Variant(K, MODE)
    rng = np.random.default_rng(seed)
    var_deg = np.concatenate([J.deg[m] for m in J.masks])
    best = None
    for _ in range(12):
        x, nr, _ = J.solve(rng, max_iter=60)
        if nr > 1e-10:
            x0 = np.concatenate([rng.normal(size=k) * 0.4 ** J.deg[m] for m, k in zip(J.masks, J.sizes)])
            x, nr, _ = trust_region(J, x0)
        if nr > 1e-10 or J.sym_B(x, 2) < 0.2:
            continue
        g = max(np.linalg.norm(x[var_deg == d]) ** (1.0 / d) for d in range(2, K + 1))
        if best is None or g < best[0]:
            best = (g, x.copy())
    if best is None:
        print(f"seed {seed}: no non-symmetric start found", flush=True)
        continue
    x = best[1]
    np.save(f"runs/local_{MODE}_seed{seed}_K{K}.npy", x)
    print(f"seed {seed}", flush=True)
    while True:
        print(f"   K={K}  sym_B {J.sym_B(x, 2):8.2e}  kill_u {J.kill_u(x):8.2e}  max|x| {np.abs(x).max():7.1e}"
              f"  residual {np.linalg.norm(J.residual(x)):7.1e}", flush=True)
        if K == Kmax:
            break
        t0 = time.time()
        steps = []
        ok = extend(K, seed, max_steps, mode=MODE, out=steps.append)
        Jn = Variant(K + 1, MODE)
        if not ok:
            xp = np.load(f"runs/local_{MODE}_seed{seed}_K{K + 1}_partial.npy")
            print(f"   K={K + 1}  NOT CONVERGED after {len(steps) - 2} steps [{time.time() - t0:.0f}s]: "
                  f"sym_B {Jn.sym_B(xp, 2):8.2e}  kill_u {Jn.kill_u(xp):8.2e}  max|x| {np.abs(xp).max():7.1e}"
                  f"  residual {np.linalg.norm(Jn.residual(xp)):7.1e}", flush=True)
            break
        J, K = Jn, K + 1
        x = np.load(f"runs/local_{MODE}_seed{seed}_K{K}.npy")
