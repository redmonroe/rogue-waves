"""Carry a saved truncated solution up one order, with a log line per Newton step.

local_jets_project.py does the same step silently, and at order 14 it ran out of time with
nothing to show. This script prints the residual, the damping and the symmetry measure after
every step and saves the iterate, so a stalled run still says how it stalled and can resume.

The Jacobian is built on all cores. Each step is the minimum-norm Gauss-Newton step with
backtracking, as in local_jets.py; if no length works, the smallest singular values are dropped.
(Levenberg-Marquardt damping, and picking the greediest of several truncated steps, were both
tried and both crawled along a narrow valley where the plain step converges in about 6.)

Usage:  python local_extend.py K [seed] [max_steps]      extends runs/local_qs_seed{seed}_K{K}.npy
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

from local_jets_count import Variant

_J = None
_x = None


def _cols(idx):
    xc = _x.astype(complex)
    out = np.empty((len(idx), _J.neq))
    for n, j in enumerate(idx):
        xc[j] += 1e-30j
        out[n] = _J.residual(xc).imag / 1e-30
        xc[j] = _x[j]
    return out


def jacobian(J, x, workers):
    global _J, _x
    _J, _x = J, x
    chunks = np.array_split(np.arange(J.nvar), workers * 8)
    with Pool(workers) as pool:  # fork: the workers see _J and _x
        parts = pool.map(_cols, chunks)
    return np.vstack(parts).T


def extend(K, seed=0, max_steps=200, workers=os.cpu_count(), out=print, ckpt_dir="runs"):
    path = lambda k, tag="": os.path.join(ckpt_dir, f"local_qs_seed{seed}_K{k}{tag}.npy")
    lower = Variant(K, "qs")
    x_lower = np.load(path(K))
    J = Variant(K + 1, "qs")
    rng = np.random.default_rng(seed)
    if os.path.exists(path(K + 1, "_partial")):
        x = np.load(path(K + 1, "_partial"))
        out(f"resuming the partial order-{K + 1} iterate")
    else:
        x = J.extend_from(lower, x_lower, rng, noise=0.0)
    J.neq = len(J.residual(x))
    r = J.residual(x)
    nr = np.linalg.norm(r)
    out(f"order {K} -> {K + 1}: {J.nvar} unknowns, {J.neq} equations, starting residual {nr:.2e}")
    for it in range(1, max_steps + 1):
        t0 = time.time()
        Jac = jacobian(J, x, workers)
        U, s, Vt = np.linalg.svd(Jac, full_matrices=False)
        Utr = U.T @ r
        # full Gauss-Newton step first, shortened until it gives a real decrease; only if that
        # fails, drop the smallest singular values and try again
        found = False
        for cut in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4):
            keep = s > cut * s[0]
            full = Vt[keep].T @ (Utr[keep] / s[keep])
            t = 1.0
            while t > 1e-3:
                xn = x - t * full
                rn = J.residual(xn)
                nn = np.linalg.norm(rn)
                if nn < (1 - 0.3 * t) * nr:
                    found = True
                    break
                t /= 2
            if found:
                break
        if not found:
            out(f"step {it:3d}  no step reduces the residual ({nr:.3e}); stopping")
            return False
        lam = (cut, t)
        moved = np.linalg.norm(xn[: J.offsets[3]] - x[: J.offsets[3]])
        x, r = xn, rn
        nr = nn
        np.save(path(K + 1, "_partial"), x)
        out(f"step {it:3d}  residual {nr:9.3e}  cutoff {lam[0]:5.0e} length {lam[1]:4.2f}  |dB| {moved:8.2e}  "
            f"sym_B {J.sym_B(x, 2):8.2e}  kill_u {J.kill_u(x):8.2e}  [{time.time() - t0:.0f}s]")
        if nr < 1e-12:
            np.save(path(K + 1), x)
            os.remove(path(K + 1, "_partial"))
            out(f"converged: saved order {K + 1}")
            return True
    out("did not converge")
    return False


if __name__ == "__main__":
    K = int(sys.argv[1])
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    extend(K, seed, steps, out=lambda s: print(s, flush=True))
