"""Carry a saved truncated solution up one order, with a log line per Newton step.

local_jets_project.py does the same step silently, and at order 14 it ran out of time with
nothing to show. This script prints the residual, the damping and the symmetry measure after
every step and saves the iterate, so a stalled run still says how it stalled and can resume.

The Jacobian is built on all cores. Each step is the minimum-norm Gauss-Newton step with
backtracking, as in local_jets.py; if no length works, the smallest singular values are dropped.
(Levenberg-Marquardt damping, and picking the greediest of several truncated steps, were both
tried and both crawled along a narrow valley where the plain step converges in about 6.)

Usage:  python local_extend.py K [seed] [max_steps] [mode] [rescale]   extends runs/local_{mode}_seed{seed}_K{K}.npy
"""
import os
import sys
import time
import multiprocessing

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
    if not hasattr(J, "neq"):
        J.neq = len(J.residual(x))
    _J, _x = J, x
    if workers <= 1:
        return _cols(np.arange(J.nvar)).T
    chunks = np.array_split(np.arange(J.nvar), workers * 8)
    # "fork" explicitly: the workers must inherit _J and _x, and macOS defaults to "spawn"
    with multiprocessing.get_context("fork").Pool(workers) as pool:
        parts = pool.map(_cols, chunks)
    return np.vstack(parts).T


def extend(K, seed=0, max_steps=200, workers=os.cpu_count(), out=print, ckpt_dir="runs", mode="qs",
           rescale=False):
    """rescale=True measures lengths in a unit in which the lower-order solution's coefficients
    are of order one (degree-d coefficients are multiplied by rho^d, degree-e equations by
    rho^e). The solution set is unchanged; what changes is which directions the double-precision
    least-squares step can see. Needed when coefficients reach 1e3 and more."""
    path = lambda k, tag="": os.path.join(ckpt_dir, f"local_{mode}_seed{seed}_K{k}{tag}.npy")
    lower = Variant(K, mode)
    x_lower = np.load(path(K))
    J = Variant(K + 1, mode)
    rng = np.random.default_rng(seed)
    if os.path.exists(path(K + 1, "_partial")):
        x = np.load(path(K + 1, "_partial"))
        out(f"resuming the partial order-{K + 1} iterate")
    else:
        x = J.extend_from(lower, x_lower, rng, noise=0.0)
    J.neq = len(J.residual(x))
    rho = 1.0
    if rescale:
        low_deg = np.concatenate([lower.deg[m] for m in lower.masks])
        growth = max(np.abs(x_lower[low_deg == d]).max() ** (1.0 / d) for d in range(2, K + 1))
        rho = 1.0 / max(1.0, growth)
    from local_jets_project import levels
    var_deg, row_top = levels(J)
    cs = rho ** (-var_deg.astype(float))        # dx = cs * dy
    rs = rho ** (row_top - 1.0)                 # scaled residual = rs * residual
    res = lambda v: rs * J.residual(v)
    r = res(x)
    nr = np.linalg.norm(r)
    out(f"order {K} -> {K + 1}: {J.nvar} unknowns, {J.neq} equations, length unit {rho:.3f}, "
        f"starting residual {nr:.2e}")
    for it in range(1, max_steps + 1):
        t0 = time.time()
        Jac = rs[:, None] * jacobian(J, x, workers) * cs[None, :]
        U, s, Vt = np.linalg.svd(Jac, full_matrices=False)
        Utr = U.T @ r
        # full Gauss-Newton step first, shortened until it gives a real decrease; only if that
        # fails, drop the smallest singular values and try again
        found = False
        for cut in (1e-12, 1e-10, 1e-8, 1e-6, 1e-4):
            keep = s > cut * s[0]
            full = cs * (Vt[keep].T @ (Utr[keep] / s[keep]))
            t = 1.0
            while t > 1e-3:
                xn = x - t * full
                rn = res(xn)
                nn = np.linalg.norm(rn)
                if nn < (1 - 0.3 * t) * nr:
                    found = True
                    break
                t /= 2
            if found:
                break
        if not found:
            if nr < 1e-11 * max(1.0, np.abs(x / cs).max()):   # stalled at rounding level: converged
                np.save(path(K + 1), x)
                os.remove(path(K + 1, "_partial"))
                out(f"step {it:3d}  residual {nr:9.3e} at rounding level; saved order {K + 1}")
                out(f"converged: saved order {K + 1}")
                return True
            out(f"step {it:3d}  no step reduces the residual ({nr:.3e}); stopping")
            return False
        lam = (cut, t)
        moved = np.linalg.norm(xn[: J.offsets[3]] - x[: J.offsets[3]])
        x, r = xn, rn
        nr = nn
        np.save(path(K + 1, "_partial"), x)
        out(f"step {it:3d}  residual {nr:9.3e}  cutoff {lam[0]:5.0e} length {lam[1]:4.2f}  |dB| {moved:8.2e}  "
            f"sym_B {J.sym_B(x, 2):8.2e}  kill_u {J.kill_u(x):8.2e}  [{time.time() - t0:.0f}s]")
        if nr < 1e-12 * max(1.0, np.abs(x / cs).max()):   # relative to the coefficient size
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
    mode = sys.argv[4] if len(sys.argv) > 4 else "qs"
    extend(K, seed, steps, out=lambda s: print(s, flush=True), mode=mode, rescale="rescale" in sys.argv)
