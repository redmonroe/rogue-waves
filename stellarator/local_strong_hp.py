"""Strong-form local test with extended-precision checking.

The double-precision runs (local_strong_seeds.py) end with coefficients of size 1e3 to 1e4 and a
symmetry measure of 1e-14, which is at the limit of what double precision can say. This script
repeats them and, at every order, polishes the solution in extended precision:

    residual evaluated with mpmath at `--digits` significant digits (default 40);
    Newton correction computed from the double-precision Jacobian (mixed-precision
    iterative refinement: each pass gains several digits, and the limit is an exact solution
    of the extended-precision equations, not of the rounded ones).

Then the symmetry measure of B and the stretch in u are computed in extended precision too. So
"symmetric to 1e-14" becomes either "symmetric to 1e-35" or a definite non-zero number.

Two modes:
  ladder (default): for each seed, best-scaled non-symmetric solution at order 4, carried upward.
  --scratch K     : for each seed, solve at order K directly from a random start (no ladder).

Needs only mpmath (already in the project venv). `pip install gmpy2` makes it several times
faster; mpmath picks it up by itself.

Usage:  python local_strong_hp.py [--seeds 100:108] [--kmax 9] [--digits 40] [--steps 150]
                                  [--scratch K] [--mode strong] [--rescale] [--out FILE]
Output: one line per seed and order, also appended to runs/hp_<mode>.txt.
"""
import argparse
import os
import time

import mpmath as mp
import numpy as np

from local_extend import extend, jacobian
from local_jets_count import Variant
from local_jets_project import trust_region


class HP(Variant):
    """Same equations as Variant, on object arrays of mpmath numbers. Only the series product
    needs replacing; everything else in Jets / Variant is written with generic arithmetic."""

    def mul(self, a, b):
        out = np.zeros(self.N, dtype=object)
        np.add.at(out, self.O, a[self.I] * b[self.J])
        return out


def to_hp(x):
    return np.array([mp.mpf(float(v)) for v in x], dtype=object)


def to_float(x):
    return np.array([float(v) for v in x])


def norm(v):
    return mp.sqrt(mp.fsum(t * t for t in v))


def hp_sym_B(J, x, dmax=2):
    """Jets.sym_B in extended precision: smallest singular value of (a, w) -> L_K B to degree dmax."""
    B, _, _, _ = J.unpack(x)
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
    M = mp.matrix(len(cols[0]), 6)
    for j in range(6):
        for i in range(len(cols[0])):
            M[i, j] = mp.mpf(cols[j][i])
    return mp.svd_r(M, compute_uv=False)[5]


def hp_kill_u(J, x):
    """Jets.kill_u in extended precision."""
    _, u, _, _ = J.unpack(x)
    G = [[u[j][J.lin[i]] for j in range(3)] for i in range(3)]
    S = [[G[i][j] + G[j][i] for j in range(3)] for i in range(3)]
    c = u[2][J.i0]
    M = [[0, 1, 0], [1, 0, c], [0, c, 0]]
    num = sum(S[i][j] * M[i][j] for i in range(3) for j in range(3))
    den = sum(M[i][j] * M[i][j] for i in range(3) for j in range(3))
    kappa = num / den
    return mp.sqrt(sum((S[i][j] - kappa * M[i][j]) ** 2 for i in range(3) for j in range(3)))


def polish(Jd, Jh, x, digits, workers, max_pass=60):
    """Mixed-precision refinement. Returns (x, residual, passes, status)."""
    r = Jh.residual(x)
    nr = norm(r)
    svd = None
    for it in range(max_pass):
        size = max(1.0, max(abs(float(v)) for v in x))
        if nr < mp.mpf(10) ** (-digits) * size:
            return x, nr, it, "ok"
        if svd is None:
            svd = np.linalg.svd(jacobian(Jd, to_float(x), workers), full_matrices=False)
            fresh = True
        U, s, Vt = svd
        proj = U.T @ to_float(r)
        done = False
        for cut in (1e-12, 1e-10, 1e-13, 1e-8, 1e-14):
            keep = s > cut * s[0]
            xn = x - Vt[keep].T @ (proj[keep] / s[keep])
            rn = Jh.residual(xn)
            nn = norm(rn)
            if nn < 0.5 * nr:
                x, r, nr, done, fresh = xn, rn, nn, True, False
                break
        if not done:
            if fresh:
                return x, nr, it, "stalled"
            svd = None  # the Jacobian has gone stale: rebuild once before giving up
    return x, nr, max_pass, "pass limit"


def report(tag, K, Jh, x, nr, passes, status, t0, log):
    line = (f"{tag}  K={K}  sym_B {mp.nstr(hp_sym_B(Jh, x), 3):>9s}  kill_u {mp.nstr(hp_kill_u(Jh, x), 3):>9s}"
            f"  max|x| {max(abs(float(v)) for v in x):7.1e}  residual {mp.nstr(nr, 3):>9s}"
            f"  ({passes} passes, {status})  [{time.time() - t0:.0f}s]")
    print(line, flush=True)
    log.write(line + "\n")
    log.flush()


def random_start(J, rng, scale):
    return np.concatenate([rng.normal(size=n) * scale ** J.deg[m] for m, n in zip(J.masks, J.sizes)])


def start_order4(seed, mode):
    """Best-scaled non-symmetric double-precision solution at order 4 (as local_strong_seeds.py)."""
    J = Variant(4, mode)
    rng = np.random.default_rng(seed)
    var_deg = np.concatenate([J.deg[m] for m in J.masks])
    best = None
    for _ in range(12):
        x, nr, _ = J.solve(rng, max_iter=60)
        if nr > 1e-10:
            x, nr, _ = trust_region(J, random_start(J, rng, 0.4))
        if nr > 1e-10 or J.sym_B(x, 2) < 0.2:
            continue
        g = max(np.linalg.norm(x[var_deg == d]) ** (1.0 / d) for d in range(2, 5))
        if best is None or g < best[0]:
            best = (g, x.copy())
    return None if best is None else best[1]


def ladder(seed, a, log, workers):
    tag = f"seed {seed}"
    x = start_order4(seed, a.mode)
    if x is None:
        print(f"{tag}: no non-symmetric start found", flush=True)
        return
    K = 4
    path = lambda k, t="": f"runs/local_{a.mode}_seed{seed}_K{k}{t}.npy"
    for k in range(4, a.kmax + 2):                       # a fresh run, not a resume
        for t in ("", "_partial"):
            if os.path.exists(path(k, t)):
                os.remove(path(k, t))
    while True:
        t0 = time.time()
        Jd, Jh = Variant(K, a.mode), HP(K, a.mode)
        xh, nr, passes, status = polish(Jd, Jh, to_hp(x), a.digits, workers)
        report(tag, K, Jh, xh, nr, passes, status, t0, log)
        if K == a.kmax:
            return
        np.save(path(K), to_float(xh))
        steps = []
        ok = extend(K, seed, a.steps, mode=a.mode, out=steps.append, workers=workers, rescale=a.rescale)
        K += 1
        if not ok:
            x = np.load(path(K, "_partial"))
            Jd = Variant(K, a.mode)
            line = (f"{tag}  K={K}  double-precision step NOT CONVERGED after {len(steps) - 2} steps: "
                    f"sym_B {Jd.sym_B(x, 2):8.2e}  kill_u {Jd.kill_u(x):8.2e}  max|x| {np.abs(x).max():7.1e}"
                    f"  residual {np.linalg.norm(Jd.residual(x)):7.1e}")
            print(line, flush=True)
            log.write(line + "\n")
            return
        x = np.load(path(K))


def scratch(seed, a, log, workers):
    """Solve at order a.scratch from a random start, no ladder."""
    K = a.scratch
    tag = f"scratch seed {seed}"
    Jd, Jh = Variant(K, a.mode), HP(K, a.mode)
    rng = np.random.default_rng(seed)
    t0 = time.time()
    np.save(f"runs/local_{a.mode}_seed{seed}_K{K}_partial.npy", random_start(Jd, rng, 0.4))
    # reuse the logged Newton step: it starts from the "partial" file of order K
    lower = f"runs/local_{a.mode}_seed{seed}_K{K - 1}.npy"
    made = not os.path.exists(lower)
    if made:
        np.save(lower, np.zeros(Variant(K - 1, a.mode).nvar))
    steps = []
    ok = extend(K - 1, seed, a.steps, mode=a.mode, out=steps.append, workers=workers)  # random start: no rescaling
    if made:
        os.remove(lower)
    if not ok:
        line = f"{tag}  K={K}  no convergence from the random start ({steps[-2].strip()[:70]})"
        print(line, flush=True)
        log.write(line + "\n")
        return
    x = np.load(f"runs/local_{a.mode}_seed{seed}_K{K}.npy")
    xh, nr, passes, status = polish(Jd, Jh, to_hp(x), a.digits, workers)
    report(tag, K, Jh, xh, nr, passes, status, t0, log)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="100:104")
    ap.add_argument("--kmax", type=int, default=9)
    ap.add_argument("--digits", type=int, default=40)
    ap.add_argument("--steps", type=int, default=150)
    ap.add_argument("--scratch", type=int, default=0)
    ap.add_argument("--mode", default="strong")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--rescale", action="store_true",
                    help="rescale the length unit in each step (for solutions with large coefficients)")
    ap.add_argument("--out", default="", help="log file (default runs/hp_<mode>.txt)")
    a = ap.parse_args()
    mp.mp.dps = a.digits + 15
    lo, hi = (int(v) for v in a.seeds.split(":"))
    os.makedirs("runs", exist_ok=True)
    with open(a.out or f"runs/hp_{a.mode}.txt", "a") as log:
        log.write(f"# {time.strftime('%Y-%m-%d %H:%M')}  seeds {a.seeds}  kmax {a.kmax}  digits {a.digits}"
                  f"  scratch {a.scratch}  rescale {a.rescale}  backend {mp.libmp.BACKEND}\n")
        for seed in range(lo, hi):
            (scratch if a.scratch else ladder)(seed, a, log, a.workers)
