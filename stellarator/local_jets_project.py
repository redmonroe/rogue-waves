"""Stable version of the freedom count in local_jets_count.py.

Question, as before: among degree-K truncated solutions of force balance + quasisymmetry, how
many free parameters does the Taylor polynomial of B of degree <= k have?

local_jets_count.py answers by computing the kernel of the whole Jacobian and projecting it,
which fails in double precision beyond K = 10: kernel vectors must carry the huge high-degree
responses to a low-degree change. This script never forms them. It eliminates the unknowns one
degree at a time from the top down:

    the equations of degree L-1, together with the constraints carried down so far, are linear in
    the coefficients of degree L; the combinations that do not involve those coefficients are
    constraints on the lower degrees; they are re-orthonormalised and carried down.

Every rank decision is on a matrix of moderate size with entries of order one, and the smallest
gap between kept and dropped singular values is reported (`min log10 gap`). A gap below about 3
means the counts at that order should not be trusted.

Usage:  python local_jets_project.py [Kmax] [modes] [seed]     modes: qs,killing (default qs)
"""
import sys
import time

import numpy as np

from local_jets_count import Variant, gap_rank


def levels(J):
    """Degree of each unknown, and for each equation the highest degree of unknown it involves."""
    var_level = np.concatenate([J.deg[m] for m in J.masks])
    n_eq = int(J.eq_mask.sum())
    m = len(J.residual(np.zeros(J.nvar)))
    eq_deg = J.deg[J.eq_mask] + 1
    if J.mode == "qs":
        row_top = np.concatenate([eq_deg] * 9 + [np.ones(6, dtype=int)])
    elif J.mode == "killing":
        row_top = np.concatenate([eq_deg] * 9 + [np.ones(6, dtype=int)] + [eq_deg] * 6)
    else:
        raise ValueError(J.mode)
    assert len(row_top) == m
    return var_level, row_top


def rank_and_left_null(M, gaps):
    """Rank of M and an orthonormal basis of its left null space (as columns)."""
    if M.shape[0] == 0:
        return 0, np.zeros((0, 0))
    if M.shape[1] == 0:
        return 0, np.eye(M.shape[0])
    U, s, _ = np.linalg.svd(M, full_matrices=True)
    if s[0] == 0:
        return 0, np.eye(M.shape[0])
    r, g = gap_rank(s)
    if r < len(s):
        gaps.append(g)
    return r, U[:, r:]


def row_basis(M, gaps):
    """Orthonormal rows spanning the row space of M."""
    if M.shape[0] == 0 or not np.any(M):
        return np.zeros((0, M.shape[1]))
    _, s, Vt = np.linalg.svd(M, full_matrices=False)
    r, g = gap_rank(s)
    if r < len(s):
        gaps.append(g)
    return Vt[:r]


def projected_counts(J, x, kmax):
    Jac = J.jacobian(x)
    var_level, row_top = levels(J)
    is_B = np.zeros(J.nvar, dtype=bool)
    is_B[: J.offsets[3]] = True
    gaps = []
    C = np.zeros((0, J.nvar))  # constraints carried down; zero on every eliminated unknown
    out = {}
    top = int(var_level.max())
    for L in range(top, 0, -1):
        high = var_level == L
        low = var_level < L
        M = np.vstack([Jac[row_top == L], C])
        _, N = rank_and_left_null(M[:, high], gaps)
        Cn = N.T @ M if N.size else np.zeros((0, J.nvar))
        Cn[:, ~low] = 0.0
        C = row_basis(Cn, gaps)
        k = L - 1
        if 1 <= k <= kmax:
            # all constraints on unknowns of degree <= k: original equations of lower degree, plus C
            Mk = np.vstack([Jac[row_top <= k], C])[:, low]
            Bk = is_B[low]
            _, N2 = rank_and_left_null(Mk[:, ~Bk], gaps)   # eliminate u, psi, p of degree <= k
            MB = N2.T @ Mk[:, Bk] if N2.size else np.zeros((0, int(Bk.sum())))
            rB = 0
            if MB.shape[0] and np.any(MB):
                s = np.linalg.svd(MB, compute_uv=False)
                rB, g = gap_rank(s)
                if rB < len(s):
                    gaps.append(g)
            out[k] = int(Bk.sum()) - rB
    return [out[k] for k in range(1, kmax + 1)], (min(gaps) if gaps else np.inf)


def ladder(mode, Kmax, kmax, seed=0, Kstart=3, out=print, ckpt_dir="runs"):
    """Carry one solution up order by order. Each order's solution is saved, so a later run
    (here or on another machine, after a git pull) resumes from the highest saved order."""
    import os
    os.makedirs(ckpt_dir, exist_ok=True)
    path = lambda K: os.path.join(ckpt_dir, f"local_{mode}_seed{seed}_K{K}.npy")
    rng = np.random.default_rng(seed)
    saved = [K for K in range(Kstart, Kmax + 1) if os.path.exists(path(K))]
    if saved:
        K = max(saved)
        J = Variant(K, mode)
        x = np.load(path(K))
        out(f"{mode}: resuming from saved order {K}")
    else:
        K = Kstart
        J = Variant(K, mode)
        for _ in range(20):
            x, nr, _ = J.solve(rng)
            if nr < 1e-10:
                break
        else:
            out(f"{mode}: no solution found at order {Kstart}")
            return
        np.save(path(K), x)
    while True:
        t0 = time.time()
        dims, gap = projected_counts(J, x, min(kmax, K - 1))
        out(f"{mode:8s} K={K:2d}  freedom in B, k=1.. : " + " ".join(f"{v:4d}" for v in dims)
            + f"   sym_B {J.sym_B(x, 2):8.2e}  kill_u {J.kill_u(x):8.2e}  residual {np.linalg.norm(J.residual(x)):7.1e}"
            + f"  min log10 gap {np.log10(gap):4.1f}   [{time.time() - t0:.0f}s]")
        if K == Kmax:
            return
        Jn = Variant(K + 1, mode)
        t0 = time.time()
        xn, nr, it = Jn.solve(rng, x0=Jn.extend_from(J, x, rng, noise=0.0))
        if nr > 1e-10:
            xn, nr, it = Jn.solve(rng, x0=Jn.extend_from(J, x, rng, noise=1e-2))
        if nr > 1e-10:
            out(f"{mode:8s} K={K + 1:2d}  could not be extended (residual {nr:.1e} after {it} iterations)")
            return
        J, x, K = Jn, xn, K + 1
        np.save(path(K), x)
        out(f"         extended to order {K} in {it} Newton steps, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    Kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    modes = sys.argv[2].split(",") if len(sys.argv) > 2 else ["qs"]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    for mode in modes:
        ladder(mode, Kmax, 9, seed=seed, out=lambda s: print(s, flush=True))
