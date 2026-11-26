"""Overnight search for a quasisymmetric member of a family of exact 3D equilibria.

Every family here is an exact MHD equilibrium for all parameter values, so the only thing being
searched is quasisymmetry. The objective is the label-independent residual

    frob = | N . curl N | / (|N| |curl N|),   N = B x (grad|B|^2 x grad(B . grad|B|^2)),

root-mean-squared over sample points. It must vanish for a weakly quasisymmetric field. Each
family has an asymmetry parameter that is held at or above --asym-min, so the optimiser cannot
escape to the axisymmetric members, where the residual is trivially zero.

Runs are deterministic (seeded), parallel, resumable, and append one JSON line per restart.
A restart whose best residual is below --flag is marked as a candidate; candidates should then
be checked at high precision with qs_residuals.py and by exact algebra. A run with no candidates
is NOT evidence that none exists.

Examples:
    python search.py --family deg2 --restarts 64 --workers 8
    python search.py --family deg2 --hours 8 --workers 8 --asym-min 0.02
    python search.py --report runs/deg2.jsonl
"""
import argparse
import json
import math
import multiprocessing as mp
import os
import time

import numpy as np
import sympy as sp
from scipy.optimize import differential_evolution, minimize

x, y, z = sp.symbols("x y z", real=True)
X = (x, y, z)


def _grad(f):
    return sp.Matrix([sp.diff(f, v) for v in X])


def _re(e):
    return (e + e.subs(sp.I, -sp.I)) / 2


def _im(e):
    return (e - e.subs(sp.I, -sp.I)) / (2 * sp.I)


# ---------------------------------------------------------------- families

def family_deg2():
    """Stretched degree-2 polynomial seeds (see seeds.py), A = diag(a, b, 1).
    Landreman's iota=2 family is c02=-1, c10=1/2, c11=0, c20=-1/4."""
    c02, c10, c11, c20, m, al = sp.symbols("c02 c10 c11 c20 m alpha", real=True)
    a, b = m * (1 + al), m * (1 - al)
    u, Z = sp.symbols("u Z")
    Psi = Z**2 * c02 + Z * c10 * c11 / (2 * c20) + Z * c11 * u + c10 * u + c20 * u**2
    d1 = 4 * c02 * c20 - c11**2
    I2 = -sp.diff(Psi, Z) ** 2 + 2 * u * (sp.diff(Psi, Z) * sp.diff(Psi, u, Z) - sp.diff(Psi, u) * sp.diff(Psi, Z, 2)) + d1 * u**2
    xi, eta = x / a, y / b
    sub = {u: xi**2 + eta**2, Z: z}
    PZ, Pu, I2x = sp.diff(Psi, Z).subs(sub), sp.diff(Psi, u).subs(sub), I2.subs(sub)
    I = sp.sqrt(I2x)
    uu = xi**2 + eta**2
    B = sp.Matrix([a * (-PZ * xi - I * eta) / uu, b * (-PZ * eta + I * xi) / uu, 2 * Pu])
    return dict(
        params=[c02, c10, c11, c20, m, al],
        bounds=[(-2, 2), (-2, 2), (-2, 2), (-2, 2), (0.5, 2.0), None],  # None = asymmetry slot
        B=B, valid=[I2x, uu], box=[(-2.0, 2.0), (-2.0, 2.0), (-1.0, 1.0)],
        # d1 > 0: field lines are harmonic oscillators (bounded). d1 <= 0 gives unbounded field
        # lines, where the optimiser otherwise finds meaningless small residuals.
        param_valid=[d1],
        reference=dict(c02=-1, c10=0.5, c11=0, c20=-0.25, m=1.0),
    )


def family_sheared():
    """Landreman's sheared-iota family, eqs. (3.1)-(3.3). Asymmetry parameter is eps."""
    lam, S, eps = sp.symbols("lambda S epsilon", positive=True)
    w, wb = x + sp.I * y, x - sp.I * y
    rad = 1 + eps / wb**2
    K = wb * sp.sqrt(rad)
    Xi = w * K + sp.pi / 2 - S
    h = sp.I * sp.exp(-sp.I * lam * z) * sp.sin(Xi) / (2 * K)
    c = _re(sp.exp(-sp.I * lam * z) * sp.cos(Xi))
    B = sp.Matrix([_re(h), _im(h), c / lam])
    return dict(
        params=[lam, S, eps], bounds=[(0.3, 3.0), (0.0, math.pi), None],
        B=B, valid=[_re(rad), x**2 + y**2], box=[(-1.5, 1.5), (-1.5, 1.5), (-1.0, 1.0)],
        reference=dict(lam=1.0, S=1.0),
    )


FAMILIES = dict(deg2=family_deg2, sheared=family_sheared)

# ---------------------------------------------------------------- objective

_STATE = {}


def _init(family_name, n_points, seed):
    fam = FAMILIES[family_name]()
    B = fam["B"]
    B2 = (B.T * B)[0]
    BdB2 = (B.T * _grad(B2))[0]
    N = B.cross(_grad(B2).cross(_grad(BdB2)))
    args = (*X, *fam["params"])
    rng = np.random.default_rng(seed)
    pts = np.column_stack([rng.uniform(lo, hi, n_points) for lo, hi in fam["box"]])
    _STATE.update(
        fam=fam, pts=pts,
        N=sp.lambdify(args, list(N), "numpy", cse=True),
        valid=sp.lambdify(args, fam["valid"], "numpy", cse=True),
        param_valid=sp.lambdify(fam["params"], fam.get("param_valid", []), "math"),
    )


def _N(P, params):
    out = _STATE["N"](P[:, 0], P[:, 1], P[:, 2], *params)
    return np.real(np.stack([np.broadcast_to(np.asarray(c, dtype=complex), P.shape[0]) for c in out], axis=1))


def objective(params, h=1e-5, margin=0.05, min_valid=40, detail=False):
    P = _STATE["pts"]
    if any(not (c > margin) for c in _STATE["param_valid"](*params)):
        return (1.0, 0) if detail else 1.0
    with np.errstate(all="ignore"):
        v = _STATE["valid"](P[:, 0], P[:, 1], P[:, 2], *params)
        ok = np.ones(P.shape[0], dtype=bool)
        for c in v:
            ok &= np.real(np.broadcast_to(np.asarray(c, dtype=complex), P.shape[0])) > margin
        if ok.sum() < min_valid:
            return (1.0, 0) if detail else 1.0
        Q = P[ok]
        N0 = _N(Q, params)
        J = np.empty((Q.shape[0], 3, 3))  # J[:, i, j] = d N_i / d x_j
        for j in range(3):
            e = np.zeros(3)
            e[j] = h
            J[:, :, j] = (_N(Q + e, params) - _N(Q - e, params)) / (2 * h)
        curlN = np.stack([J[:, 2, 1] - J[:, 1, 2], J[:, 0, 2] - J[:, 2, 0], J[:, 1, 0] - J[:, 0, 1]], axis=1)
        nN, nC = np.linalg.norm(N0, axis=1), np.linalg.norm(curlN, axis=1)
        good = np.isfinite(nN) & np.isfinite(nC) & (nN > 1e-9) & (nC > 1e-9)
        if good.sum() < min_valid:
            return (1.0, 0) if detail else 1.0  # field strength structure is degenerate here
        r = np.abs(np.sum(N0[good] * curlN[good], axis=1)) / (nN[good] * nC[good])
        val = float(np.sqrt(np.mean(r**2)))
    if not np.isfinite(val):
        val = 1.0
    return (val, int(good.sum())) if detail else val


def _restart(job):
    index, asym_min, maxiter, seed = job
    fam = _STATE["fam"]
    bounds = [(asym_min, 0.8) if b is None else b for b in fam["bounds"]]
    t0 = time.time()
    de = differential_evolution(objective, bounds, seed=seed + index, maxiter=maxiter, popsize=12,
                                tol=1e-10, polish=False, init="sobol")
    nm = minimize(objective, de.x, method="Nelder-Mead",
                  options=dict(xatol=1e-10, fatol=1e-14, maxiter=4000))
    xbest = np.clip(nm.x, [b[0] for b in bounds], [b[1] for b in bounds])
    val, n_used = objective(xbest, detail=True)
    return dict(index=index, asym_min=asym_min, residual=val, n_points=n_used,
                params={str(s): float(v) for s, v in zip(fam["params"], xbest)},
                seconds=round(time.time() - t0, 1))


# ---------------------------------------------------------------- driver

def report(path, flag):
    rows = [json.loads(line) for line in open(path)]
    rows.sort(key=lambda r: r["residual"])
    print(f"{len(rows)} restarts in {path}")
    by = {}
    for r in rows:
        by.setdefault(r["asym_min"], []).append(r["residual"])
    print("best residual by minimum asymmetry (a residual that shrinks only as asymmetry shrinks is a dead end):")
    for k in sorted(by):
        print(f"  asym >= {k:<6g}  best {min(by[k]):.3e}   median {np.median(by[k]):.3e}   ({len(by[k])} restarts)")
    cands = [r for r in rows if r["residual"] < flag]
    print(f"candidates below {flag:g}: {len(cands)}")
    for r in rows[:5]:
        print(f"  {r['residual']:.3e}  asym>={r['asym_min']:g}  {json.dumps(r['params'])}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family", choices=sorted(FAMILIES), default="deg2")
    ap.add_argument("--restarts", type=int, default=16, help="number of restarts (ignored with --hours)")
    ap.add_argument("--hours", type=float, default=None, help="keep launching restarts for this long")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--asym-min", type=float, nargs="+", default=[0.2, 0.05, 0.01])
    ap.add_argument("--maxiter", type=int, default=60)
    ap.add_argument("--points", type=int, default=3000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--flag", type=float, default=1e-6)
    ap.add_argument("--out", default=None)
    ap.add_argument("--report", default=None, metavar="JSONL")
    a = ap.parse_args()
    if a.report:
        return report(a.report, a.flag)
    out = a.out or f"runs/{a.family}.jsonl"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    done = sum(1 for _ in open(out)) if os.path.exists(out) else 0
    deadline = time.time() + a.hours * 3600 if a.hours else None

    def jobs():
        i = done
        while (deadline and time.time() < deadline) or (not deadline and i < done + a.restarts):
            yield (i, a.asym_min[i % len(a.asym_min)], a.maxiter, a.seed)
            i += 1

    print(f"family {a.family}: resuming after {done} restarts, writing {out}", flush=True)
    with mp.Pool(a.workers, initializer=_init, initargs=(a.family, a.points, a.seed)) as pool, open(out, "a") as f:
        for r in pool.imap_unordered(_restart, jobs()):
            f.write(json.dumps(r) + "\n")
            f.flush()
            tag = "  <-- CANDIDATE" if r["residual"] < a.flag else ""
            print(f"restart {r['index']:4d}  asym>={r['asym_min']:<5g} residual {r['residual']:.3e}  [{r['seconds']}s]{tag}", flush=True)
    report(out, a.flag)


if __name__ == "__main__":
    main()
