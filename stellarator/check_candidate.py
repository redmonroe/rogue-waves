"""High-precision check of a search candidate.

Rebuilds the family at the candidate's parameters (floats converted to exact rationals) and
evaluates, with 40-digit arithmetic at random valid points:
  mhs  : |div B| and the residual of (curl B) x B = -grad(Q + |B|^2/2) is NOT assumed; instead
         curl((curl B) x B) is checked to vanish, which holds iff (curl B) x B is a gradient
  frob : the label-independent quasisymmetry residual (see qs_residuals.py)
Usage:  python check_candidate.py deg2 '{"c02": -1, "c10": 0.5, "c11": 0, "c20": -0.25, "m": 1, "alpha": 0.25}'
        python check_candidate.py deg2 runs/deg2.jsonl        (checks the best line in the file)
"""
import json
import os
import random
import sys

import mpmath as mp
import sympy as sp

import search

mp.mp.dps = 40
x, y, z = search.X


def check(family, params, n_points=15, seed=1):
    fam = search.FAMILIES[family]()
    sub = {s: sp.nsimplify(params[str(s)], rational=True) for s in fam["params"]}
    B = fam["B"].subs(sub)
    g = lambda f: sp.Matrix([sp.diff(f, v) for v in search.X])
    curl = lambda V: sp.Matrix([sp.diff(V[2], y) - sp.diff(V[1], z), sp.diff(V[0], z) - sp.diff(V[2], x),
                                sp.diff(V[1], x) - sp.diff(V[0], y)])
    divB = sum(sp.diff(B[i], search.X[i]) for i in range(3))
    force = curl(B).cross(B)
    B2 = (B.T * B)[0]
    N = B.cross(g(B2).cross(g((B.T * g(B2))[0])))
    f_div = sp.lambdify(search.X, divB, "mpmath", cse=True)
    f_force = sp.lambdify(search.X, list(force), "mpmath", cse=True)
    f_N = sp.lambdify(search.X, list(N), "mpmath", cse=True)
    f_valid = sp.lambdify(search.X, [v.subs(sub) for v in fam["valid"]], "mpmath")

    def num_curl(fn, pt):
        d = lambda i, j: mp.diff(lambda t: fn(*[t if k == j else pt[k] for k in range(3)])[i], pt[j])
        return [d(2, 1) - d(1, 2), d(0, 2) - d(2, 0), d(1, 0) - d(0, 1)]

    nrm = lambda v: mp.sqrt(sum(abs(c) ** 2 for c in v))
    rng = random.Random(seed)
    worst = dict(div=mp.mpf(0), force_not_gradient=mp.mpf(0), frob=mp.mpf(0))
    n = tries = 0
    while n < n_points and tries < 200000:
        tries += 1
        pt = [mp.mpf(rng.uniform(lo, hi)) for lo, hi in fam["box"]]
        try:
            if any(mp.re(v) <= 0.05 for v in f_valid(*pt)):
                continue
            Nv = f_N(*pt)
        except (ZeroDivisionError, ValueError):
            continue
        n += 1
        worst["div"] = max(worst["div"], abs(f_div(*pt)))
        cf = num_curl(f_force, pt)
        worst["force_not_gradient"] = max(worst["force_not_gradient"], nrm(cf))
        cN = num_curl(f_N, pt)
        if nrm(Nv) > 0 and nrm(cN) > 0:
            worst["frob"] = max(worst["frob"], abs(mp.re(sum(a * b for a, b in zip(Nv, cN)))) / (nrm(Nv) * nrm(cN)))
    return {k: float(v) for k, v in worst.items()}, n


if __name__ == "__main__":
    family, arg = sys.argv[1], sys.argv[2]
    if os.path.exists(arg):
        rows = sorted((json.loads(l) for l in open(arg)), key=lambda r: r["residual"])
        params = rows[0]["params"]
        print("best line: residual", rows[0]["residual"])
    else:
        params = json.loads(arg)
    res, n = check(family, params)
    print(json.dumps(params))
    print(f"{n} points:", "  ".join(f"{k} = {v:.2e}" for k, v in res.items()))
