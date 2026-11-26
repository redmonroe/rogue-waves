"""Grid scan of the whole normalised iota=2 family (two parameters: eps, tau0).

For each asymmetry eps, prints the smallest quasisymmetry residual over tau0 in [-0.9, 0.9].
If min residual / eps stays roughly constant as eps -> 0, the residual vanishes only in the
axisymmetric limit: the family has no quasisymmetric member away from it (numerically).
"""
import numpy as np

import search

search._init("iota2", 6000, 0)
taus = np.linspace(-0.9, 0.9, 37)
print(f"{'eps':>8s} {'min residual':>14s} {'at tau0':>8s} {'max residual':>14s} {'min / eps':>10s} {'points':>7s}")
for eps in [0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
    vals = [search.objective([t, eps], detail=True) for t in taus]
    res = np.array([v[0] for v in vals])
    i = int(res.argmin())
    print(f"{eps:8.3f} {res[i]:14.3e} {taus[i]:8.2f} {res.max():14.3e} {res[i] / eps:10.3f} {vals[i][1]:7d}")
