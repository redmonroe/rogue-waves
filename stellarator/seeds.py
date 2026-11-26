"""Exact search for axisymmetric seed fields whose stretched versions are 3D equilibria.

Recipe (generalising Landreman, arXiv:2609.26742, sec. 2). Write an axisymmetric field as
    B0 = grad(Psi) x grad(phi) + I grad(phi),   Psi = Psi(u, Z),  u = R^2.
If its tension is linear,  (B0 . grad) B0 = -D x  with  D = diag(d1, d1, d3),  then for any
A = diag(a, b, 1) the stretched field  B(x) = A B0(A^-1 x)  is divergence-free and satisfies
(curl B) x B = grad p  with  p = -(x.D.x + |B|^2)/2.  So every such seed gives a family of
exact 3D equilibria, non-axisymmetric when a != b.

For polynomial Psi the tension condition is a polynomial system in the coefficients:
    E1:  -4 Psi_Z Psi_uu + 4 Psi_u Psi_uZ + d3 Z = 0          (vertical tension)
    I^2 = -Psi_Z^2 + 2u (Psi_Z Psi_uZ - Psi_u Psi_ZZ) + d1 u^2  (radial tension defines I^2)
    E2:  Psi_u d(I^2)/dZ - Psi_Z d(I^2)/du = 0                (I^2 is constant on flux surfaces)
This script solves that system exactly with SymPy for Psi of a given degree.

Usage:  python seeds.py [degree]      (default 2; degree 3 can take a long time)
"""
import json
import sys
import time

import sympy as sp

u, Z = sp.symbols("u Z")
d1, d3 = sp.symbols("d1 d3")


def I2_of(Psi, d1v):
    Pu, PZ = sp.diff(Psi, u), sp.diff(Psi, Z)
    return sp.expand(-PZ**2 + 2 * u * (PZ * sp.diff(Psi, u, Z) - Pu * sp.diff(Psi, Z, 2)) + d1v * u**2)


def seed_equations(Psi):
    Pu, PZ = sp.diff(Psi, u), sp.diff(Psi, Z)
    E1 = sp.expand(-4 * PZ * sp.diff(Psi, u, 2) + 4 * Pu * sp.diff(Psi, u, Z) + d3 * Z)
    I2 = I2_of(Psi, d1)
    E2 = sp.expand(Pu * sp.diff(I2, Z) - PZ * sp.diff(I2, u))
    eqs = []
    for E in (E1, E2):
        eqs += [c for c in sp.Poly(E, u, Z).coeffs() if c != 0]
    return eqs


def search(degree):
    coeffs = {}
    Psi = 0
    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            if i == 0 and j == 0:
                continue
            c = sp.Symbol(f"c{i}{j}")
            coeffs[(i, j)] = c
            Psi += c * u**i * Z**j
    eqs = seed_equations(Psi)
    unknowns = list(coeffs.values()) + [d1, d3]
    sols = sp.solve(eqs, unknowns, dict=True)
    out = []
    for s in sols:
        P = sp.expand(Psi.subs(s))
        I2 = sp.factor(I2_of(P, d1.subs(s) if d1 in s else d1))
        if P == 0 or sp.diff(P, u) == 0 or sp.diff(P, Z) == 0:
            continue  # no poloidal field, or not a genuine function of both u and Z
        out.append(dict(Psi=str(P), I2=str(I2), d1=str(s.get(d1, d1)), d3=str(s.get(d3, d3)),
                        free=[str(v) for v in sorted(P.free_symbols | I2.free_symbols - {u, Z}, key=str) if v not in (u, Z)]))
    return out


# Landreman's seed (eq. 2.6), in this notation.
LANDREMAN = dict(Psi="u/2 - u**2/4 - Z**2", I2="2*u - u**2 - 4*Z**2", d1="1", d3="4")

if __name__ == "__main__":
    degree = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    t0 = time.time()
    sols = search(degree)
    print(f"degree {degree}: {len(sols)} seed families with both u- and Z-dependence ({time.time() - t0:.0f}s)")
    for s in sols:
        print(json.dumps(s))
    with open(f"seeds_degree{degree}.json", "w") as f:
        json.dump(sols, f, indent=1)
