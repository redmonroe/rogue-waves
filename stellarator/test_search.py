"""Sanity checks for search.py: the objective must be ~0 on axisymmetric members and clearly
non-zero on Landreman's non-axisymmetric members. Run: python test_search.py"""
import search

for fam, zero, nonzero in [
    ("deg2", [-1, 0.5, 0, -0.25, 1.0, 0.0], [-1, 0.5, 0, -0.25, 1.0, 0.25]),
    ("sheared", [1.0, 1.0, 1e-12], [1.0, 1.0, 0.5]),
]:
    search._init(fam, 3000, 0)
    r0, n0 = search.objective(zero, detail=True)
    r1, n1 = search.objective(nonzero, detail=True)
    print(f"{fam:8s} axisymmetric: {r0:.2e} ({n0} pts)   non-axisymmetric: {r1:.2e} ({n1} pts)")
    assert n0 >= 40 and n1 >= 40, "too few valid sample points"
    assert r0 < 1e-5, "control should be ~0"
    assert r1 > 1e-2, "known non-quasisymmetric member should be clearly non-zero"
print("ok")
