"""The generator must refuse a false identity. Run: python test_leangen.py"""
import sympy as sp

from leangen import LeanGen

a, b, x, y, z = sp.symbols("a b x y z")


def make(pressure_coeff):
    g = LeanGen("Test", [a, b], [x, y, z])
    S = g.atom("S", x**2 / a**2 + y**2 / b**2)
    F = g.sqrt("F", 1 - (1 - S) ** 2 - 4 * z**2)
    Bx, By, Bz = (2 * z * x - (a / b) * F * y) / S, (2 * z * y + (b / a) * F * x) / S, 1 - S
    g.define("Bx", Bx); g.define("By", By); g.define("Bz", Bz)
    g.define("P", -(x**2 + y**2 + pressure_coeff * z**2 + Bx**2 + By**2 + Bz**2) / 2)
    J = lambda D: (D("Bz", y) - D("By", z), D("Bx", z) - D("Bz", x), D("By", x) - D("Bx", y))
    g.identity("force_z", lambda v, D: J(D)[0] * v("By") - J(D)[1] * v("Bx") - D("P", z))
    return g


make(4)  # the true pressure: must succeed
try:
    make(3)  # wrong pressure: must be rejected
except ValueError as e:
    print("rejected false identity:", str(e)[:70])
else:
    raise SystemExit("FAILED: a false identity was accepted")
print("ok")
