"""Generate the Lean certificate for Landreman's iota=2 equilibrium (arXiv:2609.26742, eqs. 2.1-2.4).

Writes RogueWaves/Stellarator/Generated/LandremanIota2.lean with:
  div_free                 d(Bx)/dx + d(By)/dy + d(Bz)/dz = 0
  force_x, force_y, force_z   ((curl B) x B)_i - d(P)/dx_i = 0
where P = -(x^2 + y^2 + 4 z^2 + |B|^2)/2 is the paper's pressure with its constant dropped.
"""
import os
import sys

import sympy as sp

from leangen import LeanGen

a, b, x, y, z = sp.symbols("a b x y z")

g = LeanGen(
    "Stellarator.Generated.LandremanIota2", [a, b], [x, y, z],
    header="# Landreman's closed-form 3D equilibrium with rotational transform 2\n\n"
           "Field and pressure from Landreman, arXiv:2609.26742, eqs. (2.1)-(2.4), with stretch factors\n"
           "`a`, `b` free (the paper takes `a = √(1+ε)`, `b = √(1-ε)`). With\n"
           "`s = x²/a² + y²/b²` and `G = 1 - (1 - s)² - 4z²`:\n\n"
           "  B = ( (2zx - (a/b) √G y)/s , (2zy + (b/a) √G x)/s , 1 - s ),\n"
           "  P = -(x² + y² + 4z² + |B|²)/2.\n\n"
           "Theorems: `∇·B = 0` and the three components of `(∇×B)×B = ∇P`, each stated with\n"
           "one-variable partial derivatives, on the set where the listed hypotheses hold.",
)
S = g.atom("S", x**2 / a**2 + y**2 / b**2)
F = g.sqrt("F", 1 - (1 - S) ** 2 - 4 * z**2)
Bx = (2 * z * x - (a / b) * F * y) / S
By = (2 * z * y + (b / a) * F * x) / S
Bz = 1 - S
P = -(x**2 + y**2 + 4 * z**2 + Bx**2 + By**2 + Bz**2) / 2
g.define("Bx", Bx, "x-component of the field")
g.define("By", By, "y-component of the field")
g.define("Bz", Bz, "z-component of the field")
g.define("P", P, "pressure, up to an additive constant")

g.identity("div_free", lambda v, D: D("Bx", x) + D("By", y) + D("Bz", z), "The field is divergence-free.")

J = lambda D: (D("Bz", y) - D("By", z), D("Bx", z) - D("Bz", x), D("By", x) - D("Bx", y))
only = sys.argv[1:] or ["x", "y", "z"]
if "x" in only:
    g.identity("force_x", lambda v, D: J(D)[1] * v("Bz") - J(D)[2] * v("By") - D("P", x), "x-component of force balance.")
if "y" in only:
    g.identity("force_y", lambda v, D: J(D)[2] * v("Bx") - J(D)[0] * v("Bz") - D("P", y), "y-component of force balance.")
if "z" in only:
    g.identity("force_z", lambda v, D: J(D)[0] * v("By") - J(D)[1] * v("Bx") - D("P", z), "z-component of force balance.")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "RogueWaves", "Stellarator", "Generated", "LandremanIota2.lean")
os.makedirs(os.path.dirname(out), exist_ok=True)
open(out, "w").write(g.render())
print("wrote", os.path.normpath(out))
