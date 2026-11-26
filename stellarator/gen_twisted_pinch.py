"""Generate the Lean certificate for the twisted elliptical pinch (see planar_family.py).

The field, with W = 1 + y^2, c = (1 - y^2)/W, s = 2y/W (cosine and sine of the twist angle
theta(y) = 2 arctan y), x' = c x + s z, z' = -s x + c z:

  B   = ( -2 beta c z' - 2 alpha s x' , 0 , -2 beta s z' + 2 alpha c x' )
  u   = ( -theta' z , 1 , theta' x ),   theta' = 2 / W
  psi = -(alpha x'^2 + beta z'^2),      P = 2 (alpha + beta) psi

Writes RogueWaves/Stellarator/Generated/TwistedPinch.lean with:
  div_free                    div B = 0
  force_x, force_y, force_z   (curl B) x B = grad P
  div_u                       div u = 0
  flux_x, flux_y, flux_z      B x u = grad psi
  qs                          u . grad |B|^2 = 0
B has no y-component, so the statements are written with that component already set to zero.
"""
import os

import sympy as sp

from leangen import LeanGen

al, be, x, y, z = sp.symbols("alpha beta x y z")

g = LeanGen(
    "Stellarator.Generated.TwistedPinch", [al, be], [x, y, z],
    header="# A quasisymmetric equilibrium with no rigid-motion symmetry: the twisted elliptical pinch\n\n"
           "With `W = 1 + y²`, `c = (1 - y²)/W`, `s = 2y/W`, `x' = c x + s z`, `z' = -s x + c z`:\n\n"
           "  B   = ( -2β c z' - 2α s x' , 0 , -2β s z' + 2α c x' ),\n"
           "  u   = ( -2z/W , 1 , 2x/W ),\n"
           "  ψ   = -(α x'² + β z'²),   P = 2(α + β) ψ.\n\n"
           "In each plane `y = const` this is the linear field with elliptical field lines\n"
           "`α x² + β z² = const`, turned by the angle `2 arctan y`. Because that angle is not linear in\n"
           "`y`, the field has no translation, rotation or screw symmetry (checked numerically in\n"
           "`stellarator/planar_family.py`, not proved here).\n\n"
           "Theorems: `∇·B = 0`; `(∇×B)×B = ∇P`; `∇·u = 0`; `B×u = ∇ψ`; `u·∇|B|² = 0`. The last three are\n"
           "weak quasisymmetry. Stated with one-variable partial derivatives and with `B_y = 0`\n"
           "substituted. Not yet linked to the operators in `Stellarator/Defs.lean`.",
)
W = g.atom("W", 1 + y**2)
c, s = (1 - y**2) / W, 2 * y / W
xp, zp = c * x + s * z, -s * x + c * z
g.define("Bx", -2 * be * c * zp - 2 * al * s * xp, "x-component of the field")
g.define("Bz", -2 * be * s * zp + 2 * al * c * xp, "z-component of the field (the y-component is zero)")
g.define("ux", -2 * z / W, "x-component of the symmetry vector field")
g.define("uz", 2 * x / W, "z-component of the symmetry vector field (the y-component is one)")
g.define("psi", -(al * xp**2 + be * zp**2), "flux label")
g.define("P", -2 * (al + be) * (al * xp**2 + be * zp**2), "pressure, up to an additive constant")
g.define("N", (-2 * be * c * zp - 2 * al * s * xp) ** 2 + (-2 * be * s * zp + 2 * al * c * xp) ** 2, "|B|^2")

g.identity("div_free", lambda v, D: D("Bx", x) + D("Bz", z), "The field is divergence-free.")
# curl B with B_y = 0
J = lambda D: (D("Bz", y), D("Bx", z) - D("Bz", x), -D("Bx", y))
g.identity("force_x", lambda v, D: J(D)[1] * v("Bz") - D("P", x), "x-component of force balance.")
g.identity("force_y", lambda v, D: J(D)[2] * v("Bx") - J(D)[0] * v("Bz") - D("P", y), "y-component of force balance.")
g.identity("force_z", lambda v, D: -J(D)[1] * v("Bx") - D("P", z), "z-component of force balance.")
g.identity("div_u", lambda v, D: D("ux", x) + D("uz", z), "The symmetry vector field is divergence-free.")
g.identity("flux_x", lambda v, D: -v("Bz") - D("psi", x), "x-component of B x u = grad psi.")
g.identity("flux_y", lambda v, D: v("Bz") * v("ux") - v("Bx") * v("uz") - D("psi", y), "y-component of B x u = grad psi.")
g.identity("flux_z", lambda v, D: v("Bx") - D("psi", z), "z-component of B x u = grad psi.")
g.identity("qs", lambda v, D: v("ux") * D("N", x) + D("N", y) + v("uz") * D("N", z),
           "The field strength is constant along u: quasisymmetry.")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "RogueWaves", "Stellarator", "Generated", "TwistedPinch.lean")
open(out, "w").write(g.render())
print("wrote", os.path.normpath(out))
