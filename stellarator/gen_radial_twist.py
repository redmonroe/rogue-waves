"""Generate the Lean certificate for the radial twisted field (see radial_twist.py).

With W = 1 + y^2, c = (1 - y^2)/W, s = 2y/W (cosine and sine of theta(y) = 2 arctan y),
q = c z + s x and R = x^2 + z^2:

  B = q^2 (x, 0, z) / R^2          (radial about the y axis, strength cos^2(phi - theta(y)) / r)
  u = ( 2z/W , 1 , -2x/W )         (move up in y while turning about the y axis at rate theta')

Writes RogueWaves/Stellarator/Generated/RadialTwist.lean with:
  div_free                       div B = 0
  div_u                          div u = 0
  frozen_x, frozen_y, frozen_z   curl(u x B) = 0
  perp                           u . B = 0
  strong_x, strong_y, strong_z   (curl B) x u + grad(u . B) = 0   (gradient written out by the product rule)
  qs                             u . grad |B|^2 = 0
  strain_xy                      d(u_x)/dy = -4 y z / W^2   (so u is not a rigid motion)
B has no y-component and u_y = 1; the statements are written with those substituted.
"""
import os

import sympy as sp

from leangen import LeanGen

x, y, z = sp.symbols("x y z")

g = LeanGen(
    "Stellarator.Generated.RadialTwist", [], [x, y, z],
    header="# A quasisymmetric field (strong sense, no force balance) whose symmetry is not a rigid motion\n\n"
           "With `W = 1 + y²`, `c = (1 - y²)/W`, `s = 2y/W`, `q = c z + s x`, `R = x² + z²`:\n\n"
           "  B = q² (x, 0, z) / R²,\n"
           "  u = ( 2z/W , 1 , -2x/W ).\n\n"
           "B points straight out from the y axis; its strength is `cos²(φ - θ(y)) / r` with\n"
           "`θ(y) = 2 arctan y`. The flow of u moves up in y while turning about the y axis by `θ`.\n"
           "Because `θ` is not linear in `y`, u is not a rigid motion (`strain_xy`), and B has no\n"
           "translation, rotation or screw symmetry (checked numerically in\n"
           "`stellarator/radial_twist.py`, argued by hand there, not proved here).\n\n"
           "Theorems: `∇·B = 0`; `∇·u = 0`; `∇×(u×B) = 0`; `(∇×B)×u + ∇(u·B) = 0`; `u·∇|B|² = 0`.\n"
           "These are the quasisymmetry conditions of Burby, Kallinikos and MacKay (2020), without\n"
           "force balance. Force balance does not hold for this field. Stated with one-variable\n"
           "partial derivatives and with `B_y = 0`, `u_y = 1` substituted. Not yet linked to the\n"
           "operators in `Stellarator/Defs.lean`.",
)
W = g.atom("W", 1 + y**2)
R = g.atom("R", x**2 + z**2)
c, s = (1 - y**2) / W, 2 * y / W
q = c * z + s * x
Bx, Bz = q**2 * x / R**2, q**2 * z / R**2
ux, uz = 2 * z / W, -2 * x / W
g.define("Bx", Bx, "x-component of the field")
g.define("Bz", Bz, "z-component of the field (the y-component is zero)")
g.define("ux", ux, "x-component of the symmetry vector field")
g.define("uz", uz, "z-component of the symmetry vector field (the y-component is one)")
g.define("Ey", uz * Bx - ux * Bz, "y-component of u x B (the others are B_z and -B_x)")
g.define("N", Bx**2 + Bz**2, "|B|^2")

g.identity("div_free", lambda v, D: D("Bx", x) + D("Bz", z), "The field is divergence-free.")
g.identity("div_u", lambda v, D: D("ux", x) + D("uz", z), "The symmetry vector field is divergence-free.")
# u x B = (Bz, Ey, -Bx)
g.identity("frozen_x", lambda v, D: -D("Bx", y) - D("Ey", z), "x-component of curl(u x B) = 0.")
g.identity("frozen_y", lambda v, D: D("Bz", z) + D("Bx", x), "y-component of curl(u x B) = 0.")
g.identity("frozen_z", lambda v, D: D("Ey", x) - D("Bz", y), "z-component of curl(u x B) = 0.")
# curl B with B_y = 0
J = lambda D: (D("Bz", y), D("Bx", z) - D("Bz", x), -D("Bx", y))
# grad(u . B) written out with the product rule (u_y = 1, B_y = 0)
dC = lambda v, D, w: D("ux", w) * v("Bx") + v("ux") * D("Bx", w) + D("uz", w) * v("Bz") + v("uz") * D("Bz", w)
g.identity("perp", lambda v, D: v("ux") * v("Bx") + v("uz") * v("Bz"), "u . B = 0.")
g.identity("strong_x", lambda v, D: J(D)[1] * v("uz") - J(D)[2] + dC(v, D, x),
           "x-component of (curl B) x u + grad(u . B) = 0.")
g.identity("strong_y", lambda v, D: J(D)[2] * v("ux") - J(D)[0] * v("uz") + dC(v, D, y),
           "y-component of (curl B) x u + grad(u . B) = 0.")
g.identity("strong_z", lambda v, D: J(D)[0] - J(D)[1] * v("ux") + dC(v, D, z),
           "z-component of (curl B) x u + grad(u . B) = 0.")
g.identity("qs", lambda v, D: v("ux") * D("N", x) + D("N", y) + v("uz") * D("N", z),
           "The field strength is constant along u.")
g.identity("strain_xy", lambda v, D: D("ux", y) + 4 * y * z / W**2,
           "d(u_x)/dy = -4yz/W^2, not identically zero: u is not a rigid motion (u_y = 1 is constant).")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "RogueWaves", "Stellarator", "Generated", "RadialTwist.lean")
open(out, "w").write(g.render())
print("wrote", os.path.normpath(out))
