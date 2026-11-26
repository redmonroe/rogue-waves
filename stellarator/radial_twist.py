"""Strong quasisymmetry WITHOUT force balance, with a symmetry field that is not a rigid motion.

Burby, Kallinikos, MacKay, Perrella and Pfefferle (arXiv:2403.03352, J. Plasma Phys. 2025) define
B quasisymmetric with generator u when

    div u = 0,   div B = 0,   curl(u x B) = 0,   (curl B) x u + grad(u . B) = 0,

show that a non-rigid ("non-Killing") u must satisfy  tr S = 0, det S = 0, [e3, u] = 0  with
S the symmetric part of grad u and e3 its null direction, and that then B is parallel to e3. They
leave open whether any non-Killing u satisfying these exists.

Candidate (same symmetry field as the twisted pinch in planar_family.py):

    u = d/dy + theta'(y) * (rotation about the y axis),      any function theta(y),
    B = G(phi - theta(y)) * (x, 0, z) / (x^2 + z^2),         any function G,  phi = angle about the y axis.

This script checks every equation symbolically for arbitrary theta and G, checks the three
conditions on u, and computes the strain of u (non-zero wherever theta'' is non-zero).
It also checks what fails (force balance), and for one explicit rational case (tan(theta/2) = y,
G = cos^2) measures at 40 digits whether any rigid motion leaves B unchanged.

Why B has no rigid symmetry when theta'' != 0 and G' != 0 (by hand): the direction of B is the
radial direction about the y axis, so a rigid motion preserving B preserves that axis, so it is
a d/dy + b d/dphi; that changes B by G' (b - a theta'(y)) / r, zero only if theta' is constant.

Usage:  python radial_twist.py
"""
import mpmath as mp
import sympy as sp

from planar_family import rigid_symmetry

x, y, z = X = sp.symbols("x y z", real=True)
theta = sp.Function("theta")(y)
G = sp.Function("G")

grad = lambda f: sp.Matrix([sp.diff(f, v) for v in X])
div = lambda V: sum(sp.diff(V[i], X[i]) for i in range(3))
curl = lambda V: sp.Matrix([sp.diff(V[(i + 2) % 3], X[(i + 1) % 3]) - sp.diff(V[(i + 1) % 3], X[(i + 2) % 3])
                            for i in range(3)])
simp = lambda M: sp.simplify(M.doit()) if hasattr(M, "doit") else sp.simplify(M)

r2 = x**2 + z**2
phi = sp.atan2(x, z)                       # angle about the y axis, increasing from z towards x
w = sp.Matrix([z, 0, -x])                  # rotation about the y axis that DEcreases phi ... sign fixed below
s = sp.simplify(w.dot(grad(phi)))          # = +1 or -1
u = sp.Matrix([0, 1, 0]) + sp.diff(theta, y) * s * w   # u . grad(phi - theta) = 0
arg = phi - theta
B = G(arg) * sp.Matrix([x, 0, z]) / r2

print("u . grad(phi - theta) :", simp(u.dot(grad(arg))))
print("div u                 :", simp(div(u)))
print("div B                 :", simp(div(B)))
print("curl(u x B)           :", list(simp(curl(u.cross(B)))))
strong = curl(B).cross(u) + grad(u.dot(B))
print("(curl B) x u + grad(u.B):", list(simp(strong)))
print("u . grad |B|^2        :", simp(u.dot(grad(B.dot(B)))))

# conditions on u alone
Ju = u.jacobian(X)
S = (Ju + Ju.T) / 2
print("strain S              :", S.applyfunc(sp.simplify).tolist())
print("tr S, det S           :", simp(S.trace()), simp(S.det()))
e3 = sp.Matrix([x, 0, z])                  # claimed null direction of S (radial)
print("S . e3                :", list(simp(S * e3)))
lie = e3.jacobian(X) * u - Ju * e3          # [u, e3]
print("[u, e3]               :", list(simp(lie)))
print("S : S                 :", sp.factor(simp(sum(S[i, j]**2 for i in range(3) for j in range(3)))))

# what fails: force balance. J x B must be a gradient, so its curl must vanish.
JxB = curl(B).cross(B)
print("curl(J x B)           :", [sp.factor(c) for c in simp(curl(JxB))])
print("J = lambda u, lambda  :", sp.factor(simp(curl(B)[1] / u[1])))

# explicit rational case and controls: rigid-symmetry measure (1e-40 = symmetric, 0.01+ = not)
mp.mp.dps = 40
W = 1 + y**2
c, sn = (1 - y**2) / W, 2 * y / W


def field(cos_t, sin_t, power):
    q = cos_t * z + sin_t * x                # r cos(phi - theta)
    return (q**2 / r2) ** power * sp.Matrix([x, 0, z]) / r2


off = {x: x + sp.Rational(3, 2), z: z + sp.Rational(1, 5)}   # sample box away from the axis
k = sp.Rational(7, 10)
print()
for name, Bc in [("tan(theta/2) = y, G = cos^2          ", field(c, sn, 1)),
                 ("theta = 0.7 y (control: screw)       ", field(sp.cos(k * y), sp.sin(k * y), 1)),
                 ("tan(theta/2) = y, G = 1 (control)    ", field(c, sn, 0))]:
    print(f"rigid-symmetry measure, {name}: {mp.nstr(rigid_symmetry(Bc.subs(off, simultaneous=True)), 3)}")


# ---------------------------------------------------------------------------------------------
# Second family, same recipe. Recipe: u = K0 + f K1 with K0, K1 commuting rigid motions and f a
# function with K1 . grad f = 0 and K0 . grad f a function of f. Then S is the symmetrised
# product of grad f and K1, its null direction is grad f x K1, and B must point along it.
# Here K0 = d/dx, K1 = d/dz, f = a(y) exp(kappa x):
#     u = (1, 0, a(y) exp(kappa x)),   B = Phi(a(y) exp(kappa x) - kappa z) * (-a'(y)/a(y), kappa, 0)
# for any functions a and Phi. Field lines lie in the planes z = const.
print("\nsecond family: u = d/dx + a(y) exp(kappa x) d/dz")
kap = sp.symbols("kappa", positive=True)
a = sp.Function("a")(y)
Phi = sp.Function("Phi")
f = a * sp.exp(kap * x)
u2 = sp.Matrix([1, 0, f])
B2 = Phi(f - kap * z) * sp.Matrix([-sp.diff(a, y) / a, kap, 0])
print("div u, div B          :", simp(div(u2)), simp(div(B2)))
print("curl(u x B)           :", list(simp(curl(u2.cross(B2)))))
print("(curl B) x u + grad(u.B):", list(simp(curl(B2).cross(u2) + grad(u2.dot(B2)))))
J2 = u2.jacobian(X)
S2 = (J2 + J2.T) / 2
print("S : S                 :", sp.factor(simp(sum(S2[i, j]**2 for i in range(3) for j in range(3)))))
conc = lambda av, ph: sp.Matrix([v.subs(kap, 1).replace(Phi, ph).subs(a, av).doit() for v in B2])
shift = {x: x + sp.Rational(1, 5), y: y + sp.Rational(3, 2)}
for name, Bc in [("a = y, Phi = square                 ", conc(y, lambda t: t**2 + 1)),
                 ("a = exp(y) (control: one direction) ", conc(sp.exp(y), lambda t: t**2 + 1))]:
    print(f"rigid-symmetry measure, {name}: {mp.nstr(rigid_symmetry(Bc.subs(shift, simultaneous=True)), 3)}")
