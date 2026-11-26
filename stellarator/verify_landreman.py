"""Check Landreman's closed-form 3D equilibrium (arXiv:2609.26742, eqs. 2.1-2.4).

1. Symbolic: div B = 0 and (curl B) x B = grad p.
2. Quasisymmetry test: T = grad(psi) x grad(|B|^2) . grad(B . grad |B|^2).
   T vanishes identically iff the field is (weakly) quasisymmetric.
"""
import random
import sympy as sp

x, y, z = sp.symbols("x y z", real=True)
# a, b are kept as free positive symbols: the construction B = A B0(A^-1 x) works for any
# A = diag(a, b, 1). The paper's family is a = sqrt(1+eps), b = sqrt(1-eps).
a, b = sp.symbols("a b", positive=True)
X = (x, y, z)

s = x**2 / a**2 + y**2 / b**2
F = sp.sqrt(1 - (1 - s) ** 2 - 4 * z**2)
B = sp.Matrix([(2 * z * x - (a / b) * F * y) / s, (2 * z * y + (b / a) * F * x) / s, 1 - s])

grad = lambda f: sp.Matrix([sp.diff(f, v) for v in X])
div = lambda V: sum(sp.diff(V[i], X[i]) for i in range(3))
curl = lambda V: sp.Matrix([
    sp.diff(V[2], y) - sp.diff(V[1], z),
    sp.diff(V[0], z) - sp.diff(V[2], x),
    sp.diff(V[1], x) - sp.diff(V[0], y),
])

B2 = (B.T * B)[0]
psi = (x**2 + y**2 + 4 * z**2 + B2 - 2) / 4  # constant offset dropped
p = -2 * psi  # p_a = 0

print("div B            =", sp.simplify(div(B)))
force = curl(B).cross(B) - grad(p)
print("(curl B)xB - grad p =", [sp.simplify(c) for c in force])
print("B . grad psi     =", sp.simplify((B.T * grad(psi))[0]))

# Quasisymmetry triple product, evaluated numerically at random interior points.
BdB2 = (B.T * grad(B2))[0]
T = grad(psi).cross(grad(B2)).dot(grad(BdB2))
T_fn = sp.lambdify((x, y, z, a, b), T, "mpmath")
F2_fn = sp.lambdify((x, y, z, a, b), 1 - (1 - s) ** 2 - 4 * z**2, "mpmath")
scale_fn = sp.lambdify((x, y, z, a, b),
                       sp.sqrt(grad(psi).dot(grad(psi)) * grad(B2).dot(grad(B2)) * grad(BdB2).dot(grad(BdB2))),
                       "mpmath")

random.seed(0)
for e in (0, sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1, 2)):
    worst = 0
    n = 0
    while n < 200:
        pt = [random.uniform(-1.3, 1.3), random.uniform(-1.3, 1.3), random.uniform(-0.45, 0.45)]
        av, bv = sp.sqrt(1 + e).evalf(30), sp.sqrt(1 - e).evalf(30)
        if F2_fn(*pt, av, bv) <= 0.05:
            continue
        n += 1
        worst = max(worst, abs(T_fn(*pt, av, bv)) / scale_fn(*pt, av, bv))
    print(f"eps = {str(e):>5}: max |T| / (|grad psi||grad B^2||grad(B.grad B^2)|) over 200 points = {float(worst):.3e}")
