import Mathlib

/-!
# Statement layer: vector calculus on ℝ³, MHD equilibrium, quasisymmetry

This file contains definitions only. It is the file a reader has to trust.
Everything is stated on an open set `U ⊆ ℝ³`, so these are local (pointwise) notions:
nothing here says that the pressure surfaces are nested tori.

Sources:
* MHD equilibrium: Landreman, arXiv:2609.26742, eq. (1.1).
* Weak quasisymmetry: Constantin, Drivas, Ginsberg, arXiv:2009.08860, Definition 1.1 (weak form);
  Rodríguez, Helander, Bhattacharjee, arXiv:2004.11431.
* Triple-product form: Rodríguez, Helander, Bhattacharjee, arXiv:2004.11431.

Regularity is `∞` (smooth), written with `open scoped ContDiff`. In current Mathlib `⊤` would
mean real-analytic, which is a different and stronger condition.
-/

open scoped ContDiff

noncomputable section

namespace Stellarator

/-- Points and vectors of ℝ³. -/
abbrev V3 := Fin 3 → ℝ

/-- Partial derivative of `f` in coordinate direction `i`. -/
def pd (i : Fin 3) (f : V3 → ℝ) (x : V3) : ℝ := fderiv ℝ f x (Pi.single i 1)

/-- Gradient of a scalar field. -/
def grad (f : V3 → ℝ) (x : V3) : V3 := fun i => pd i f x

/-- Divergence of a vector field. -/
def div (B : V3 → V3) (x : V3) : ℝ := ∑ i, pd i (fun y => B y i) x

/-- Curl of a vector field. -/
def curl (B : V3 → V3) (x : V3) : V3 :=
  ![pd 1 (fun y => B y 2) x - pd 2 (fun y => B y 1) x,
    pd 2 (fun y => B y 0) x - pd 0 (fun y => B y 2) x,
    pd 0 (fun y => B y 1) x - pd 1 (fun y => B y 0) x]

/-- Dot product. -/
def dot (u v : V3) : ℝ := ∑ i, u i * v i

/-- Cross product. -/
def cross (u v : V3) : V3 :=
  ![u 1 * v 2 - u 2 * v 1, u 2 * v 0 - u 0 * v 2, u 0 * v 1 - u 1 * v 0]

/-- Squared field strength `|B|²`. -/
def normSq (B : V3 → V3) (x : V3) : ℝ := dot (B x) (B x)

/-- Magnetohydrostatic equilibrium on an open set `U`: `B` and `p` are smooth on `U`,
`∇·B = 0`, and `(∇×B)×B = ∇p`. -/
structure IsMHSOn (B : V3 → V3) (p : V3 → ℝ) (U : Set V3) : Prop where
  isOpen : IsOpen U
  smooth_B : ContDiffOn ℝ ∞ B U
  smooth_p : ContDiffOn ℝ ∞ p U
  div_free : ∀ x ∈ U, div B x = 0
  force : ∀ x ∈ U, cross (curl B x) (B x) = grad p x

/-- `ψ` is a flux label for `B` on `U`: smooth, constant along `B`, with nonvanishing gradient. -/
structure IsFluxLabelOn (B : V3 → V3) (ψ : V3 → ℝ) (U : Set V3) : Prop where
  smooth : ContDiffOn ℝ ∞ ψ U
  tangent : ∀ x ∈ U, dot (B x) (grad ψ x) = 0
  grad_ne_zero : ∀ x ∈ U, grad ψ x ≠ 0

/-- Weak quasisymmetry on `U` with symmetry field `u` and flux label `ψ`:
`∇·u = 0`, `B × u = ∇ψ`, and `u·∇|B|² = 0`.
(`u·∇|B|² = 0` is equivalent to `u·∇|B| = 0` wherever `B ≠ 0`, and avoids a square root.) -/
structure IsWeakQSOn (B u : V3 → V3) (ψ : V3 → ℝ) (U : Set V3) : Prop where
  smooth_u : ContDiffOn ℝ ∞ u U
  div_free : ∀ x ∈ U, div u x = 0
  cross_eq : ∀ x ∈ U, cross (B x) (u x) = grad ψ x
  strength : ∀ x ∈ U, dot (u x) (grad (normSq B) x) = 0

/-- Triple-product form of quasisymmetry for a flux label `ψ`:
`∇ψ × ∇|B|² · ∇(B·∇|B|²) = 0` on `U`. -/
def IsQSTripleOn (B : V3 → V3) (ψ : V3 → ℝ) (U : Set V3) : Prop :=
  ∀ x ∈ U, dot (cross (grad ψ x) (grad (normSq B) x))
    (grad (fun y => dot (B y) (grad (normSq B) y)) x) = 0

/-- `B` is invariant under rotation about the z-axis, in the sense used for axisymmetric
equilibria: its cylindrical components do not depend on the toroidal angle. Stated through the
rotation generator `e_z × x = (-y, x, 0)`: the Lie derivative of `B` along it vanishes. -/
def IsAxisymmetricOn (B : V3 → V3) (U : Set V3) : Prop :=
  ∀ x ∈ U, ∀ i : Fin 3,
    (-(x 1)) * pd 0 (fun y => B y i) x + (x 0) * pd 1 (fun y => B y i) x
      = (![-(B x 1), B x 0, 0] : V3) i

end Stellarator
