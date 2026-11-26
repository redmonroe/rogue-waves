import Mathlib

/-!
# Landreman's closed-form 3D equilibrium with rotational transform 2

Field from Landreman, arXiv:2609.26742, eqs. (2.1)-(2.2), with the two stretch factors `a`, `b`
kept as free nonzero parameters (the paper takes `a = √(1+ε)`, `b = √(1-ε)`).

  s = x²/a² + y²/b²,   G = 1 - (1 - s)² - 4z²,   F = √G,
  B = ( (2zx - (a/b) F y)/s , (2zy + (b/a) F x)/s , 1 - s ).

Proved here: `∇·B = 0` wherever `s > 0` and `G > 0`, with the divergence written as the sum of
the three one-variable partial derivatives.
Not yet proved here: force balance, and the link to `Stellarator.div` from the statement layer.
-/

noncomputable section

namespace Stellarator.Landreman

/-- `s = x²/a² + y²/b²`. -/
def s (a b x y : ℝ) : ℝ := x ^ 2 / a ^ 2 + y ^ 2 / b ^ 2

/-- `G = 1 - (1 - s)² - 4z²`; the field is defined where `G > 0`. -/
def G (a b x y z : ℝ) : ℝ := 1 - (1 - s a b x y) ^ 2 - 4 * z ^ 2

def Bx (a b x y z : ℝ) : ℝ := (2 * z * x - a / b * Real.sqrt (G a b x y z) * y) / s a b x y
def By (a b x y z : ℝ) : ℝ := (2 * z * y + b / a * Real.sqrt (G a b x y z) * x) / s a b x y
def Bz (a b x y : ℝ) : ℝ := 1 - s a b x y

/-- `∇·B = 0` on the domain `s > 0`, `G > 0`. -/
theorem div_free (a b x y z : ℝ) (ha : a ≠ 0) (hb : b ≠ 0)
    (hs : 0 < s a b x y) (hG : 0 < G a b x y z) :
    deriv (fun t => Bx a b t y z) x + deriv (fun t => By a b x t z) y
      + deriv (fun _ : ℝ => Bz a b x y) z = 0 := by
  have hI : Real.sqrt (G a b x y z) ≠ 0 := (Real.sqrt_pos.mpr hG).ne'
  have hGne : G a b x y z ≠ 0 := hG.ne'
  have hsne : s a b x y ≠ 0 := hs.ne'
  have ha2 : a ^ 2 ≠ 0 := pow_ne_zero 2 ha
  have hb2 : b ^ 2 ≠ 0 := pow_ne_zero 2 hb
  -- x-derivatives
  have hsx : HasDerivAt (fun t : ℝ => s a b t y) (2 * x / a ^ 2) x := by
    have h := ((hasDerivAt_pow 2 x).div_const (a ^ 2)).add_const (y ^ 2 / b ^ 2)
    simpa [s] using h
  have hGx : HasDerivAt (fun t : ℝ => G a b t y z) (2 * (1 - s a b x y) * (2 * x / a ^ 2)) x := by
    have h := (((hsx.const_sub 1).pow 2).const_sub 1).sub_const (4 * z ^ 2)
    have h' : HasDerivAt (fun t : ℝ => G a b t y z) _ x := h
    refine h'.congr_deriv ?_
    simp only [Nat.cast_ofNat, Nat.add_one_sub_one, pow_one]
    ring
  have hBx : HasDerivAt (fun t : ℝ => Bx a b t y z)
      (((2 * z - a / b * ((2 * (1 - s a b x y) * (2 * x / a ^ 2)) / (2 * Real.sqrt (G a b x y z))) * y)
          * s a b x y
        - (2 * z * x - a / b * Real.sqrt (G a b x y z) * y) * (2 * x / a ^ 2)) / s a b x y ^ 2) x := by
    have h1 : HasDerivAt (fun t : ℝ => 2 * z * t) (2 * z) x := by
      simpa using (hasDerivAt_id x).const_mul (2 * z)
    have h2 := ((hGx.sqrt hGne).const_mul (a / b)).mul_const y
    exact (h1.sub h2).div hsx hsne
  -- y-derivatives
  have hsy : HasDerivAt (fun t : ℝ => s a b x t) (2 * y / b ^ 2) y := by
    have h := ((hasDerivAt_pow 2 y).div_const (b ^ 2)).const_add (x ^ 2 / a ^ 2)
    simpa [s] using h
  have hGy : HasDerivAt (fun t : ℝ => G a b x t z) (2 * (1 - s a b x y) * (2 * y / b ^ 2)) y := by
    have h := (((hsy.const_sub 1).pow 2).const_sub 1).sub_const (4 * z ^ 2)
    have h' : HasDerivAt (fun t : ℝ => G a b x t z) _ y := h
    refine h'.congr_deriv ?_
    simp only [Nat.cast_ofNat, Nat.add_one_sub_one, pow_one]
    ring
  have hBy : HasDerivAt (fun t : ℝ => By a b x t z)
      (((2 * z + b / a * ((2 * (1 - s a b x y) * (2 * y / b ^ 2)) / (2 * Real.sqrt (G a b x y z))) * x)
          * s a b x y
        - (2 * z * y + b / a * Real.sqrt (G a b x y z) * x) * (2 * y / b ^ 2)) / s a b x y ^ 2) y := by
    have h1 : HasDerivAt (fun t : ℝ => 2 * z * t) (2 * z) y := by
      simpa using (hasDerivAt_id y).const_mul (2 * z)
    have h2 := ((hGy.sqrt hGne).const_mul (b / a)).mul_const x
    exact (h1.add h2).div hsy hsne
  have hBz : deriv (fun _ : ℝ => Bz a b x y) z = 0 := by simp
  rw [hBx.deriv, hBy.deriv, hBz]
  have hsdef : s a b x y = x ^ 2 / a ^ 2 + y ^ 2 / b ^ 2 := rfl
  field_simp
  rw [hsdef]
  field_simp
  ring

end Stellarator.Landreman

#print axioms Stellarator.Landreman.div_free
