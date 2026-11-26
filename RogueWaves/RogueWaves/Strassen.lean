import Mathlib

/-!
Strassen's scheme: a 2×2 matrix product using 7 multiplications instead of 8.
Proved over any ring, commutative or not, so it applies recursively to block matrices.
-/

theorem strassen {R : Type*} [Ring R] (a11 a12 a21 a22 b11 b12 b21 b22 : R) :
    let m1 := (a11 + a22) * (b11 + b22)
    let m2 := (a21 + a22) * b11
    let m3 := a11 * (b12 - b22)
    let m4 := a22 * (b21 - b11)
    let m5 := (a11 + a12) * b22
    let m6 := (a21 - a11) * (b11 + b12)
    let m7 := (a12 - a22) * (b21 + b22)
    (m1 + m4 - m5 + m7 = a11 * b11 + a12 * b21) ∧
    (m3 + m5 = a11 * b12 + a12 * b22) ∧
    (m2 + m4 = a21 * b11 + a22 * b21) ∧
    (m1 - m2 + m3 + m6 = a21 * b12 + a22 * b22) := by
  dsimp only
  refine ⟨?_, ?_, ?_, ?_⟩ <;> noncomm_ring

#print axioms strassen
