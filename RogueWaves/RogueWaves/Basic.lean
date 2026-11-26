import Mathlib

/-- Smoke test: Mathlib loads and tactics work. -/
example (a b : ℝ) (h : 0 < a) (hb : 0 < b) : 0 < a * b := by positivity

/-- A rogue wave: height exceeds twice the significant wave height. -/
def IsRogue (H Hs : ℝ) : Prop := H > 2 * Hs

example : IsRogue 5 2 := by unfold IsRogue; norm_num
