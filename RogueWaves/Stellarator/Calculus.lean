import Stellarator.Defs

/-! # Basic identities for the hand-rolled operators: `∇·(∇×B) = 0` and `∇×(∇f) = 0`. -/

noncomputable section

namespace Stellarator

lemma pd_comm {f : V3 → ℝ} (hf : ContDiff ℝ 2 f) (i j : Fin 3) (x : V3) :
    pd i (pd j f) x = pd j (pd i f) x := by
  have hd : DifferentiableAt ℝ (fderiv ℝ f) x :=
    ((hf.fderiv_right (m := 1) (by norm_num)).differentiable (by norm_num)).differentiableAt
  have hsymm : IsSymmSndFDerivAt ℝ f x :=
    hf.contDiffAt.isSymmSndFDerivAt (by simp)
  unfold pd
  have h1 : ∀ k : Fin 3, fderiv ℝ (fun y => fderiv ℝ f y (Pi.single k 1)) x
      = (fderiv ℝ (fderiv ℝ f) x).flip (Pi.single k 1) := by
    intro k
    ext v
    rw [fderiv_clm_apply hd (differentiableAt_const _)]
    simp
  rw [h1 j, h1 i]
  simp only [ContinuousLinearMap.flip_apply]
  exact hsymm _ _

theorem div_curl {B : V3 → V3} (hB : ContDiff ℝ 2 B) (x : V3) : div (curl B) x = 0 := by
  have hc : ∀ k : Fin 3, ContDiff ℝ 2 (fun y => B y k) := fun k => (contDiff_pi.mp hB) k
  have e01 := pd_comm (hc 2) 0 1 x
  have e02 := pd_comm (hc 1) 0 2 x
  have e12 := pd_comm (hc 0) 1 2 x
  have hd : ∀ (k i : Fin 3), Differentiable ℝ (pd i (fun y => B y k)) := by
    intro k i
    have : ContDiff ℝ 1 (fun y => fderiv ℝ (fun y => B y k) y (Pi.single i 1)) :=
      ((hc k).fderiv_right (m := 1) (by norm_num)).clm_apply contDiff_const
    exact this.differentiable (by norm_num)
  simp only [div, curl, Fin.sum_univ_three, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons]
  have sub : ∀ (i : Fin 3) (f g : V3 → ℝ), Differentiable ℝ f → Differentiable ℝ g →
      pd i (fun y => f y - g y) x = pd i f x - pd i g x := by
    intro i f g hf hg
    unfold pd
    rw [fderiv_fun_sub hf.differentiableAt hg.differentiableAt]
    rfl
  rw [sub 0 _ _ (hd 2 1) (hd 1 2), sub 1 _ _ (hd 0 2) (hd 2 0), sub 2 _ _ (hd 1 0) (hd 0 1)]
  linarith

theorem curl_grad {f : V3 → ℝ} (hf : ContDiff ℝ 2 f) (x : V3) : curl (grad f) x = 0 := by
  ext i
  fin_cases i <;> simp [curl, grad] <;>
    first
      | (rw [sub_eq_zero]; exact pd_comm hf _ _ x)
      | (rw [sub_eq_zero]; exact (pd_comm hf _ _ x).symm)

end Stellarator
