# rogue-waves

Working repository for a Lean 4 + AI attempt at an open problem. The current target is
**exact quasisymmetry**: is there a non-axisymmetric magnetic field that is in MHD force balance
and exactly quasisymmetric? (The repository name is from an earlier candidate project on rogue
wave prediction, which was set aside; its data script and notes are still here.)

## Current target: exact quasisymmetry

Background and the case for and against are in
`reports/Stellarator quasisymmetry Lean viability.md`. Short version: Grad's conjecture was
refuted in September 2026 with closed-form 3D equilibria; whether any such equilibrium can be
exactly quasisymmetric is still open, and experts expect not. The plan is a time-boxed search
for a quasisymmetric closed-form equilibrium, with Lean certifying anything found.

### Status

| Item | State |
|---|---|
| Residual tool for closed-form fields (`stellarator/qs_residuals.py`) | Working. Checks MHD equilibrium and two quasisymmetry residuals at 40-digit precision |
| Landreman's iota=2 family (arXiv:2609.26742, sec. 2) | Exact equilibrium confirmed; **not quasisymmetric** at eps = 1/4, 1/2 (numerical, not a proof) |
| Landreman's sheared family (sec. 3) | Exact equilibrium confirmed; **not quasisymmetric** at the two parameter sets tried (numerical) |
| Other published closed-form families (vacuum, force-free, broken stellarator symmetry, any field periods) | Not yet transcribed or tested |
| Lean statement layer (`RogueWaves/Stellarator/Defs.lean`) | First draft, builds. **Not reviewed by anyone, and no definition has been tested against an example yet** |
| Lean: `div curl = 0`, `curl grad = 0` (`Calculus.lean`) | Proved |
| Lean: Landreman iota=2 field is divergence-free (`Landreman.lean`) | Proved, in one-variable partial-derivative form; not yet linked to `Stellarator.div` |
| Lean: force balance for that field | Not started |
| Lean: axisymmetric equilibrium satisfies the quasisymmetry definition (the "should pass" test) | Not started |
| Calibration run (can the model rediscover a known family from Landreman's prompt, closed book) | See `stellarator/calibration/` if present; otherwise not yet recorded |
| Search for a quasisymmetric closed-form equilibrium | Not started |

Output of the residual table is in `stellarator/qs_residuals_output.txt`. Columns: `mhs` should be
about zero (the formulas are transcribed correctly); `triple` is zero only if the field is weakly
quasisymmetric with its flux label; `frob` must be zero if the field is weakly quasisymmetric
with any flux label. The axisymmetric member of each family is a control and gives zeros.

### Next steps, in order

1. Transcribe and test the remaining published closed-form families (arXiv:2610.07304,
   2610.29921, 2611.07672, 2611.31603).
2. Lean: force balance for Landreman's iota=2 field; link both results to the operators in
   `Defs.lean`, giving a first `IsMHSOn` instance (this also shows the definition is not empty).
3. Lean: prove the axisymmetric member satisfies `IsQSTripleOn` and `IsWeakQSOn`; confirm by
   computer algebra that the stretched member fails. This is the two-sided test of the definitions.
4. Have the statement file reviewed adversarially, then by an outside reader.
5. Search: deformation families of axisymmetric seeds, with the residuals as the test.
6. Hard stop at the end of the second usage reset unless a candidate's residual vanishes identically.

### Known traps

- In current Mathlib `ContDiff ℝ ⊤` means real-analytic. Smooth is `∞` with `open scoped ContDiff`.
- Landreman's iota=2 family has closed field lines, so its flux label is not unique. Use the
  `frob` residual, which does not depend on the label, before concluding anything from `triple`.
- A failed search is not evidence of non-existence.

## Layout

- `RogueWaves/` — Lake project (Mathlib, toolchain pinned in `lean-toolchain`)
  - `Stellarator/` — statement layer and certificates
  - `RogueWaves/Strassen.lean` — Strassen's 2×2 scheme, an early test of the setup
- `stellarator/` — SymPy/mpmath scripts: `fields.py`, `qs_residuals.py`, `verify_family1_symbolic.py`
- `research_notes/`, `reports/` — background research for both candidate projects
- `fetch_data.sh` — rogue wave dataset download (6.2 GB into `data/`, not tracked)
- `setup.sh` — one-time setup on a new machine
- `code/`, `data/`, `.venv/` — not tracked

## Setup on a new machine

    bash setup.sh

Installs the Lean toolchain and prebuilt Mathlib and creates `.venv` with SymPy and mpmath.
Add `--data` to also download the rogue wave dataset (6.2 GB).

Then:

    cd RogueWaves && lake build
    cd ../stellarator && ../.venv/bin/python qs_residuals.py 40

## Notes

- `lakefile.toml` pulls Mathlib by git URL instead of through the Lean package index.
- The setup script was written for macOS/Linux and has not been run end to end.
