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
| Residual tool (`stellarator/qs_residuals.py`) | Working. Checks equilibrium and two quasisymmetry residuals at 40-digit precision |
| Landreman's iota=2 family (arXiv:2609.26742, sec. 2) | Exact equilibrium; **not quasisymmetric** at eps = 1/4, 1/2 (numerical, not a proof) |
| Landreman's sheared family (sec. 3) | Exact equilibrium; **not quasisymmetric** at the two parameter sets tried (numerical) |
| Other published closed-form families (arXiv:2610.07304, 2610.29921, 2611.07672, 2611.31603) | Not yet transcribed or tested |
| Seed search (`stellarator/seeds.py`) | Working. Exact algebra: finds every polynomial axisymmetric seed of a given degree whose stretched versions are 3D equilibria. Degree 2 gives a 4-parameter family containing Landreman's seed; degrees 2 and 3 only give frequency ratio 2 |
| Overnight search (`stellarator/search.py`) | Working: parallel, seeded, resumable. Only short test runs so far (3 restarts, best residual 3e-3, no candidates). **No real search has been run yet** |
| Candidate checker (`stellarator/check_candidate.py`) | Working. Re-checks a search result at 40 digits |
| Lean proof generator (`stellarator/leangen.py`) | Working for fields built from + - * /, powers, named rational sub-expressions and one square root. No sin/cos/exp yet, so it cannot do the sheared family |
| Lean: Landreman iota=2 field is divergence-free and in force balance (`Stellarator/Generated/LandremanIota2.lean`) | **Proved**, generated, standard axioms only. Stated with one-variable partial derivatives; not yet linked to the operators in `Defs.lean` |
| Lean statement layer (`Stellarator/Defs.lean`) | First draft, builds. **Not reviewed, and no definition has been tested against an example yet** |
| Lean: `div curl = 0`, `curl grad = 0` (`Calculus.lean`) | Proved |
| Lean: axisymmetric equilibrium satisfies the quasisymmetry definition (the "should pass" test) | Not started |
| Calibration run (rediscover a known family from Landreman's prompt, closed book) | Started once, lost when the cloud machine restarted; not yet rerun |

Output of the residual table is in `stellarator/qs_residuals_output.txt`. Columns: `mhs` should be
about zero (the formulas are transcribed correctly); `triple` is zero only if the field is weakly
quasisymmetric with its flux label; `frob` must be zero if the field is weakly quasisymmetric
with any flux label. The axisymmetric member of each family is a control and gives zeros.

### Running the tools

From `stellarator/`, with `PY=../.venv/bin/python`:

    $PY test_search.py                      # sanity check of the search objective (seconds)
    $PY test_leangen.py                     # the generator must reject a false identity
    $PY seeds.py 2                          # exact seed search, polynomial degree 2
    $PY search.py --family deg2 --hours 8   # overnight search, all cores, appends to runs/deg2.jsonl
    $PY search.py --report runs/deg2.jsonl  # summary: best residual at each asymmetry level
    $PY check_candidate.py deg2 runs/deg2.jsonl   # 40-digit check of the best result
    $PY gen_landreman_iota2.py              # regenerate the Lean certificate
    cd ../RogueWaves && lake build          # check every Lean proof

How to read a search report: every family member is an exact equilibrium, so only quasisymmetry
is being searched. The asymmetry parameter is held at or above a floor (0.2, 0.05, 0.01 by
default). If the best residual falls only as the floor falls, the family has no quasisymmetric
member away from the round case. A residual below 1e-6 at a fixed floor is a candidate and must
then pass `check_candidate.py`; the float search bottoms out near 1e-8.

### Next steps, in order

1. Run the overnight search on `deg2` and read the report.
2. Transcribe and test the remaining published closed-form families.
3. Lean: link the generated certificate to the operators in `Defs.lean`, giving a first
   `IsMHSOn` instance (this also shows the definition is not empty).
4. Lean: prove the axisymmetric member satisfies `IsQSTripleOn` and `IsWeakQSOn`; confirm by
   computer algebra that the stretched member fails. This is the two-sided test of the definitions.
5. Extend the proof generator to sin/cos/exp so it can certify the sheared family.
6. New recipes for the search (the current families have at most six parameters).
7. Have the statement file reviewed adversarially, then by an outside reader.
8. Hard stop at the end of the second usage reset unless a candidate's residual vanishes identically.

### Known traps

- In current Mathlib `ContDiff ℝ ⊤` means real-analytic. Smooth is `∞` with `open scoped ContDiff`.
- Landreman's iota=2 family has closed field lines, so its flux label is not unique. Use the
  `frob` residual, which does not depend on the label, before concluding anything from `triple`.
- A failed search is not evidence of non-existence.
- In the `deg2` family, parameters with `4*c02*c20 - c11**2 <= 0` give unbounded field lines and
  spurious small residuals (about 1e-6, not zero). The search excludes them.

## Layout

- `RogueWaves/` — Lake project (Mathlib, toolchain pinned in `lean-toolchain`)
  - `Stellarator/` — statement layer and certificates
  - `RogueWaves/Strassen.lean` — Strassen's 2×2 scheme, an early test of the setup
- `stellarator/` — Python tools (SymPy, mpmath, NumPy, SciPy); see "Running the tools"
- `research_notes/`, `reports/` — background research for both candidate projects
- `fetch_data.sh` — rogue wave dataset download (6.2 GB into `data/`, not tracked)
- `setup.sh` — one-time setup on a new machine
- `code/`, `data/`, `.venv/` — not tracked

## Setup on a new machine

    bash setup.sh

Installs the Lean toolchain and prebuilt Mathlib and creates `.venv` with SymPy, mpmath, NumPy and SciPy.
Add `--data` to also download the rogue wave dataset (6.2 GB).

Then:

    cd RogueWaves && lake build
    cd ../stellarator && ../.venv/bin/python qs_residuals.py 40

## Notes

- `lakefile.toml` pulls Mathlib by git URL instead of through the Lean package index.
- The setup script was written for macOS/Linux and has not been run end to end.
