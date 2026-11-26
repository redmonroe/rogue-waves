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

**Status, results so far, and next steps are in [`LOG.md`](LOG.md).** Read it first and add to it
before finishing a session.

### Running the tools

From `stellarator/`, with `PY=../.venv/bin/python`:

    $PY test_search.py                      # sanity check of the search objective (seconds)
    $PY test_leangen.py                     # the generator must reject a false identity
    $PY seeds.py 2                          # exact seed search, polynomial degree 2
    $PY local_jets_project.py 12 qs 0       # local Taylor-series test, stable count, orders 3 to 12
    $PY local_extend.py 13                  # carry the saved order-13 solution to order 14, with a log
    $PY local_planar_check.py               # run the local test on the exact non-symmetric solution
    $PY local_jets_project.py 11 strong 0   # local test with the strong form of quasisymmetry
    $PY planar_family.py                    # exact non-symmetric quasisymmetric families: check all equations
    $PY gen_twisted_pinch.py                # regenerate the Lean certificate for the twisted pinch
    $PY scan_iota2.py                       # scan the whole two-parameter iota=2 family (seconds)
    $PY qs_residuals.py 40                  # residual table for the explicit families (minutes)
    $PY qs_other_families.py                # same for the vacuum and N-field-period families
    $PY search.py --family sheared4 --hours 8   # overnight search, all cores, appends to runs/sheared4.jsonl
    $PY search.py --report runs/sheared4.jsonl  # summary: best residual at each asymmetry level
    $PY check_candidate.py sheared4 runs/sheared4.jsonl   # 40-digit check of the best result
    $PY gen_landreman_iota2.py              # regenerate the Lean certificate
    cd ../RogueWaves && lake build          # check every Lean proof

How to read a search report: every family member is an exact equilibrium, so only quasisymmetry
is being searched. The asymmetry parameter is held at or above a floor (0.2, 0.05, 0.01 by
default). If the best residual falls only as the floor falls, the family has no quasisymmetric
member away from the round case. A residual below 1e-6 at a fixed floor is a candidate and must
then pass `check_candidate.py`; the float search bottoms out near 1e-8.

### Known traps

- In current Mathlib `ContDiff ℝ ⊤` means real-analytic. Smooth is `∞` with `open scoped ContDiff`.
- Landreman's iota=2 family has closed field lines, so its flux label is not unique. Use the
  `frob` residual, which does not depend on the label, before concluding anything from `triple`.
- A failed search is not evidence of non-existence.
- A search family must have a fixed overall size. If the optimiser can shrink the sampled
  region relative to the torus, it will, and the residual falls without the field being symmetric.
- In the `deg2` family, parameters with `4*c02*c20 - c11**2 <= 0` give unbounded field lines and
  spurious small residuals (about 1e-6, not zero). The search excludes them.

## Layout

- `RogueWaves/` — Lake project (Mathlib, toolchain pinned in `lean-toolchain`)
  - `Stellarator/` — statement layer and certificates
  - `RogueWaves/Strassen.lean` — Strassen's 2×2 scheme, an early test of the setup
- `LOG.md` — what has been tried, results, next steps
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
