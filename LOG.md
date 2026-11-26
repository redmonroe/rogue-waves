# Project log

One place for what has been tried, what it showed, and what comes next. Newest entries at the
bottom of the log. Any session (cloud or Mac) should read this first and add to it before
finishing.

How to add an entry: date, what was tried, the result in one or two lines, where the evidence is
(file or commit), and how far the result can be trusted. If an entry is later found to be wrong,
do not delete it: add a correction line under it.

## Goal

Find a magnetic field that is in force balance and exactly quasisymmetric without being
axisymmetric (or helically symmetric), or find out why none exists. Anything found gets a Lean
certificate. Agreed limit: stop at the end of the second usage reset unless a candidate's
symmetry error vanishes identically.

## Next steps

In order. Move an item to the log when it is done.

1. **Local test to higher order.** Finish order 14 in the cloud session and order 16 on the Mac
   (`local_jets_project.py`). Read two things: whether the counts stop falling, and whether
   `min log10 gap` stays above about 3.
2. **If the gap closes before the counts settle:** redo the elimination in extended precision.
   Not built yet.
3. **Work out which hidden conditions are appearing.** The count at degree 5 loses 2 or 3 per
   order. Identifying the condition by hand would say whether the fall must continue.
4. **Recipe type 2: symmetry by construction.** Build fields in which every field line on a
   surface is the same curve moved rigidly, so the field-strength pattern is shared by
   construction, then search for force balance. Not started.
5. **Remaining published families.** Sengupta et al. vacuum family (arXiv:2610.28466) not yet
   read. Their force-free family (arXiv:2611.31603) cannot be tested until the companion paper
   with its formulas appears.
6. **Lean.** Link the generated certificate to the operators in `Stellarator/Defs.lean`; prove
   the axisymmetric member satisfies the quasisymmetry definitions and check the stretched member
   fails. This is the two-sided test of the definitions, which are otherwise unreviewed.
7. **Calibration run.** Give Landreman's original prompt to a fresh session, closed book, and see
   whether it rediscovers his field. Started once, lost in a machine restart.
8. **Outside reader** for the Lean statement file.

## Open questions

- Is the local question (force balance + quasisymmetry near a point, no torus) already settled
  in the literature? One arXiv search returned nothing usable. Klotz, Duignan, Burby and Meiss
  (arXiv:2609.09517) treat the case without force balance.
- Which family was the Mac run labelled "sengupta"? Its output was pasted but not pushed.

## Log

### 2026-11-26

**Rogue wave prediction: set aside.**
First candidate project. Research verdict: a real open problem with no identifiable user; worth
at most a short audit. Data script and notes remain in the repo.
Evidence: `reports/Rogue wave prediction viability.md`.

**Stellarator quasisymmetry: chosen.**
Research verdict: Grad's conjecture was refuted in September 2026 (with a Lean certification);
exact quasisymmetry is what remains open; experts expect the answer is no; no reactor design
depends on it. Agreed as a bounded bet.
Evidence: `reports/Stellarator quasisymmetry Lean viability.md`.

**Quasisymmetry residual tool built.**
Checks any closed-form field for equilibrium and for two symmetry residuals at 40 digits. The
axisymmetric member of each family is used as a control and reads zero.
Evidence: `stellarator/qs_residuals.py`, `fields.py`.

**Published closed-form families tested: none is quasisymmetric.**

| Family | Source | Result |
|---|---|---|
| Landreman, rotational transform 2 | arXiv:2609.26742 | exact equilibrium, not quasisymmetric |
| Landreman, sheared | arXiv:2609.26742 | exact equilibrium, not quasisymmetric at the settings tried |
| Issan et al., transform 2, broken stellarator symmetry | arXiv:2610.07304 | not quasisymmetric anywhere in its two parameters |
| Issan et al., sheared, broken stellarator symmetry | arXiv:2610.07304 | not quasisymmetric at the settings tried |
| Lopez vacuum field | arXiv:2610.29921 | curl-free and divergence-free, not quasisymmetric |
| Issan et al., N field periods | arXiv:2611.07672 | exact equilibrium, not quasisymmetric in 7 cases |
| Sengupta et al., force-free | arXiv:2611.31603 | not testable: formulas deferred to an unpublished paper |

Trust: numerical at 40 digits, not proofs. The N-field-period cases were sampled within 3% of the
magnetic axis.
Evidence: `stellarator/qs_residuals_output.txt`, `qs_other_families_output.txt`,
`scan_iota2_output.txt`.

**Seed search by exact algebra.**
Finds every polynomial axisymmetric seed (degree 2 and 3) whose stretched versions are 3D
equilibria. Degree 2 gives a four-coefficient family containing Landreman's seed; both degrees
only give field-line frequency ratio 2. The extra coefficient it found turned out to be the
symmetry-breaking parameter of arXiv:2610.07304, found independently.
Evidence: `stellarator/seeds.py`, `seeds_degree2.json`, `seeds_degree3.json`.

**Overnight search, `deg2` family, on the Mac: no candidates in 531 restarts.**
Best residual about 1.3e-4.
Evidence: `stellarator/runs/deg2.jsonl`.
- Correction: those small residuals are an artifact. `deg2` is the two-parameter Issan family
  with four redundant scaling numbers, and the optimiser used them to shrink the sampled region
  to a sliver of a very large torus, where any field looks nearly symmetric. The direct scan of
  the two real parameters (`scan_iota2.py`) gives a residual of 0.23 or more whenever eps >= 0.2.
  Lesson recorded in the README under known traps: a search family must have a fixed size.

**Overnight search, four-parameter sheared family, on the Mac (run labelled "sengupta").**
About 40 restarts pasted into the conversation; residuals 0.019 to 0.031 at every asymmetry
floor; no candidates. Family not confirmed and results not pushed.

**Lean: statement layer, proof generator, first certificate.**
- `Stellarator/Defs.lean`: gradient, divergence, curl, equilibrium, two forms of quasisymmetry.
  First draft; unreviewed; no definition tested against an example yet.
- `stellarator/leangen.py`: turns SymPy formulas into a complete Lean proof. Handles + - * /,
  powers, named rational sub-expressions and one square root. Refuses false identities.
- `Stellarator/Generated/LandremanIota2.lean`: Landreman's transform-2 field is divergence-free
  and in force balance, all three components, standard axioms only. Stated with one-variable
  partial derivatives; not yet linked to `Defs.lean`.

**Local Taylor-series test (recipe type 3).**
Question: near a point, with no torus, do force balance plus quasisymmetry force a Euclidean
symmetry? Method: Taylor coefficients of B, u, psi, p as unknowns; all equations imposed up to
order K; solutions found numerically and carried up one order at a time.
- Non-symmetric solutions exist at every order from 2 to 11. The symmetry measure reads 1e-15 on
  solutions required to be symmetric and 0.2 to 1 on these.
- Free parameters in the Taylor polynomial of B of degree <= k (two independent solutions agree):

  | | k=1 | k=2 | k=3 | k=4 | k=5 |
  |---|---|---|---|---|---|
  | force balance only | 7 | 19 | 35 | 55 | |
  | force balance + quasisymmetry | 7 | 18 | 33 | 49 | 68, 65, 63, 61, 59, 57, 54 at K = 5 .. 11 |
  | same, with u a rigid motion | 7 | 12 | 16 | 20 | 24 |

- New conditions are still arriving at order 11, so this does not yet show whether a
  non-symmetric family survives.
- Correction 1: an early summary said "no local obstruction". Too strong: constraints were still
  accumulating.
- Correction 2: a run to order 13 seemed to show a collapse toward the symmetric counts. That
  was a numerical artifact of a fixed cutoff in `local_jets_count.py`; the separating singular
  values shrink about tenfold per order.
- Correction 3: `local_jets_count.py` reported the k=4 count falling to 47 or 48 at order 11.
  The stable method (`local_jets_project.py`, top-down elimination) gives 49. Use the stable
  method from now on.
Trust: numerical, truncated series, double precision. Reliable while `min log10 gap` is above
about 3 (it is 5.2 at order 11).
Evidence: `stellarator/local_jets_project_output.txt`, `local_jets_count_output.txt`.

**Local test, orders 12 and 13: the freedom collapses.**
Stable method, same solution carried upward:

| K | k=3 | k=4 | k=5 | k=6 | k=7 | symmetry measure of the solution | min log10 gap |
|---|---|---|---|---|---|---|---|
| 10 | 33 | 49 | 57 | 62 | 71 | 0.22 | 6.4 |
| 11 | 33 | 49 | 54 | 56 | 60 | 0.20 | 5.2 |
| 12 | 33 | 47 | 47 | 47 | 49 | 0.16 | 4.4 |
| 13 | 33 | 33 | 33 | 33 | 33 | 0.07 | 3.7 |

- At order 13 the count is flat from degree 3 to degree 7: every coefficient of B above degree 3
  is fixed by the lower ones. That is a finite-parameter family, and it lost 14 parameters
  between orders 12 and 13.
- The carried solution's symmetry measure is falling at the same time (0.22 to 0.07).
- Reading: the non-symmetric solutions seen at low order look like an effect of truncation that
  the higher-order equations are removing. Whether anything non-symmetric survives is not yet
  known; the trend points toward the equations forcing symmetry near a point.
- Caveats: one carried solution (seed 0) at orders 12 and 13; the gap at order 13 (3.7) is close
  to the limit of what double precision can resolve; no proof.
In progress: orders 14 and 15 in the cloud session, with checkpoints in `stellarator/runs/` so
any machine can resume after a `git pull`. Order 16 on the Mac.

**Not completed today.**
- Calibration run (lost in a restart).
- Sengupta vacuum family.
- Recipe type 2.
