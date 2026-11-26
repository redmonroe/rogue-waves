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

1. **Check the literature on the planar family** (log entry "Exact quasisymmetric equilibria
   with no rigid symmetry"). Read Sato, Sci. Rep. 2022 (weakly quasisymmetric fields in
   asymmetric tori, no rotational transform, anisotropic pressure) and the infinite-cylinder
   constructions cited by Constantin, Drivas and Ginsberg (arXiv:2009.08860). Is our family
   already written down?
2. **Recipe type 2, continued: look for a version with rotational transform.** The planar
   family works because the symmetry flow is a rigid motion on each plane and B lies in the
   planes. Try: symmetry fields of the form u = K0 + f K1 + ... (K rigid motions, f functions),
   with B not confined to the leaves. Condition to satisfy: B . (grad u + grad u^T) . B = 0.
   Two toroidal variants are already ruled out (see the log entry).
3. **Lean.** Prove the twisted pinch has no rigid-motion symmetry (finite linear algebra at a
   few rational points). Link both certificates to the operators in `Stellarator/Defs.lean` and
   state "IsMHSOn and IsWeakQSOn and not symmetric" in one theorem. This also tests the
   definitions on a real example.
4. **Local test.** The order-14 step is heading into the symmetric family (log). To make that
   solid: repeat from a second seed, and start the test on the planar family (a known
   non-symmetric solution) to confirm the method keeps a real branch alive. Lower priority: the
   planar family shows this test cannot prove non-existence by itself.
5. **Remaining published families.** Sengupta et al. vacuum family (arXiv:2610.28466) not yet
   read. Their force-free family (arXiv:2611.31603) cannot be tested until the companion paper
   with its formulas appears.
6. **Calibration run.** Give Landreman's original prompt to a fresh session, closed book, and see
   whether it rediscovers his field. Started once, lost in a machine restart.
7. **Outside reader** for the Lean statement file and for the planar family.

## Open questions

- Is the planar family new? Not known. One search found related but different results (Sato
  2022: toroidal, anisotropic pressure; ours: ordinary pressure, not toroidal).
- Can the planar construction be given rotational transform, or bent into a torus?
- The local test's carried solution sits on a branch with finitely many parameters (33 at order
  13). Is that branch real, or does it die at higher order?
- Is the local question (force balance + quasisymmetry near a point, no torus) treated in the
  literature? Klotz, Duignan, Burby and Meiss (arXiv:2609.09517) treat the case without force
  balance.
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
- Correction: "the trend points toward the equations forcing symmetry near a point" is wrong as
  a general statement. See "Exact quasisymmetric equilibria with no rigid symmetry" below: a
  closed-form non-symmetric family exists. The collapse seen here is on a different branch.
Checkpoints for orders 3 to 13 are in `stellarator/runs/`, so any machine can resume.

**Local test, order 14: first attempt ran out of time.**
The two-hour run spent 100 minutes trying to carry the order-13 solution to order 14 and printed
nothing. Rewritten as `local_extend.py`: one log line per Newton step, iterate saved each step,
Jacobian on all cores (45 s per step when the machine is otherwise idle).
- The step is converging, slowly: residual 1.2 to 2e-3 in 19 steps, with shortened steps.
- Along the way the solution's symmetry measure keeps falling: 0.071 (order 13) to 0.016.
- Result after 80 steps (the limit): not converged. Residual 4.4e-5 and still falling about 4%
  per step; symmetry measure 0.0030, down from 0.071; steps shortened to 6% of full length.
- Reading: at order 14 the solution is being pulled onto the symmetric family, and the slow,
  shortened steps are what Newton's method does near a point where two solution sets meet.
  So the branch followed since order 3 appears to end in the symmetric family at order 14.
  Not proved: one solution, one seed, not converged.
- Evidence: `stellarator/local_extend_output.txt`. Partial iterate saved as
  `runs/local_qs_seed0_K14_partial.npy`; `local_extend.py 13` resumes from it.
- Two other solvers were tried for this step and were worse: Levenberg-Marquardt damping, and
  picking the greediest of several truncated steps. Both crawl along a narrow valley.

**Local test: two side checks, both inconclusive.**
- Fields whose strength is constant on each surface (isodynamic) are quasisymmetric near a point
  for a trivial reason. The test's normalisation (field strength varies along the field line at
  the origin) excludes them, so they were never in the counts. Added as mode `iso`: counts 5, 12,
  18, 28 at order 5, falling with order; the solver stalls at order 6.
- `local_rigidity.py`: count first-order quasisymmetric changes of a symmetric solution. If the
  count equalled the symmetric count, nothing non-symmetric could branch off. Result: larger
  through order 6 (48 against 20 at degree 4); beyond that the rank decisions fail (gap near 1).

**Exact quasisymmetric equilibria with no rigid symmetry (recipe type 2, first result).**
Found by hand while setting up "symmetry by construction", then checked by machine.
Construction: take a 2D field in the (x, z) plane whose field-line force is the gradient of a
function unchanged by shifts in x (or by rotations about the origin). In the plane at height y,
use that field shifted by m(y) (or turned by theta(y)), for any function m or theta. Then:
- force balance holds, with pressure a function of the flux label;
- u = "go up in y and apply the motion" is divergence-free, B x u = grad psi, and |B| is
  constant along u: weak quasisymmetry, exactly;
- if m (or theta) is not linear in y, u is not a rigid motion and B has no translation,
  rotation or screw symmetry.

| Family | Field | Surfaces |
|---|---|---|
| shift | B = exp(kappa z/2) (cosh(mu X), 0, -(2 mu/kappa) sinh(mu X)), X = x - m(y) | open sheets |
| twist | B = R(y) M R(y)^T (x, z), M = [[0, -2 beta], [2 alpha, 0]], R(y) = rotation by theta(y) | nested elliptical tubes around the y axis, turning at any rate |

Checks:
- All five equations vanish identically in SymPy for both families (`planar_family.py`).
- Rigid-symmetry measure at 40 digits: 1e-41 for the controls (m or theta linear, or circular
  cross-section), 0.007 to 0.22 for the others (`planar_family_output.txt`).
- Lean: for the twist family with tan(theta/2) = y the field is rational, and
  `Stellarator/Generated/TwistedPinch.lean` proves div B = 0, force balance (3 components),
  div u = 0, B x u = grad psi (3 components) and u . grad |B|^2 = 0. Standard axioms only.
  Non-symmetry is not proved in Lean yet. Statements use one-variable partial derivatives and
  are not yet linked to `Defs.lean`.

What it is not: an answer to the open problem. Every field line lies in a plane y = const (no
rotational transform), and the surfaces are sheets or infinite tubes, not tori in ordinary
space. The twist family is a pinch with no field along its axis; making y periodic gives tori
only in a periodic box.
What it changes here: force balance plus quasisymmetry do not force a rigid symmetry near a
point, so no purely local argument can prove non-existence; a proof would have to use the torus
or the rotational transform.
Two toroidal variants already ruled out by hand: (a) poloidal field shifted vertically by
m(toroidal angle) needs the vertical field component constant along closed field lines,
impossible; (b) turning the pattern within each poloidal plane is not volume-preserving.
Trust: the equations are exact (SymPy and Lean). Novelty unknown; not checked against the
literature beyond one search.

**Recipe type 2 in coordinates (derived by hand, not machine-checked).**
A field with div B = 0, div u = 0, B x u = grad psi and |B| constant along u is the same thing
as a map x(psi, a, b) with unit Jacobian and |x_b|^2 independent of a; then B = x_b, u = x_a.
Force balance adds: x_a . x_b independent of b, and d/db (x_psi . x_b) - d/dpsi |x_b|^2 a
function of psi only. Three unknown functions, four conditions.

**Not completed today.**
- Calibration run (lost in a restart).
- Sengupta vacuum family.
- Local test at order 14: step run to its limit, not converged (see above).
