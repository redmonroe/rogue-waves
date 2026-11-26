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

1. **Strong form: from evidence to proof.** State of play: seven numerical solutions all end
   symmetric by order 9 (checked to 40 digits where polished), and exact arithmetic shows
   nothing branches off a Solov'ev equilibrium to first order. Next, in order:
   (a) finish the eight-seed extended-precision batch, and solve from scratch at order 8 and 9
   (`local_strong_hp.py --scratch 8`);
   (b) first-order count at other symmetric equilibria: helical symmetry, non-polynomial flux
   function, other points (does it hold everywhere on the symmetric family?);
   (c) the same first-order question over the rationals with certificates (a non-zero minor mod
   p proves a rank lower bound; explicit kernel vectors prove the upper bound), then Lean for
   the statement at one equilibrium;
   (d) literature: is the local first-order statement known? Read Plunk and Helander 2018 and
   Plunk 2020 (finite pressure) in full.
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

- Is the weak-quasisymmetry property of the planar equilibria stated anywhere? Not found; the
  1995 and 1997 papers could not be read.
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
Trust: the equations are exact (SymPy and Lean).

**Literature check on the planar family.**
Read in full text on arXiv unless noted.
- The twist family, as an equilibrium, is known: it is the constant-ellipticity case of Salat and
  Kaiser, Phys. Plasmas 2, 3777 (1995) ("field lines lie in planes orthogonal to the axis ...
  elliptical cross sections, whose ellipticity and orientation are arbitrary functions along
  the axis"; abstract only, paper not readable from here). It also appears in Gomez-Serrano,
  Liehr and Taylor, arXiv:2609.24739, section 2.3.
- The shift family was not found, but it is of the same known plane-stacked type (Kaiser and
  Salat 1997; Low 1982).
- That these equilibria are weakly quasisymmetric with a non-rigid u was not found stated
  anywhere. Related: Jacobson, arXiv:2408.05919, has non-rigid weak quasisymmetries of a
  symmetric field; Sato, arXiv:2204.09241, has non-symmetric weakly quasisymmetric tori that are
  not equilibria with ordinary pressure.
- Important limit: the family is quasisymmetric only in the weak sense. The strong form adds
  (curl B) x u + grad(u . B) = 0, and here that is proportional to m'' or theta'', not zero.
  Burby, Kallinikos and MacKay (arXiv:2010.08968, Thm 10.1) show weak implies strong for an
  equilibrium with irrational surfaces, so on a real torus the strong form is the one that
  counts. Their open question (arXiv:1912.06468, section 12): "Is every quasi-symmetry a
  Killing field? ... Or at least if one requires magnetohydrostatics? ... It is not even clear
  whether these questions are global or local in nature." That is about the strong form; the
  planar family does not answer it.
- Not read: Salat-Kaiser 1995, Kaiser-Salat 1997, Lortz 1970, Low 1982 (publisher sites
  blocked). Confidence in "not stated before": moderate.

**Local test validated on the exact solution, and what the falling counts mean.**
`local_planar_check.py` puts the Taylor series of the twisted pinch into the local test.
- Every residual of `local_jets.py` is 2.5e-15 at orders 3 to 13. First check of that code
  against a non-symmetric exact solution.
- Free-parameter counts around the exact solution (well-scaled base point):

  | K | k=3 | k=4 | k=5 | k=6 | k=7 | min log10 gap |
  |---|---|---|---|---|---|---|
  | 10 | 33 | 49 | 57 | 62 | 71 | 7.1 |
  | 11 | 33 | 49 | 52 | 56 | 60 | 4.9 |
  | 12 | 33 | 42 | 43 | 45 | 48 | 4.2 |
  | 13 | 27 | 27 | 28 | 29 | 30 | 3.4 |

- Through order 10 these are the same numbers as around the solution carried in the main test.
  They collapse at orders 12 and 13 here too, although this solution exists to all orders.
- So a collapse in the counts does not mean the solution is dying. It means its neighbours are.
  At order 13 the count grows by one per degree, which is what a family with one free function
  (theta) looks like; around the carried solution it is flat at 33 (no free function).
- Correction to the order 12-13 entry above: its "Reading" took the collapse as a sign that
  non-symmetric solutions were a truncation effect. The evidence for that is only the failed
  order-14 step, not the counts.
Evidence: `stellarator/local_planar_check_output.txt`.

**Local test, strong form (started).**
Mode `strong` adds (curl B) x u + grad(u . B) = 0 (imposed one degree lower than the rest;
checked: zero on a symmetric solution, 0.012 on the planar family). First counts, one seed:
7, 16, 24, 29, 33, 39 at order 7, against 7, 18, 33, 49, 63, 80 for the weak form and
7, 12, 16, 20, 24, 28 for symmetric fields. The order-6 to order-7 step already needed 60
Newton steps and the count at order 7 has no clean gap. Running further.

**Local test, strong form: results so far.**
- Independent random solutions get more symmetric with order. Symmetry measure of B over 2 to 8
  solutions per order: 1 to 20 at order 3; 0.3 to 0.7 at order 4; 0.14 to 0.34 at order 5; 0.13 to
  0.15 at order 6.
- Seed 0 was a badly scaled solution (largest coefficient 1.6e4 at order 7); dropped.
- Seed 1, picked at order 4 as the best-scaled of 39 non-symmetric solutions, carried upward:

  | K | freedom in B, k = 1.. | symmetry of B | stretch in u | largest coefficient | min log10 gap |
  |---|---|---|---|---|---|
  | 4 | 7 16 24 | 0.36 | 0.35 | 1.0 | 10.2 |
  | 5 | 7 16 24 30 | 0.34 | 0.23 | 1.0 | 8.8 |
  | 6 | 7 16 24 29 33 | 0.25 | 0.14 | 1.0 | 8.7 |
  | 7 | 7 16 24 29 33 42 | 0.07 | 0.04 | 2.4e2 | 1.3 (not reliable) |
  | 8 | step not converged | 1e-4 | 4e-5 | 2.4e3 | |

- Reading: under the strong form the solution is pushed onto the symmetric family by order 8,
  with its coefficients running away. This is the behaviour the weak form showed at order 14,
  six orders earlier. It is what one would see if force balance plus strong quasisymmetry forced
  a rigid symmetry near a point. One solution; double precision; no proof.
- A first-order version of the same question at a symmetric solution (does any non-symmetric
  change survive?) could not be read: symmetric solutions found by Newton iteration are badly
  scaled (largest coefficient 7e3 at order 6) and the rank decisions fail from order 5.
Evidence: `stellarator/local_strong_ladder_output.txt`.

**Local test, strong form, six more solutions: all end symmetric by order 9.**
`local_strong_seeds.py`: for each seed, the best-scaled non-symmetric solution at order 4,
carried upward. Symmetry measure of B at each order:

| seed | K=4 | K=5 | K=6 | K=7 | K=8 | K=9 |
|---|---|---|---|---|---|---|
| 10 | 0.27 | 0.22 | 0.74 | 0.88 | 0.09 | 4e-14 |
| 11 | 0.23 | 0.22 | 0.19 | 0.11 | 1e-14 | (stays, step not converged) |
| 12 | 0.26 | 0.16 | 0.13 | 0.32 | 2e-14 | 2e-14 (step not converged) |
| 13 | 0.22 | 0.93 | 1.20 | 0.44 | 0.20 | 4e-13 |
| 14 | 0.36 | 0.30 | 0.15 | 0.05 | 0.007 (step not converged) | |
| 15 | 0.39 | 0.55 | 0.37 | 0.21 | 9e-4 | 1e-14 (step not converged) |

- Five of six reach a symmetry measure at rounding level (1e-13 or below) at order 8 or 9; the
  stretch in u goes to 1e-15 with it. The sixth was at 0.007 and falling when its step limit ran
  out. With seed 1 above, that is seven solutions with the same ending.
- A measure of 1e-14 is not "small": it is zero to the precision of the arithmetic. The
  iteration lands exactly on the symmetric family, not near it.
- Coefficients grow to 1e3 to 1e4 on the way, so these runs are near the limit of double
  precision. Residuals relative to coefficient size stay at 1e-15.
- Fixed along the way: `local_extend.py` used an absolute convergence test, which stopped seeds
  10 and 13 early; now relative to coefficient size.
- Reading: evidence that force balance plus strong quasisymmetry force a rigid symmetry near a
  generic point, with the forcing complete at order 8 or 9. Not a proof. Not checked: solutions
  found from scratch at order 9 (only solutions carried up from order 4), other normalisations,
  higher precision.
Evidence: `stellarator/local_strong_seeds_output.txt`.

**Extended-precision check of the strong-form runs (`local_strong_hp.py`).**
Residual evaluated with mpmath at 40 digits; Newton corrections from the double-precision
Jacobian (mixed-precision refinement, a few passes, seconds per solution). Fast enough to run in
the cloud session; no Mac run needed.
- Seed 11: at order 7 a genuine non-symmetric solution (residual 4e-45, symmetry measure 0.105);
  at order 8 symmetric to 39 digits (2.6e-39; double precision could only say 1e-14).
- Seed 13 at order 9: symmetric to 40 digits (1.0e-40).
- So the non-symmetric solutions at orders up to 7 are real solutions of the truncated equations,
  and the symmetric endings are exact, not rounding.
- Batch of eight new seeds (100 to 107), every entry a solution polished to 37 to 53 digits
  unless marked. Symmetry measure of B:

  | seed | K=4 | K=5 | K=6 | K=7 | K=8 | K=9 |
  |---|---|---|---|---|---|---|
  | 100 | 0.31 | 0.23 | 0.13 | 0.077 | step stalled at 0.018 | |
  | 101 | 0.21 | 0.16 | 0.045 | 0.0092 | step stalled at 1e-4 | |
  | 102 | 0.36 | 0.21 | 0.088 | 0.88 | step stalled at 0.05 | |
  | 103 | 0.23 | 0.12 | 0.097 | 0.025 | 1.8e-42 | 4.0e-45 |
  | 104 | 0.22 | 0.82 | 0.80 | 0.83 | 0.0033 | step stalled at 3e-6 |
  | 105 | 0.42 | 0.28 | 0.082 | 0.021 | 2.9e-40 | 6.4e-44 |
  | 106 | 0.30 | 0.16 | 0.14 | 0.13 | 0.0029 | 3.5e-42 |
  | 107 | 0.27 | 0.11 | 0.088 | 0.046 | 0.036 | 4.4e-43 |

- Order 8 does not force symmetry: seeds 104, 106 and 107 are non-symmetric solutions there,
  verified to 37 to 44 digits. Order 9 does in every run that got there: four of four are
  symmetric to 42 digits or better.
- Four runs stalled in the double-precision step (coefficients 6e3 to 1e5). They are neither
  for nor against. Rescaling the length unit (`rescale` option of `local_extend.py`) was tried
  on seed 100 and did not help in the few steps run; left off by default.
- A longer batch for the Mac: `local_strong_hp.py --seeds 200:216 --kmax 10 --steps 400
  --out runs/hp_strong_mac.txt`.
Evidence: `stellarator/runs/hp_strong.txt`.

**Exact arithmetic: nothing non-symmetric branches off an axisymmetric equilibrium, to first
order, under the strong form (`local_modp.py`).**
Setup: a Solov'ev equilibrium (polynomial flux function, pressure linear in psi, F^2 linear in
psi) with random parameters, expanded about a random point; u is the rotation. All arithmetic
with integers modulo a prime: the Taylor coefficients, the residuals (zero mod p for all three
systems: checked), the Jacobian, and Gaussian elimination. No rounding, no cutoffs.
Count = number of first-order changes of the Taylor polynomial of B of degree <= k that keep the
field a solution (no normalisation, so moving, turning and rescaling the field are included in
every row).

| order K | system | k=1 | k=2 | k=3 | k=4 | k=5 | k=6 | k=7 | k=8 |
|---|---|---|---|---|---|---|---|---|---|
| any | symmetric | 11 | 17 | 21 | 25 | 29 | 33 | 37 | 41 |
| 6 | strong | 11 | 21 | 28 | 32 | 38 | | | |
| 7 | strong | 11 | 19 | 23 | 27 | 31 | 37 | | |
| 8 | strong | 11 | 17 | 21 | 25 | 29 | 33 | 39 | |
| 10, 11, 12 | strong | 11 | 17 | 21 | 25 | 29 | 33 | 37 | 41 |
| 12 | weak | 11 | 23 | 30 | 34 | 38 | 42 | 46 | 51 |

- From order 10 on, the strong count equals the symmetric count at every degree up to 8, and
  stays there at orders 11 and 12.
- Same numbers with a second prime and a second set of random parameters, and in three special
  cases: vacuum, force-free (no pressure), and constant F.
- The weak form keeps extra first-order changes (9 more at order 12, possibly still falling),
  as expected from the planar family.
- What this shows: at these equilibria, every first-order change that keeps force balance and
  strong quasisymmetry is itself symmetric, through degree 8. No non-symmetric family can leave
  these equilibria smoothly.
- What it does not show: (1) it is first order only, and only at symmetric equilibria of this
  one type (axisymmetric, Solov'ev); (2) non-symmetric solutions far from any symmetric one are
  not excluded by it (the seven-seed runs above are the only evidence there); (3) it is exact
  arithmetic on random instances mod p, not a proof over the rationals: a proof needs rank
  certificates.
- Known related result: Plunk and Helander, J. Plasma Phys. 84 (2018), show by expansion that
  non-axisymmetric, globally quasi-axisymmetric vacuum fields close to axisymmetry do not exist
  (they exist on a single surface only). Ours is the same conclusion at a point instead of on a
  torus, with pressure and current, in exact arithmetic. Whether the local version is in the
  literature: not checked.
Evidence: `stellarator/local_modp_output.txt`.

**The strong system, reduced by hand (not machine-checked).**
With force balance, strong quasisymmetry, and u . grad p = 0: curl B = -C'(psi) B - p'(psi) u
with C = u . B. So u is determined by B, and the whole system is
  div B = 0,   (curl B) x B = p'(psi) grad psi,   (curl B) . B = -C'(psi) |B|^2 - p'(psi) C(psi):
five equations for B and psi, plus two free functions of one variable. This is the standard
form in which quasisymmetric equilibria are studied (the near-axis expansions of Garren and
Boozer use it), so the strong local test is the same question asked at a generic point instead
of on the axis. Weak quasisymmetry allows the coefficient of B to vary across field lines within
a surface; that extra room is what the planar family uses.

**Recipe type 2 in coordinates (derived by hand, not machine-checked).**
A field with div B = 0, div u = 0, B x u = grad psi and |B| constant along u is the same thing
as a map x(psi, a, b) with unit Jacobian and |x_b|^2 independent of a; then B = x_b, u = x_a.
Force balance adds: x_a . x_b independent of b, and d/db (x_psi . x_b) - d/dpsi |x_b|^2 a
function of psi only. Three unknown functions, four conditions.

**Not completed today.**
- Calibration run (lost in a restart).
- Sengupta vacuum family.
- Local test at order 14: step run to its limit, not converged (see above).
