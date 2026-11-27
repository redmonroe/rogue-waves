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

1. **Strong form with force balance ("bottom-right cell"): from evidence to proof.** No proof
   yet. The claim is a local form of the Garren-Boozer conjecture (1991) and is open (log,
   "Bottom-right cell"). In order:
   (a) in progress: sixteen-seed batch on the Mac (`runs/hp_strong_mac.txt`);
   (b) smallest open sub-case, to try by hand: vacuum with u . B = 0. Statement: Phi harmonic,
   psi constant along grad Phi, and |grad Phi| a function of (psi, Phi) only. Show Phi has a
   rigid symmetry. Then vacuum with u . B = C not zero;
   (c) first-order question over the rationals with certificates, then Lean, at one
   equilibrium. Note the first-order test fails at u . B = 0 (two extra directions, removed
   only at second order), so a proof by linearisation alone cannot cover that case;
   (d) a note to the people working on this (Burby, MacKay, Sengupta and co-authors) with the
   two exact families and the order-9 evidence: they are attacking the same question with
   Cartan-Kahler methods.
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
7. **Outside reader** for the Lean statement file, the planar family and the radial twisted
   field. The last one is the easiest to check and bears on a question stated as open in a 2025
   paper; a one-page note to its authors is the natural route.

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
- Are the radial twisted field and its relatives known? They are elementary. Four searches
  found no explicit strongly quasisymmetric field without a rigid symmetry, and the 2025 paper
  says existence is open, but the authors may mean toroidal fields only.
- Is there a strongly quasisymmetric field without force balance and without rigid symmetry
  whose field lines are not planar? Both families found here have planar field lines.

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
- Two more types of symmetric equilibrium give the same result (strong count = symmetric count
  at every degree up to 8 from order 10; checked to order 11):
  - translation symmetry, flux function a sum of four exponentials (non-polynomial), pressure
    and axial field non-linear in the flux;
  - screw symmetry: the elliptical pinch twisting at a uniform rate (the symmetric member of
    the twisted-pinch family). Looked for a screw-symmetric equilibrium with an axial field and
    quadratic flux function: SymPy finds none, only the pinch.
  The weak counts differ slightly between the three types at orders 10 and 11; the strong and
  symmetric counts are identical in all of them.
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

**Strong form from random starts at order 9: no solution found, so no information.**
Six random starts (seeds 300 to 305), 300 Newton steps each, no ladder from order 4.
- None converged: residuals 0.3 to 2.1 after 300 steps, step lengths 0.00 to 0.03.
- So this says nothing for or against: the solver did not reach any solution at order 9,
  symmetric or not. At order 6 the same procedure reached non-symmetric solutions (0.98, 0.78)
  from two of three starts.
- The gap stays open: every order-9 solution we have was carried up from order 4.
Evidence: `stellarator/runs/hp_strong_scratch9.txt`.

**Literature: what is known about "does quasisymmetry force a rigid symmetry".**
Read in full or in the relevant parts: Plunk and Helander 2018 (arXiv:1801.02990), Plunk 2020
(arXiv:2005.02981), Burby, Kallinikos and MacKay 2020 (arXiv:1912.06468), Burby, Kallinikos,
MacKay, Perrella and Pfefferle 2024 (arXiv:2403.03352, J. Plasma Phys. 2025), Sato, Qu,
Pfefferle and Dewar 2021 (arXiv:2108.01818), Klotz et al. 2026 (arXiv:2609.09517, abstract and
citations only).
- Vacuum, near axisymmetry, first order: done by Plunk and Helander (Appendix F). Their argument
  is pointwise in the poloidal plane, one toroidal mode number at a time, so the vacuum case of
  our exact first-order result is in effect already there.
- Finite pressure, near axisymmetry: Plunk 2020 imposes the condition on one surface and says
  "we will not consider here the question of whether this condition might, under special
  circumstances, be satisfied globally". So our first-order result with pressure and current is
  not in that paper. Not found elsewhere either.
- Strong form, general: Burby, Kallinikos and MacKay 2020, section 7: "We have not managed to
  prove or disprove that a quasi-symmetry is a Killing field". Section 12 asks it again, "at
  least if one requires magnetohydrostatics", and whether it is local or global.
- Strong form without force balance: Burby et al. 2024 prove that if u is not a rigid motion,
  the strain S of u has zero trace and zero determinant, B points along the null direction of
  S, and that direction is carried by the flow of u. They say the existence of such u "needs to
  be addressed" and announce a follow-up. Klotz et al. (September 2026) still cite it as a
  proposal.
- Weak form without force balance: non-symmetric local solutions are known (Sato et al. 2021,
  Sato 2022; Klotz et al. 2026).
- Limits of this check: four web searches and six papers. One 2025 paper could not be opened
  (osti.gov is not reachable from here); from its abstract (arXiv:2507.06480) it is about |B| on
  one surface. Not exhaustive.

**New exact fields: strongly quasisymmetric, no force balance, no rigid symmetry
(`radial_twist.py`).**
Found while reading Burby et al. 2024: the symmetry field of our twisted pinch satisfies their
three conditions, so their theorem says which B goes with it.
- Field: B = G(phi - theta(y)) r_hat / r, pointing straight out from the y axis; phi is the
  angle about that axis; theta and G are any functions. Symmetry field: u = d/dy + theta'(y)
  d/dphi (move along the axis while turning about it).
- SymPy, for arbitrary theta and G: div u = 0, div B = 0, curl(u x B) = 0,
  (curl B) x u + grad(u . B) = 0, u . grad|B|^2 = 0, all identically zero. The strain of u has
  S:S = r^2 theta''^2 / 2, so u is not a rigid motion where theta'' is not zero.
- B itself has no rigid symmetry when theta'' and G' are not zero. By hand: a rigid motion
  preserving the direction of B preserves the y axis, so it is a d/dy + b d/dphi, which changes
  B by G' (b - a theta') / r. Numerically at 40 digits: measure 0.016 for tan(theta/2) = y,
  G = cos^2; 1e-41 for the two controls (theta linear; G constant).
- Lean, for tan(theta/2) = y and G = cos^2 (a rational field): `Stellarator/Generated/RadialTwist.lean` proves div B = 0, div u = 0, curl(u x B) = 0 (3 components), u . B = 0, (curl B) x u + grad(u . B) = 0 (3 components), u . grad|B|^2 = 0, and d(u_x)/dy = -4yz/(1+y^2)^2. Hypothesis: x^2 + z^2 is not zero. Standard axioms only; `lake build` passes. Non-symmetry of B is not proved in Lean. Not linked to `Defs.lean`.
- Second family by the same recipe (u = K0 + f K1, with K0 and K1 commuting rigid motions and f
  constant along K1): u = d/dx + a(y) exp(kappa x) d/dz, B = Phi(a(y) exp(kappa x) - kappa z)
  (-a'/a, kappa, 0), any functions a and Phi. SymPy: all equations zero; measure 0.034 for
  a = y, 1e-41 for the control a = exp(y).
- What they are not: equilibria (the current is parallel to u and J x B is not a gradient
  unless G G' = 0); toroidal (field lines are straight rays, or planar curves in the second
  family); an answer to the main problem.
- What they show: in flat space, strong quasisymmetry by itself does not force a rigid
  symmetry of B, even locally. With the planar family (weak form, with force balance) the
  picture is now:

  | | weak quasisymmetry | strong quasisymmetry |
  |---|---|---|
  | no force balance | non-symmetric fields known (Sato; Klotz et al.) | non-symmetric fields exist (this entry) |
  | force balance | non-symmetric fields exist (planar family, ours) | none found; evidence that none exist with \|B\| varying along B |

  So the rigidity seen in the strong-form runs needs both ingredients, force balance and the
  strong condition. Each alone leaves room.
- Trust: the equations are exact (SymPy for the general case, Lean for one case). Novelty: not
  known. The examples are elementary and specialists may know them or regard them as trivial;
  see the open questions.
Evidence: `stellarator/radial_twist_output.txt`, `RogueWaves/Stellarator/Generated/RadialTwist.lean`.

**Bottom-right cell (force balance + strong quasisymmetry): what is to be proved, a smaller
test, and where a proof stands.**
Checked first that nothing published contradicts the claim: Sengupta et al. 2021
(arXiv:2011.08149) have exact quasisymmetry on one surface only; Sengupta et al. 2026
(arXiv:2610.28466) have exact vacuum fields with nested surfaces and no symmetry, but not
quasisymmetric (they say so).

Reductions, by hand (not machine-checked unless stated):
- Strong quasisymmetry gives S . B = 0 and trace S = 0 (S = strain of u; Burby et al. 2024). So
  u is a rigid motion as soon as S . V = 0 for one vector V not along B. Equivalent: grad psi,
  or grad |B|, is carried along by the flow of u.
- With force balance, curl B = -C'(psi) B - p'(psi) u (logged before). Without it, u . B need
  not be constant on flux surfaces.
- Vacuum: B = grad Phi, and the whole system is
      Laplace(Phi) = 0,   grad psi . grad Phi = 0,   Hess Phi (u, B) = 0,
      u = (C grad Phi + grad psi x grad Phi) / |grad Phi|^2,   C a constant.
  Two unknown functions and a constant. Machine check: the residuals vanish on an exact
  translation-symmetric solution (`local_vacuum.py --check`).
- Counting at order K in vacuum (free coefficients against equations): 116 / 84 at order 7,
  144 / 120 at 8, 175 / 165 at 9, 209 / 220 at 10. So the system runs out of freedom at orders
  9 to 10, which is where the runs turn symmetric. This is Garren and Boozer's counting
  argument, at an ordinary point instead of on the magnetic axis. Symmetric solutions survive
  because for them the equations are in two variables. The open question is whether anything
  else survives by accident.

Vacuum test (`local_vacuum.py`; about 500 unknowns at order 9 against 1700 in the general test):
- From random starts: orders 6 and 7, 13 of 24 converge, all non-symmetric (asymmetry 0.02 to
  2.5). Orders 8 and 9, 0 of 24 converge.
- Ladders from order 6 (8 of 16 starts converged there): at order 8 the six that got there are
  still non-symmetric (0.03 to 0.26). At order 9 five converge and all five are symmetric; one
  did not converge (asymmetry 0.06 and falling), one diverged. They stay symmetric at orders
  10 and 11.
- Two of the five have u . B not zero and u a rigid motion (strain 1e-12 or less).
- Three have u . B = 0. There u is a function of psi times a rigid motion, so its strain is not
  zero; that is the known relabelling freedom, not a new solution.
- Two of those three stopped at asymmetry 1e-7, not 1e-13. Checked at 50 digits
  (`local_vacuum_c0.py`): the asymmetry halves at every Newton step and the residual falls by
  four, down to 4e-10 where the method runs out of resolution. So the exact solution is
  symmetric, and it is a double root: 1e-7 is the square root of the rounding level.
- Consequence: at symmetric solutions with u . B = 0 the first-order test is not rigid (two
  extra first-order directions, removed at second order). The exact first-order results logged
  earlier were all at u . B not zero.
Evidence: `stellarator/runs/vacuum_scratch.txt`, `runs/vacuum_ladder.txt`,
`runs/vacuum_c0_check.txt`.

Where a proof stands: not found.
- By elimination: the conditions on the lowest coefficients only appear after all coefficients
  up to order 8 or 9 have been eliminated (384 unknowns in vacuum). Exact elimination at that
  size is out of reach by brute force, and the modular method needs exact points on the
  solution set, which we only have for symmetric solutions.
- By an identity: none found. Tried: the invariant forms d psi, d|B|, B-flat and what force
  balance adds to them (everything first-order is already used; the constraint is in the
  second-order, flatness, conditions); ellipsoidal equipotentials as a possible counterexample
  (fails, as expected).
- By hand, unchecked: the problem is equivalent to classifying flat metrics whose inverse, in
  the coframe (d psi, d|B|, B-flat), has the form [[X, Y, 0], [Y, Z, f], [0, f, |B|^2]] with f a
  function of (psi, |B|) and fixed determinant; the claim is that X and Y cannot depend on the
  third coordinate. Two unknown functions, three independent flatness equations.

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
