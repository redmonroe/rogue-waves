# Lean 4 + Mathlib + Claude + one non-mathematician: fit for stellarator mathematics (quasisymmetry, MHS equilibria, near-axis expansion), as of 26 November 2026

Evidence labels used throughout:
- **[FH-local]** first-hand: read or executed on this machine (Mathlib checkout at commit `31ce06ea`, dated 19 Nov 2026, toolchain `v4.37.0-rc1`; scratch Lean files compiled with `lake env lean`).
- **[FH-repo]** first-hand: public GitHub repository cloned in this session and files read.
- **[FH-paper]** first-hand: arXiv PDF or abstract page downloaded and read in this session.
- **[snippet]** from a web-search result summary only; the underlying page was NOT opened (most hosts, including quantamagazine.org, en.wikipedia.org, scientificamerican.com, leanprover-community.github.io, are blocked in this environment).
- **[relayed]** told to me by the coordinator, not seen by me.
- **[estimate]** my own judgement, not a sourced fact.

Access note: WebFetch worked only for anthropic.com; arXiv, `git clone`, and raw.githubusercontent.com worked from Bash. Everything about Quanta, Scientific American, Wikipedia, Hacker News and Lean Zulip is [snippet] unless stated.

Headline corrections to the brief that the report writer must not miss:
1. **Grad's conjecture is no longer open in its standard formulation.** It was disproved in September 2026, and one of the disproofs ships with a 567,000-line Lean 4 certification. The stellarator-specific open problem that remains is the existence of exact non-axisymmetric quasisymmetric equilibria.
2. **Magnetohydrostatics IS now formalized in Lean** (that same certification), and plasma kinetic theory is too (Vlasov-Maxwell-Landau; relativistic Vlasov-Maxwell). The claim "nothing on magnetic confinement is formalized anywhere" was true until 2026 and is now false.
3. **The 2026 pattern is not "extend Mathlib"**. Every large AI-era formalization I inspected states its theorem in a short Mathlib-only file with hand-rolled `curl`/`div` built on `fderiv`, and builds missing analysis privately. Missing Mathlib infrastructure is a cost, not a blocker.

---

## 1. Infrastructure audit: what this mathematics needs versus what Mathlib has

### Takeaway
Mathlib (19 Nov 2026) has the calculus substrate (Fréchet derivatives, exterior derivative on normed spaces, Lie bracket of vector fields, Picard-Lindelöf, integral curves on manifolds, formal power series with derivations, Fourier series on tori, Riemannian metrics and a Levi-Civita connection) but has **no named curl/divergence operators, no differential forms on manifolds, no Lie derivative, no Stokes theorem beyond boxes, no symplectic/Hamiltonian mechanics, no KAM, no Cauchy-Kovalevskaya, no bounded-domain Sobolev spaces and no elliptic PDE theory**. Two scratch probes I wrote compiled on the first attempt, so the pointwise vector-calculus layer is cheap to build locally.

### Cited Findings

Audit table. All rows are [FH-local] from grepping `/home/user/rogue-waves/RogueWaves/.lake/packages/mathlib/Mathlib` and from `#check` in a scratch file; paths are relative to that directory.

| Topic | Status | What I found |
|---|---|---|
| Cross product in R^3 | **Present** | `LinearAlgebra/CrossProduct.lean`: `crossProduct` (notation `⨯₃`) on `Fin 3 → R`, with `cross_anticomm`, `triple_product_permutation`, `triple_product_eq_det`, `cross_dot_cross`, `jacobi_cross`, `cross_cross_eq_smul_sub_smul`. Norm lemmas for `EuclideanSpace ℝ (Fin 3)` in `Geometry/Euclidean/Angle/Unoriented/CrossProduct.lean`. |
| Gradient | **Present** | `gradient` in `Analysis/Calculus/Gradient/Basic.lean` (Hilbert-space gradient of a scalar function). |
| Curl, divergence as named operators | **Absent** | `grep '\bcurl\b'` over all of Mathlib hits only `Tactic/Polyrith.lean` (the HTTP tool). Divergence appears only inline as `∑ i, f' x eᵢ i` in the divergence theorem. |
| Vector-calculus identities (div curl = 0, curl grad = 0) | **Absent as lemmas; trivially derivable** | Not in Mathlib. Second-derivative symmetry is there (`ContDiffAt.isSymmSndFDerivAt`). See probe below. |
| Laplacian, harmonic functions | **Present (basic)** | `Analysis/InnerProductSpace/Laplacian.lean`, `Analysis/InnerProductSpace/Harmonic/{Basic,Constructions,HarmonicContOnCl}.lean`. |
| Exterior derivative, differential forms | **Partial: normed spaces only** | `Analysis/Calculus/DifferentialForm/Basic.lean` (266 lines): `extDeriv`, `extDeriv_extDeriv` (d∘d = 0). `.../VectorField.lean`: `extDeriv_apply_vectorField` (the invariant formula with Lie brackets). The module docstring's TODO says forms on manifolds are "not defined yet". Only the bundle of alternating maps exists on manifolds (`Geometry/Manifold/VectorBundle/ContinuousAlternatingMap.lean`). |
| Wedge product, interior product, Hodge star | **Partial (algebraic level)** | `LinearAlgebra/ExteriorPower/{WedgePairing,HodgeStar}.lean`, `Analysis/Normed/Module/Alternating/Curry.lean`. No calculus API (no Cartan formula). |
| Lie derivative (of forms, metrics, functions) | **Absent** | Grep for `LieDeriv|lieDeriv|Lie derivative` finds only Lie-algebra files and one mention in `CovariantDerivative/Metric.lean`. |
| Lie bracket of vector fields | **Present** | `Analysis/Calculus/VectorField.lean` (713 lines): `VectorField.lieBracket`, `pullback_lieBracket`, `leibniz_identity_lieBracket`. Manifold version `VectorField.mlieBracket` in `Geometry/Manifold/VectorField/LieBracket.lean`. |
| Stokes / divergence theorem | **Partial: rectangular boxes only** | `MeasureTheory/Integral/DivergenceTheorem.lean` (`integral_divergence_of_hasFDerivAt_off_countable`), `Analysis/BoxIntegral/DivergenceTheorem.lean`. No Stokes on manifolds or general domains. |
| Poincaré lemma / de Rham | **Partial / absent** | `MeasureTheory/Integral/CurveIntegral/Poincare.lean`: closed 1-forms on convex sets are exact. No de Rham cohomology (the only "de Rham" hit is the p-adic period ring). |
| Smooth manifolds, tangent bundle, vector bundles | **Present** | `Geometry/Manifold/` (charted spaces, `IsManifold`, `ContMDiff`, `mfderiv`, immersions, submersions, Whitney embedding, partitions of unity). |
| Riemannian manifolds, Levi-Civita | **Present (recent)** | `Geometry/Manifold/Riemannian/Basic.lean`, `Geometry/Manifold/VectorBundle/CovariantDerivative/{Basic,Torsion,Metric,LeviCivita}.lean` with the Koszul formula; smoothness of the connection is listed as future work. No Killing fields. |
| ODE existence/uniqueness | **Present** | `Analysis/ODE/{Basic,PicardLindelof,ExistUnique,Gronwall,Transform}.lean`: `IsPicardLindelof`, `IsIntegralCurve(On/At)`, `ODE_solution_unique`, local flow Lipschitz and continuous in the initial point. |
| Smooth dependence of flows on initial data / parameters | **Absent or not found** | Grep for C^k regularity of flows in `Analysis/ODE`, `Geometry/Manifold/IntegralCurve`, `Dynamics` found only "a C¹ vector field admits a local flow" (existence). No differentiable-dependence theorem located. |
| Integral curves and flows on manifolds | **Partial** | `Geometry/Manifold/IntegralCurve/{Basic,ExistUnique,Transform,UniformTime}.lean`: `IsMIntegralCurve`, local existence/uniqueness, uniform-time lemma giving global curves. No flow as a one-parameter group of diffeomorphisms; boundary case is a TODO. |
| Frobenius theorem, foliations | **Absent** | No relevant hits. |
| Tori and periodic functions | **Partial** | `AddCircle`, `UnitAddTorus d` (topological/measure level). Manifold structure exists for `Circle` ⊂ ℂ (`Geometry/Manifold/Instances/Sphere.lean:570`); no `ChartedSpace` instance on `AddCircle` found. |
| Fourier series | **Present** | `Analysis/Fourier/AddCircle.lean` (`fourierCoeff`, `hasSum_fourier_series_of_summable`), `AddCircleMulti.lean` (`UnitAddTorus.mFourier`, `mFourierCoeff`, `mFourierBasis`, `hasSum_mFourier_series_L2`). No Fourier-side regularity (decay of coefficients vs smoothness) found. |
| Rotation number | **Present (circle maps)** | `Dynamics/Circle/RotationNumber/TranslationNumber.lean`. Nothing on rotational transform of flows on tori. |
| Formal power series, multivariate | **Present and rich** | `RingTheory/PowerSeries/` (21 files: `Derivative`, `Substitution`, `Inverse`, `Trunc`, `Order`, `WeierstrassPreparation`), `RingTheory/MvPowerSeries/` (17 files: `MvPowerSeries.pderiv` as a derivation, `subst`, `Trunc`, `Order`). |
| Analytic functions | **Present** | `Analysis/Analytic/` (19 files), analytic inverse function theorem (`InverseFunctionTheorem/Analytic.lean`), implicit function theorem (`Analysis/Calculus/Implicit*.lean`, `ImplicitFunction/`). |
| Cauchy-Kovalevskaya | **Absent** | Zero hits for `Kovalevsk|Kowalevsk`. |
| Sobolev spaces | **Partial: whole space only** | `Analysis/Distribution/Sobolev.lean` (`TemperedDistribution.MemSobolev`, Bessel-potential definition via Fourier multipliers), `Analysis/FunctionalSpaces/{BesselPotentialSpace,SobolevInequality}.lean` (Gagliardo-Nirenberg-Sobolev for compactly supported C¹). No W^{k,p} on bounded domains, no traces. |
| Elliptic PDE theory (Lax-Milgram applications, Schauder, regularity) | **Absent** | Nothing found beyond harmonic-function basics. |
| KAM theory | **Absent** | Zero hits. |
| Symplectic / Hamiltonian mechanics | **Absent** | Only `LinearAlgebra/SymplecticGroup.lean` (matrices) and Lie-algebra classical types; "Hamiltonian" hits are graph theory. |
| Nash-Moser / tame Fréchet calculus | **Absent** | Not in Mathlib (it had to be built privately for the Grad counterexample, see section 3). |

- Scratch probe 1 (mine, compiled first try, about 3 s wall-clock including imports): hand-rolled `pd`, `grad`, `div`, `curl` on `Fin 3 → ℝ`; proved `div (curl B) = 0` and `curl (grad f) = 0` for C² fields in about 45 lines; stated `IsMHS` and a triple-product quasisymmetry predicate. Axioms: `propext`, `Classical.choice`, `Quot.sound`. A deliberately false `example` in the same file was rejected, confirming the check is live. [FH-local] — [Probe1.lean](file:///tmp/claude-0/-home-user/f7c08234-634f-5a50-bfc5-78c71b1b878c/scratchpad/lean/Probe1.lean)
- Scratch probe 2 (mine, compiled first try): the axisymmetric seed field of Landreman's September 2026 explicit equilibria (his eq. 2.6, written in Cartesian components with a square root) is divergence-free wherever R² > 0 and I² > 0; about 45 lines using `HasDerivAt` combinators, `field_simp`, `ring`. [FH-local] — [Probe2.lean](file:///tmp/claude-0/-home-user/f7c08234-634f-5a50-bfc5-78c71b1b878c/scratchpad/lean/Probe2.lean)
- **I reproduced a documented failure mode in my own probe.** My `IsMHS` structure wrote `ContDiffOn ℝ ⊤ B U` intending "smooth". In current Mathlib `(⊤ : WithTop ℕ∞) = ω` holds by `rfl` and `ContDiff ℝ ⊤ f` implies `AnalyticOnNhd` (checked: `ContDiff.analyticOnNhd`). So my statement silently demanded real-analyticity. The correct spelling is `∞` under `open scoped ContDiff`. The file compiled without complaint. [FH-local] — [Check2.lean](file:///tmp/claude-0/-home-user/f7c08234-634f-5a50-bfc5-78c71b1b878c/scratchpad/lean/Check2.lean)
- Outside Mathlib, an AI-written downstream library now covers much of the missing list. Tau Ceti (commit of 26 Nov 2026, same toolchain `v4.37.0-rc1` as the local project) has 43,025 Lean files, including `Geometry/Manifold` (2,691 files, with `DifferentialForm/{Integration,Pullback,Wedge,Flow}`, `DeRham/…`, `Riemannian/Curl.lean`, `Sobolev`), `Geometry/Symplectic` (1,037 files: Hamiltonian, Liouville, moment maps, reduction), `Analysis/PDE` (646: Hölder, Schauder, Krylov-Safonov, parabolic), `Analysis/Sobolev` (449: W^{k,p}, trace, extension, Morrey, compactness). Path-name counts: `DifferentialForm` 172, `DeRham` 64, `Hamiltonian` 228, `Stokes` 8, `Curl` 4, `LieDerivative` 1 (Kähler differentials only), `Kovalev` 0, `NashMoser` 0, `NavierStokes` 0. I inspected file names and one file, not quality. [FH-repo] — [TauCetiProject/TauCeti](https://github.com/TauCetiProject/TauCeti)
- Physlib (ex-PhysLean) has named operators: `Space.grad`, `Space.div`, `Space.curl` with `div_of_curl_eq_zero`, `curl_of_grad_eq_zero`, `curl_of_curl`, `exists_curl_of_div_zero`, `exists_grad_of_curl_zero` in `Physlib/SpaceAndTime/Space/Derivatives/Curl.lean`. [FH-repo] — [leanprover-community/physlib](https://github.com/leanprover-community/physlib)

### Inferences
- For anything **pointwise and local** (MHS force balance at a point, the Burby-Kallinikos-MacKay conditions written in R³ vector-calculus form, triple products, checking an explicit field), Mathlib already suffices; the cost is a few hundred lines of local definitions. [estimate]
- For anything **global** (nested toroidal flux surfaces, Boozer coordinates, rotational transform, circle actions, integrals over flux surfaces) Mathlib has no ready vocabulary. The Grad-counterexample team encoded these with explicit parametrizations and quotient types in a 242-line statement file rather than with manifold theory (section 3), which is the model to copy. [estimate]
- For **existence theory** (elliptic regularity, Nash-Moser, KAM) nothing is in Mathlib; Tau Ceti has pieces (Sobolev, Schauder) but not Nash-Moser or KAM. [FH-repo for the path counts; inference for the consequence]

### Gaps
- I did not test whether Tau Ceti builds or whether its PDE/differential-form files are usable from a downstream project; only file names and one file were read.
- Smooth dependence of ODE flows on parameters: my greps found nothing, but absence by grep is weaker than a maintainer's statement; a Mathlib overview page could not be fetched (host blocked).

---

## 2. What the stellarator mathematics needs, target by target

### Takeaway
The four families of statements separate cleanly: (iii) local equivalence of quasisymmetry definitions and (iv) "axisymmetric implies quasisymmetric" are pointwise calculus and are feasible now; (ii) the Garren-Boozer near-axis expansion is formal-power-series algebra whose *derivation* is feasible but whose "counting argument" is not a theorem in the literature; (i) stating the existence question is now easy (a template exists), while proving it is research mathematics.

### Cited Findings
- Burby-Kallinikos-MacKay define a quasisymmetry as a vector field `u` such that `(u,0)` is a continuous symmetry of first-order guiding-centre motion for all magnetic moments (Def. 4.1), and prove it is equivalent to `L_u|B| = 0`, `L_u β = 0`, `L_u b♭ = 0` (Thm 4.2), and to three alternative sets of Lie-derivative conditions (Thm 4.4). They state that "any axisymmetric B in Euclidean space" is quasisymmetric with `u = ∂_φ`. They do not assume MHS. [FH-paper] — [arXiv:1912.06468](https://arxiv.org/abs/1912.06468)
- The same paper's open question: "Is every quasi-symmetry a Killing field? At least in the Euclidean case? Or at least if one requires magnetohydrostatics?" Theorem 7.2 constrains `L_u g` to a codimension-4 subspace of the 6-dimensional space of symmetric matrices ("two-thirds of the way"). They suggest the question may be local and attackable by "a prolongation analysis based on the Cartan-Kuranishi theorem". [FH-paper] — [arXiv:1912.06468](https://arxiv.org/abs/1912.06468)
- Rodríguez-Helander-Bhattacharjee derive necessary and sufficient conditions from the single-particle Lagrangian without assuming MHD equilibrium, state their conditions are "less restrictive" than Burby et al.'s, relate them to the triple-product form, and show quasisymmetry implies existence of flux surfaces. [FH-paper] — [arXiv:2004.11431](https://arxiv.org/abs/2004.11431)
- Landreman-Sengupta (2019): "In Garren and Boozer's original work, it was argued quasisymmetry can be achieved (in the absence of axisymmetry) to O((r/R)²) but not to O((r/R)³)"; at second order "the number of scalar ϕ-dependent unknowns exceeds the number of ϕ-dependent" equations; the O(r¹) equations are nonlinear with a unique solution and the O(r²) equations are linear in the higher-order unknowns. [FH-paper] — [arXiv:1908.10253](https://arxiv.org/abs/1908.10253)
- Constantin-Drivas-Ginsberg's formal statement of Grad's conjecture: "Any non-isolated and non-vanishing (away from the 'magnetic axis') smooth unforced MHS equilibrium on a (topologically toroidal or cylindrical) domain T ⊂ R³ that has a pressure possessing nested level sets which foliate T has either plane-reflection, axial or helical symmetry." They add that addressing it "seems to be out of reach of existing techniques". [FH-paper] — [arXiv:2007.09103](https://arxiv.org/abs/2007.09103)
- That conjecture has since been stated and refuted in Lean: the Grad-counterexample repository's `Showcase.lean` (242 lines, Mathlib imports only) defines `Vec := EuclideanSpace ℝ (Fin 3)`, `gradient`, `cross`, `curl`, `divergence` via `fderiv`, `TangentTo`, `IsMHS`, `IsEmbeddedTorus`, `IsPressureFoliation`, `NestedPressureTori`, `SignedStabilizes`, `ExactCyclicSymmetry`, a moduli space as a `Quot`, and the theorem `equilibria_with_exact_cyclic_symmetry`. [FH-repo] — [lukasliehr/Grad-Conjecture](https://github.com/lukasliehr/Grad-Conjecture)

What each item needs (my analysis; all [estimate] except where a Mathlib fact from section 1 is cited):

| Item | Mathematics needed | In Mathlib? | Character |
|---|---|---|---|
| (i-a) State MHS equilibrium with nested toroidal pressure surfaces | Pointwise curl/div; embedded tori; foliation by level sets | Operators hand-rolled; topology via explicit parametrizations | **Done already** in `Showcase.lean`; reusable |
| (i-b) State exact quasisymmetry on such an equilibrium | Add a vector field `u` with `div u = 0`, `B × u = ∇ψ`, `u·∇|B| = 0` (Euclidean form of the Burby et al. conditions), or the triple-product form; plus "not a Killing field" for non-triviality | Yes, pointwise | Easy to write, **hard to get right** (which definition; where `B ≠ 0`; regularity `∞` vs `ω`; what "non-axisymmetric" means) |
| (ii-a) Near-axis expansion: derive the order-by-order equations | Formal power series in `r` whose coefficients are trigonometric polynomials in θ with coefficients smooth periodic in φ; Frenet frame of the axis; Boozer-coordinate identities | `PowerSeries`/`MvPowerSeries` with derivations: yes. Coefficient ring of smooth periodic functions with a derivation: must be defined. Existence of Boozer coordinates: **no** | Mostly algebra **if Boozer coordinates are taken as a hypothesis/structure**; heavy bookkeeping |
| (ii-b) Solvability at O(r), O(r²) | Periodic solution of a Riccati-type ODE for σ(φ); linear periodic systems | Picard-Lindelöf yes; periodic boundary-value theory no | Moderate analysis, self-contained |
| (ii-c) "Overdetermined at O(r³)" | A genuine theorem would need a genericity notion and a proof that the system is inconsistent for generic axes | n/a | **Not a theorem in the sources**; counting unknowns is a heuristic. Formalizing it means first proving something new |
| (iii) Equivalence of quasisymmetry definitions | Local: product rule, second-derivative symmetry, cross-product identities, Lie bracket. Global (Boozer form): flux coordinates, straightening of field lines on tori, circle action | Local: yes. Global: no (no Arnold-Liouville, no flows as diffeos, no Lie derivative) | Local equivalences **feasible**; Boozer-coordinate equivalence needs hypotheses or a large build |
| (iv) Axisymmetric (Grad-Shafranov) equilibria are MHS and quasisymmetric | Chain rule in cylindrical coordinates away from the axis; one explicit example (Solov'ev polynomial) | Yes | **Feasible**, pure calculus |
| Existence / non-existence of exact non-axisymmetric quasisymmetric MHS equilibria | Unknown: either an explicit construction, or a rigidity proof | n/a | **Open research problem** |

### Inferences
- "Mostly algebra / formal power series, plausibly feasible": (iii-local), (iv), (ii-a) with Boozer coordinates assumed, checking explicit solutions. "Needs absent analysis": Boozer/flux coordinates from first principles, any existence theorem, convergence of the near-axis series, anything KAM-flavoured.
- The Burby et al. Killing-field question is attractive for AI+Lean because the authors themselves say it may be local: a local statement about smooth germs in R³ can be stated with what Mathlib has. It is still an open problem, not a formalization task.

### Gaps
- I did not read Garren-Boozer 1991 itself (not on arXiv); the counting argument is reported via Landreman-Sengupta's description.
- I did not check whether anyone has turned the O(r³) obstruction into a rigorous genericity theorem; if one exists, (ii-c) changes from "new mathematics" to "formalize a paper".

---

## 3. Existing formalizations (physics, fluids, PDE, differential geometry, plasma)

### Takeaway
Until 2025 there was essentially nothing on plasma or MHD in any proof assistant. In 2026 that changed completely: MHS equilibria, Vlasov-Maxwell-Landau and relativistic Vlasov-Maxwell all have Lean 4 formalizations, each built privately on Mathlib with AI writing the code, and none of it lives in Mathlib or Physlib. Physlib itself has electromagnetism and a small fluid module but zero plasma/MHD content.

### Cited Findings
- **MHS / Grad's conjecture in Lean.** Gómez-Serrano, Liehr, Taylor, "Counterexamples to Grad's conjecture" (21 Sep 2026, 147 pages): smooth MHS equilibria on embedded solid tori with nested toroidal pressure levels whose Euclidean symmetry group is exactly cyclic C_N; "A Lean 4 certification of the above results is also provided." [FH-paper] — [arXiv:2609.24739](https://arxiv.org/abs/2609.24739)
- Their LLM disclosure: work "began in July 2026"; the authors "set up the problem and provided a detailed roadmap"; they used "GPT-5.6 Sol, Claude Fable 5, and Claude Opus 5 to fill in the technical details"; "The Lean code was produced by the same models with continuous guidance from the authors"; later "GPT-6 Astra and Claude Fable 5.1". [FH-paper] — [arXiv:2609.24739](https://arxiv.org/abs/2609.24739)
- The repository: 5,847 Lean files, 566,654 lines (my count), Lean 4.33.1. `LeanCode/` contains `NashMoser/` (10 files), `SmoothDependence/`, `FourierAnalysis/` (68 files, 14,940 lines), `AnalyticAlgebras/`, `BanachCalculus/`, `Topology/`, `Applications/MHS/` (180 files, 26,714 lines) and `Internal/` (5,530 files, 517,453 lines: weighted function spaces, weak derivatives, extension operators, tame estimates, the linear inverse). The statement file is 242 lines. [FH-repo] — [lukasliehr/Grad-Conjecture](https://github.com/lukasliehr/Grad-Conjecture)
- **Explicit counterexamples, no Lean.** Landreman (22 Sep 2026): explicit non-axisymmetric MHD equilibria with nested toroidal surfaces "given explicitly in Cartesian coordinates using elementary functions"; "These solutions were discovered using the artificial intelligence model GPT-6 Astra Pro … All equations were confirmed manually by the author." He names as open: vacuum or force-free analogues, and "the existence of non-axisymmetric equilibria satisfying exact quasisymmetry". [FH-paper] — [arXiv:2609.26742](https://arxiv.org/abs/2609.26742)
- de Lacy and Sato (12 Nov 2026): real-analytic toroidal MHD equilibrium with trivial Euclidean symmetry group; "discovered with the assistance of GPT-6 Astra Ultra"; Mathematica code, no Lean. [FH-paper] — [arXiv:2611.16688](https://arxiv.org/abs/2611.16688)
- Burby, Cho, Ruth (26 Oct 2026) write that recent counterexamples "show that Grad's conjecture is false. However, a more practical weakening of" it still stands and the counterexamples "fall far short of addressing the weakened conjecture". [FH-paper] — [arXiv:2610.36241](https://arxiv.org/abs/2610.36241)
- **Plasma kinetics in Lean.** Ilin, "Semi-Autonomous Formalization of the Vlasov-Maxwell-Landau Equilibrium" (Mar 2026): 34 files, 10K+ lines; "The theorem cannot even be stated with existing Mathlib – the Landau operator, torus differential structure, and VML system had to be defined from scratch"; torus represented as `Fin 3 → AddCircle 1` with hand-rolled `torusCurlX`, `torusDivX`; roughly half the code is supporting infrastructure. [FH-paper] — [arXiv:2603.15929](https://arxiv.org/abs/2603.15929)
- OpenAI's public `openai/math` repository lists family 362, "Global smoothness for relativistic Vlasov-Maxwell", with a Lean comparator challenge `VlasovMaxwell.lean` that imports only Mathlib and hand-defines `cross`, `divergence`, `curl` via `fderiv` on `EuclideanSpace ℝ (Fin 3)`; also family 376 (computation in forced Navier-Stokes flows) and 363 (Boltzmann non-uniqueness), each with Lean. Keyword counts in its catalogue: `magnetohydr` 0, `stellarator` 0, `plasma` 0, `Beltrami` 0, `Grad's` 0. [FH-repo] — [openai/math](https://github.com/openai/math)
- A further 2026 paper ships Lean 4 for a thermoviscous MHD Lagrangian of "sheared magnetic tori" (Ren, 27 Oct 2026). I read only the title page; the search summary says the Lean part excludes the full Maxwell equations and the toroidal embedding. [FH-paper for existence; snippet for scope] — [arXiv:2610.37305](https://arxiv.org/abs/2610.37305)
- **Physlib (formerly PhysLean/HepLean)**, commit of 25 Nov 2026, Lean v4.36.0. Line counts by directory (mine): Relativity 30,032; QFT 20,576; Particles 19,473; SpaceAndTime 17,830; Mathematics 17,504; QuantumMechanics 14,189; ClassicalMechanics 13,833; **Electromagnetism 9,575 (36 files: potentials, field strength, Lagrangian, Maxwell's equations in 3D, plane waves, infinite wire, circular coil)**; **FluidDynamics 1,361 (17 files: continuity, Cauchy momentum, Euler, Navier-Stokes as definitions plus `euler_iff_convectiveEuler`)**. A case-insensitive grep of all `.lean` and `.md` files for `plasma|magnetohydro|MHD|stellarator|tokamak|Grad-Shafranov|Beltrami|force.free` returns nothing. [FH-repo] — [leanprover-community/physlib](https://github.com/leanprover-community/physlib)
- PhysLean was merged with Lean-QuantumInfo and renamed Physlib under leanprover-community. [snippet; consistent with the repo layout I cloned, where `HEPLean/PhysLean` resolves to the same commit] — [github.com/leanprover-community/physlib](https://github.com/leanprover-community/physlib)
- **PDE in Lean, 2026.** Armstrong and Kempe formalized interior De Giorgi-Nash-Moser theory, "to our knowledge, the first machine-checked formalization of a major theorem in modern PDE theory", requiring "substantial new infrastructure for Sobolev spaces on bounded domains". [FH-paper, abstract] — [arXiv:2604.05984](https://arxiv.org/abs/2604.05984)
- Soto Franco and Marshall-Stevens formalized the Dirichlet problem for divergence-form elliptic operators (Poincaré inequality, Lax-Milgram, Rellich-Kondrachov, Fredholm alternative, interior regularity, Sobolev embedding) with "a self-contained theory of Sobolev spaces developed independently of existing formalisations". [FH-paper, abstract] — [arXiv:2609.32561](https://arxiv.org/abs/2609.32561)
- Cârstea and Oksanen formalized a special case of an anisotropic Calderón uniqueness theorem that OpenAI had released without Lean; "Our formalization uses results from Mathlib and Tau Ceti"; "We used OpenAI Codex, and to lesser extent Claude Code". [FH-paper] — [arXiv:2611.10537](https://arxiv.org/abs/2611.10537)
- Stokes' theorem for smooth singular cubes in Lean 4, with a bridge to Mathlib's `extDeriv` (May 2026). [FH-paper, abstract] — [arXiv:2605.01028](https://arxiv.org/abs/2605.01028)
- **Other provers.** Isabelle/HOL has a formalization of Green's theorem (Abdulaziz and Paulson; Archive of Formal Proofs entry "Green"); HOL Light has complex vector analysis and electromagnetic optics (law of reflection for plane waves, libraries of optics rules). [snippet] — [AFP Green](https://www.isa-afp.org/browser_info/current/AFP/Green/document.pdf), [arXiv:1403.3039](https://arxiv.org/abs/1403.3039)

### Inferences
- The answer to "is anything on magnetic confinement formalized?" is now yes, and by a wide margin the largest item is exactly on this topic. Any new project should start from `Showcase.lean`'s definitions rather than inventing its own `IsMHS`.
- No one appears to have formalized quasisymmetry, the near-axis expansion, Grad-Shafranov, or Boozer coordinates: zero hits in Physlib, in OpenAI's catalogue, and in my searches. That is the unoccupied niche. (Absence of evidence from keyword searches, so moderate confidence.)
- Duplication risk is real and fast-moving: three independent Grad counterexamples appeared within eight weeks.

### Gaps
- I did not verify Isabelle ODE/flow libraries (Immler's work) or Coq PDE work in this session; my search returned nothing specific on them. Treat any statement about them as unverified background.
- I did not check whether the Grad-Conjecture `Showcase.lean` has been independently reviewed for definition alignment (for example its notion of "embedded torus" or of the moduli topology).

---

## 4. What AI-assisted formalization can actually do (state of November 2026)

### Takeaway
Proof-writing is no longer the bottleneck: 2026 systems produce Lean developments of 10⁴ to 10⁷ lines. The documented failures are all in **statements and definitions** (wrong definition, weakened theorem, vacuous hypotheses, unusable API), and every credible project had a domain expert reading the statement file. I found no documented case of a non-expert on a consumer plan producing a Lean-verified research-level result end to end.

### Cited Findings

Scale and steering, by project:

- **Ilin, Vlasov-Maxwell-Landau (Mar 2026).** One mathematician, 10 days, about 50 hours of active supervision, 229 prompts, 213 commits, zero lines of code written by hand; Claude Code (Opus 4 era) on a $200/month Claude Max plan, plus Harmonic's Aristotle (220 submissions, 111 proved, 28 disproved) and Gemini DeepThink for the informal proof. API-equivalent cost about $6,300 with caching, about $42,700 without. [FH-paper] — [arXiv:2603.15929](https://arxiv.org/abs/2603.15929)
- Its failure modes, in the author's words: "hypothesis creep" (the main theorem "had accumulated 42 hypotheses", later cut to 12); "premature sorry", "excessive heartbeats", "dead code accumulation", "over-engineering", "no-op cycles". After release, Mathlib maintainer Sébastien Gouëzel "immediately identified a definition-alignment bug: the original theorem used `ContDiff ℝ ⊤`, which in Mathlib means analytic, not C^∞ as intended. Neither the AI agent nor the mathematician caught this during 10 days", and remarked "if there is one in the only line I've checked, how many others are there?" The author's conclusion: "compilation is not correctness. Expert review of definitions and theorem statements remains essential". [FH-paper] — [arXiv:2603.15929](https://arxiv.org/abs/2603.15929)
- The same paper's comparison table: Math Inc.'s Gauss, 181K lines of Lean of which about 22K (definitions, statements) were written by human experts and about 160K lines of `sorry` were eliminated by the AI in 3 weeks, team of 6+; AxiomProver, fully autonomous proofs of 2.7K to 4.3K lines; Numina, about 8K lines with human-AI pair programming. [FH-paper] — [arXiv:2603.15929](https://arxiv.org/abs/2603.15929)
- **Ilin and Nugent, "Sorries Are Not the Hard Part" (Jun 2026).** A semi-autonomous formalization of Grothendieck's vanishing theorem "compiles with no sorries, but an expert review found serious problems in definitions, theorem generality, file organization, and the API"; agents "adapted well to local, mechanically checkable feedback, but remained weak at choosing definitions and designing APIs". [FH-paper, abstract] — [arXiv:2606.13925](https://arxiv.org/abs/2606.13925)
- **Clason, Helin, Nursultanov, Oksanen (Nov 2026).** "AI agents have recently made writing such proofs cheap … the remaining human task is to confirm that each formal statement says what was intended … the risks lie in the statements rather than in the proofs", and automated tests on statements "can expose problems even in the established literature". [FH-paper, abstract] — [arXiv:2611.06962](https://arxiv.org/abs/2611.06962)
- **Survey (Jul 2026):** "checking a formal proof does not by itself establish that it captures the intended problem or makes a new mathematical contribution". [FH-paper, abstract] — [arXiv:2607.07779](https://arxiv.org/abs/2607.07779)
- **Google DeepMind formal proof search (May 2026):** "Our most capable agent autonomously resolved 9 of 353 open Erdős problems at the per-problem cost of a few hundred dollars, proved 44/492 OEIS conjectures". [FH-paper, abstract] — [arXiv:2605.22763](https://arxiv.org/abs/2605.22763)
- **Formal Conjectures benchmark:** 7,559 Lean statements including 1,642 open conjectures; "AI-generated proofs and disproofs serve as a valuable auditing mechanism" for mis-stated problems. [FH-paper, abstract] — [arXiv:2605.13171](https://arxiv.org/abs/2605.13171)
- **LeanMarathon (Jun 2026):** long-horizon autoformalization "fails not only at hard lemmas, but at scale: statements drift, dependencies tangle, context decays, and local repairs corrupt distant work"; its harness formalized seven theorems from two papers (258 lemmas) autonomously. [FH-paper, abstract] — [arXiv:2606.05400](https://arxiv.org/abs/2606.05400)
- **Rethlas/Archon (Apr 2026):** resolved "an open problem in commutative algebra" and verified it in Lean 4 "with essentially no human involvement". [FH-paper, abstract] — [arXiv:2604.03789](https://arxiv.org/abs/2604.03789)
- **Anthropic, zeta zeros (Aug 2026).** An unreleased research version of Claude, prompted by Jarred Sumner, "a non-mathematician" on staff; two Claude Code sessions, about 31 million output tokens; after a re-prompt, "a day and a half coordinating about 60 subagents", 2,400 shell commands; bound raised from 41.6% to 67.2%. Anthropic mathematicians Levent Alpöge and Ralph Furman examined and validated; outside experts Conrey and Goldston looked at it. Lean formalization by Claude with Eric Easley. [FH, fetched page summary] — [anthropic.com/research/riemann-zeta](https://www.anthropic.com/research/riemann-zeta)
- The formalization's own metadata: "All Lean code was written by Claude; the trusted statement file was read against the paper by a paper author"; 326 Lean files under `zeta23/`; statements "made against Mathlib's `riemannZeta` and `DirichletCharacter.LFunction`"; 17 theorem statements in a Mathlib-only `Challenge.lean`; parts adapted from the PrimeNumberTheoremAnd project. [FH-repo] — [anthropics/zeta-23-lean (now "formal-math")](https://github.com/anthropics/zeta-23-lean)
- **OpenAI `openai/math`.** README: manuscripts "produced by an unnamed internal OpenAI model"; "730 manuscripts organized into 375 families"; "On average, each result used three hours of ChatGPT Pro thinking compute with that model"; "the model was posed approximately 4,000 problems"; "~42% top-line results formalized"; "Some of the unformalized results could have issues". Exceptions to the fixed procedure "include work on a zero-free region for the Riemann zeta function and proof of the Hodge Conjecture for CM abelian varieties". The repository holds 189,309 Lean files (my count). [FH-repo] — [openai/math](https://github.com/openai/math)
- Its quasi-Riemann entry: "every Dirichlet L-function, including ζ(s), is zero-free in Re s > 7/8", with comparator challenge files. [FH-repo] — [openai/math lean/docs/003.md](https://github.com/openai/math/blob/main/lean/docs/003.md)
- **RH/GRH.** Sernau's repository README: "A Lean 4 proof of the Riemann Hypothesis and of the Generalized Riemann Hypothesis for Dirichlet L-functions, built on Mathlib … based on Levent Alpöge's, with Claude … This repository is a shorter version, produced in collaboration with my colleagues at Anthropic, and Claude"; "readable edition, about 300,000 lines"; compact edition "about 76,000 lines"; build "one to two hours on 12 cores"; "The proof has passed machine checks only … No mathematician has read it through, and no paper comes with it." [FH-repo] — [lsernau/grh-lean](https://github.com/lsernau/grh-lean)
- Alpöge's original release: 20 Nov 2026, 175,187 Lean files, about 17 million lines, produced with an unreleased frontier Claude model, development 6 Oct to 16 Nov; Lean FRO and Mathlib Initiative issued a joint statement on 21 Nov certifying statement and proof; compute not disclosed, "mathematicians say the effort likely cost millions of dollars" (an estimate). [snippet] — [Scientific American](https://www.scientificamerican.com/article/the-ai-riemann-hypothesis-proof-is-almost-certainly-true-independent-mathematicians-attest/), [Wikipedia](https://en.wikipedia.org/wiki/Alp%C3%B6ge's_proof_of_the_Riemann_hypothesis)
- **Amateur case.** Liam Price, 23, no advanced training, pasted Erdős problem #1196 into GPT-5.4 Pro and got a proof in about 80 minutes (April 2026); Lichtman said the raw output "was quite poor – it took an expert to understand what the model was trying to say"; Tao and Lichtman reviewed and shortened it; Math Inc. formalized a version in Lean; a ChatGPT Pro subscription was gifted to him afterwards. [snippet] — [Scientific American](https://www.scientificamerican.com/article/amateur-armed-with-chatgpt-vibe-maths-a-60-year-old-problem/), [heise](https://www.heise.de/en/news/Creative-solution-AI-solves-60-year-old-Erd-s-problem-11276442.html)
- **Jacobian conjecture (Jul 2026).** Counterexample in three dimensions found by Alpöge with Claude Fable 5; Tao published a human-readable derivation; formalized in Lean. Prompting details not public. [snippet] — [ScienceDaily](https://www.sciencedaily.com/releases/2026/08/260804034634.htm)

### Inferences
- The realistic division of labour in late 2026: the model writes all Lean; a human who understands the mathematics owns a short statement file and reads every definition in it. A "technically capable non-mathematician" can run the loop but cannot supply the second half, so statement review must be bought from somewhere (Zulip, a plasma theorist, adversarial model review, non-vacuity theorems).
- My own `⊤` versus `∞` slip (section 1) shows this is not a hypothetical risk for this exact topic and this exact model family.
- Cheap mechanical guards that the sources converge on: (1) Mathlib-only `Challenge.lean` plus comparator; (2) `#print axioms`; (3) a non-vacuity theorem exhibiting an object that satisfies all hypotheses; (4) a ban on new hypotheses in the main theorem without sign-off; (5) known-false sanity statements that must fail.

### Gaps
- No primary access to Lean Zulip, Tao's or Buzzard's blogs, Quanta or Scientific American in this session; I have no direct quotes from Tao, Buzzard, Massot or Commelin beyond what is embedded in the papers above (Gouëzel's quote is from the VML paper).
- Compute and steering for the Alpöge RH proof are undisclosed in everything I could reach.

---

## 5. Would it be used, and how would an outsider contribute?

### Takeaway
Mathlib is effectively closed to this project as described (AI-written, non-expert author). The open doors are PhyslibAlpha, a standalone repository with a comparator-style statement file, and a roadmap contribution to Tau Ceti. Demand from the plasma side is real but new and thin: the people doing rigorous stellarator mathematics have just started attaching Lean to papers.

### Cited Findings
- Mathlib contribution guide (source file, current): "As of mid-2026, code written by an AI without the supervision of a Lean subject expert fails to meet that bar by a large margin. Members of the review team will summarily close without comment any low quality PR produced using LLMs, especially if the author has made little effort to directly engage in the community"; repeated offences lead to suspension or ban; "Getting code to mathlib's standards requires understanding and writing Lean code by hand"; AI use "must" be disclosed in the PR description and substantial LLM code labelled `LLM-generated`; "Using an LLM when writing comments on GitHub or Zulip is not allowed"; "being able to justify each decision to reviewers without the use of an AI". [FH, raw file from GitHub] — [leanprover-community.github.io/contribute](https://github.com/leanprover-community/leanprover-community.github.io/blob/lean4/templates/contribute/index.md)
- Mathlib has roughly 300 PRs in backlog and median waits around two weeks. [snippet; the search summary did not make clear which of its results this figure came from, the "Growing Mathlib" paper being the likeliest] — [arXiv:2508.21593](https://arxiv.org/abs/2508.21593)
- Physlib `AI-POLICY.md`: "Physlib generally welcomes AI-assisted contributions"; "The human author(s) are fully responsible for every line … A clean Lean build proves what was written, but only the human can certify what was written is what they meant"; bibliographic references must be verified by the human; "All communication with human reviewers must be conducted by humans". `AGENTS.md` bans `axiom`, `sorry`, `True`-valued fields, and hypotheses that assume the conclusion, and wants proofs under 50 lines. [FH-repo] — [physlib AI-POLICY.md](https://github.com/leanprover-community/physlib/blob/master/AI-POLICY.md)
- Physlib README: "PhyslibAlpha exists for the rapid development of physics digitalizations, enabled by a lighter review process built to handle large-scale, human- or AI-generated contributions" (849 files at the commit I cloned). [FH-repo] — [leanprover-community/physlib](https://github.com/leanprover-community/physlib)
- Tau Ceti: "human-directed and human-reviewed roadmaps, implemented and maintained by AI contributors"; "AIs own the code in this repository"; "Humans can raise issues against the code"; incubated by the Lean FRO and the Mathlib Initiative. Its README also points to Lean Pool, "an arXiv/AfP-like archive of independent formalization projects, whether human- or AI-written". [FH-repo] — [TauCetiProject/TauCeti](https://github.com/TauCetiProject/TauCeti)
- Anthropic's formal-math README describes the Palomar registry convention: a Mathlib-only `Challenge.lean`, an untrusted `Solution.lean`, `comparator.json`, `formalization.yaml`, checked by `leanprover/comparator` with an independent kernel. [FH-repo] — [anthropics/zeta-23-lean](https://github.com/anthropics/zeta-23-lean)
- Evidence of demand in the plasma-mathematics community: the Grad counterexample paper devotes a section to explaining its Lean statement to non-Lean readers; Landreman cites that paper and says exact solutions are useful "for testing numerical codes"; Burby et al. cite the counterexamples within five weeks. [FH-paper] — [arXiv:2609.24739](https://arxiv.org/abs/2609.24739), [arXiv:2609.26742](https://arxiv.org/abs/2609.26742), [arXiv:2610.36241](https://arxiv.org/abs/2610.36241)

### Inferences
- Realistic contribution path, in order of friction: own GitHub repository in Palomar/comparator layout (none); PhyslibAlpha PR (low; a human must vouch for every statement); Tau Ceti roadmap for "MHS equilibria and quasisymmetry" (medium; needs a human-reviewed roadmap, then AIs implement); Physlib core (high); Mathlib (do not attempt).
- "Would it be used?" depends on the artifact. A definitions-and-equivalences library for quasisymmetry would be used only if the next rigorous paper in this area builds on it; a machine check of a specific published claim (Landreman's solutions) is self-contained and immediately citable.

### Gaps
- I could not read Physlib's Zulip or TODO list to see whether plasma/MHD has been requested there.
- I found no statement from any stellarator group (Simons collaboration, IPP, PPPL) about wanting formal verification; demand is inferred from the three papers above.

---

## 6. Candidate formalization targets, ranked by feasibility

### Takeaway
Four targets are feasible for Claude plus a non-mathematician on a subscription, provided statement review is obtained externally; one is feasible only as bookkeeping with assumed coordinates; the remaining open problem is not a formalization task at all. All effort figures are my **estimates**, calibrated on: my two one-shot probes; the VML project (10K lines, 10 days, 50 supervised hours with an older model); and the Grad counterexample (567K lines, three expert authors, about ten weeks).

### Cited Findings
- Calibration data are cited in sections 1, 3 and 4. — [arXiv:2603.15929](https://arxiv.org/abs/2603.15929), [arXiv:2609.24739](https://arxiv.org/abs/2609.24739), [Probe2.lean](file:///tmp/claude-0/-home-user/f7c08234-634f-5a50-bfc5-78c71b1b878c/scratchpad/lean/Probe2.lean)
- Landreman's first family is `B(x) = A·B₀(A⁻¹x)` with `A = diag(a,b,1)` and `B₀ = (I/R)e_φ + (2Z/R)e_R + (1−R²)e_Z`, `I = sqrt(1−4Ψ₀)`, `Ψ₀ = (R²−1)²/4 + Z²`; the key identity is `(B₀·∇)B₀ = −Dx` with `D = diag(1,1,4)`; field lines have an explicit parametrization `r(u,v,ζ)`; "All equations were confirmed manually by the author". [FH-paper] — [arXiv:2609.26742](https://arxiv.org/abs/2609.26742)

Ranked targets (class labels a/b/c are defined in section 7):

**T1. Lean certificate for Landreman's explicit non-axisymmetric equilibria. Class (c), light; the most feasible and most immediately useful.**
- Statement: for the explicit `B`, `p`, `ψ` of arXiv:2609.26742 family 1, on the explicit open set where the formulas are defined: `div B = 0`, `(curl B) × B = grad p`, `B·grad ψ = 0`, `grad p ≠ 0` off the axis, the explicit curve `r(u,v,ζ)` is an integral curve of `B` (closed field lines, ι = 2), and `B` is not invariant under rotations about the z-axis for `ε ≠ 0`. Stretch: level sets of ψ are embedded tori (via the explicit parametrization), matching `IsMHS`/`NestedPressureTori` from the Grad `Showcase.lean`.
- Prerequisites: all present in Mathlib (`fderiv`/`HasDerivAt`, `Real.sqrt`, trig). Probe 2 already did the divergence of the seed field in one shot.
- Effort [estimate]: local identities, 3 to 10 days; integral curves and non-axisymmetry, another 1 to 2 weeks; nested-tori topology, 2 to 6 further weeks and the main risk.
- Verifiably correct and useful if: statement file imports only Mathlib and is diffable against the paper's equations; Landreman (or a co-reader) confirms the transcription; non-vacuity is automatic (the object is explicit). It converts "confirmed manually by the author" into a machine check of an AI-discovered result, in a field that now expects that.
- Unknown: whether someone has already done it (two months have passed).

**T2. Axisymmetric Grad-Shafranov equilibria are MHS and quasisymmetric. Class (c), light.**
- Statement: for `ψ(R,Z)`, `F(ψ)`, `p(ψ)` satisfying the Grad-Shafranov equation on an open set with `R > 0`, the field `B = ∇ψ×∇φ + F∇φ` satisfies `div B = 0`, `(curl B)×B = ∇p`, and `u = ∂_φ` satisfies the Euclidean quasisymmetry conditions; instantiate with the Solov'ev polynomial `Ψ₀` (Landreman's eq. 2.6 to 2.7) as a non-vacuity witness.
- Prerequisites: present. Pain point: cylindrical-to-Cartesian chain rule.
- Effort [estimate]: 1 to 3 weeks. Useful as the base case and regression test for any quasisymmetry definition.

**T3. Local theory of quasisymmetry in Euclidean R³: definitions and equivalences. Class (c), moderate.**
- Statement: define quasisymmetry of `(B,u)` on an open set where `B ≠ 0` in vector-calculus form; prove the equivalence of the Burby-Kallinikos-MacKay condition sets (their Thm 4.4 (i) to (iii)) in that form, the consequences `[u,B] = 0`, `u·∇(u·B) = 0` (Thm 4.5), `B×u = ∇ψ` implies `B·∇ψ = 0 = u·∇ψ`, and that the triple-product form follows; prove Killing fields of R³ that preserve `B` are quasisymmetries.
- Prerequisites: `VectorField.lieBracket` present; Lie derivative of 1-forms and of the volume form must be hand-rolled in coordinates (absent in Mathlib).
- Effort [estimate]: 3 to 8 weeks, 3K to 10K lines. Risk is almost entirely definitional: the guiding-centre definition (Def. 4.1) versus derived conditions, the Rodríguez et al. "weak" versus Burby et al. "strong" notions, regularity classes.
- Verifiable if T2 and T1 are instances/non-instances and a plasma theorist signs the statement file. The Boozer-coordinate formulation is explicitly out of scope (needs flux coordinates).

**T4. A reusable MHS/quasisymmetry statement layer on top of the Grad `Showcase.lean`, plus a trusted statement of the open problem. Class (a)/(c) hybrid.**
- Content: port the 242-line showcase definitions to current Mathlib (it is pinned to Lean 4.33.1; the local project is 4.37.0-rc1), add a quasisymmetry predicate, state "there exists a smooth non-axisymmetric quasisymmetric MHS equilibrium with nested toroidal pressure surfaces" and its negation as `Challenge`-style files, prove non-vacuity lemmas (T2 satisfies everything except non-axisymmetry; T1 satisfies MHS).
- Effort [estimate]: 2 to 4 weeks. This is the prerequisite for any full-ambition attempt and has standalone value (it could be offered to Formal Conjectures or PhyslibAlpha).

**T5. Garren-Boozer near-axis expansion through second order as formal-series algebra. Class (c), heavy; feasible only with assumptions.**
- Statement: *assuming* a Boozer-coordinate structure (as a Lean `structure` of hypotheses), derive the O(r) and O(r²) constraint equations of Landreman-Sengupta and prove the O(r) system has a unique periodic solution for the σ equation given axis data.
- Prerequisites: `PowerSeries` with derivations present; a coefficient algebra of smooth periodic functions with `d/dφ` must be built; periodic ODE solvability must be built.
- Effort [estimate]: 2 to 4 months; high risk that the result is a faithful transcription of hypotheses with little content ("hypothesis creep" by design). The O(r³) "obstruction" cannot be included as a theorem because the literature gives a counting argument, not a proof.
- Useful mainly as a cross-check of the published equations (transcription errors in long expansions are plausible) and of codes such as pyQSC.

**T6. NOT feasible as a formalization target: existence or non-existence of exact non-axisymmetric quasisymmetric MHS equilibria (and the "weakened Grad conjecture"). Class (b).**
- Why not: there is no proof to formalize. Constantin-Drivas-Ginsberg called even the original conjecture "out of reach of existing techniques"; its refutation needed a Nash-Moser scheme, three professional analysts and 567K lines of Lean. Landreman lists exact quasisymmetry as still open after finding explicit non-symmetric equilibria with a frontier model, which suggests the obvious explicit-ansatz search has already been tried by the best-placed person. The Lean side is the easy part; the missing ingredient is the mathematics.
- What would change this verdict: a candidate explicit solution (then it becomes T1-shaped), or a new paper proving rigidity (then it becomes class (c)).

**Also not recommended:** convergence of the near-axis series, KAM persistence of flux surfaces, or building general differential-forms-on-manifolds infrastructure "for" this project. The first two are open or deep; the third duplicates Tau Ceti.

### Inferences
- A sensible sequence is T2 → T1 → T4 → T3, with T5 optional. Each earlier item is a test instance for the next one's definitions.
- T1 and T4 are the two with a plausible external audience within weeks.

### Gaps
- No estimate here has been tested beyond the two small probes. The nested-tori part of T1 could be much harder than guessed.
- I did not open Landreman's supplemental material; it may already contain symbolic verification that makes the Lean check less novel.

---

## 7. Scope update: a full-ambition attempt on an open problem with AI + Lean

### Takeaway
The 2026 breakthroughs split into three classes with very different evidence for individuals. (a) Re-formalizing or compressing an existing machine-checked proof: plausible for one person on a consumer plan. (b) Discovering a new proof of an open statement: done by labs with unreleased models, by experts steering released frontier models, and once by an amateur on a short Erdős problem whose output experts had to decode; never, that I can document, by a non-expert through to Lean verification. (c) Formalizing a known paper proof where Mathlib lacks infrastructure: done repeatedly, always with an expert reading the statements. Exact quasisymmetry existence is class (b), and unlike the Riemann-zeta results its statement is not a one-liner over existing Mathlib definitions.

### Cited Findings

Who did what, with what (rows cite the sources already given in section 4 unless a link is shown):

| Result (2026) | Who | Model access | Human steering | Time | Lean size | Label |
|---|---|---|---|---|---|---|
| Zeta zeros, 41.6% → 67.2% (Aug) | Anthropic; prompted by a non-mathematician staff member; validated by Anthropic number theorists Alpöge, Furman | Unreleased research Claude; 31M output tokens; about 60 subagents | Prompts were mostly encouragement; expert validation and a paper by experts afterwards | Two sessions; second about 1.5 days | 326 Lean files; all Lean by Claude; statement file read by a paper author | [FH page + FH-repo] |
| Quasi-RH, zero-free for Re s > 7/8 (30 Sep / 6 Oct) | OpenAI | "Unnamed internal OpenAI model"; explicitly an exception to the 3-hours-of-Pro-compute procedure | Not disclosed | Not disclosed | Part of a 189,309-file Lean library | [FH-repo] — [openai/math](https://github.com/openai/math) |
| 375 result families (Oct to Nov) | OpenAI | Same internal model; about 3 hours of ChatGPT Pro thinking compute per result; about 4,000 problems posed | Fixed procedure; curation by significance | Weeks | about 42% of top-line results formalized | [FH-repo] |
| RH and GRH for Dirichlet L-functions (20 Nov) | Alpöge (Anthropic mathematician) with Claude | Unreleased frontier Claude | "Directed the project toward a promising approach" | 6 Oct to 16 Nov | about 17M lines, 175,187 files | [snippet] |
| RH/GRH compressed (by 25 Nov) | Sernau and Anthropic colleagues with Claude | Not stated | Not stated | Days | about 300K readable / 76K compact | [FH-repo] — [lsernau/grh-lean](https://github.com/lsernau/grh-lean) |
| Same compression, independent | Hacker News user "jychang" | "GPT-6.2 Sol" on a $100 ChatGPT plan, one agent, no subagents, human nudges | Self-described | Days | "about 70% done", 35K lines | **[relayed]; self-reported and unverified; I did not see it** |
| Counterexamples to Grad's conjecture (21 Sep) | Gómez-Serrano, Liehr, Taylor (professional analysts) | Released commercial models: GPT-5.6 Sol, Claude Fable 5, Claude Opus 5, later GPT-6 Astra, Claude Fable 5.1; plan tier not stated | Authors wrote the roadmap; "continuous guidance" on Lean | July to 21 Sep | 566,654 lines, 5,847 files | [FH-paper + FH-repo] |
| Explicit non-symmetric equilibria (22 Sep) | Landreman (leading stellarator theorist), alone | GPT-6 Astra Pro (a consumer-tier product name) | Expert; "All equations were confirmed manually" | Not stated | None | [FH-paper] |
| Trivial-symmetry analytic equilibrium (12 Nov) | de Lacy, Sato | GPT-6 Astra Ultra | Experts; manual verification | Not stated | None (Mathematica) | [FH-paper] |
| Vlasov-Maxwell-Landau steady states (Mar) | Ilin (mathematician), alone | $200/month Claude Max + free Aristotle + Gemini | 229 prompts, about 50 hours | 10 days | 10K+ lines | [FH-paper] |
| Erdős #1196 (Apr) | Liam Price (amateur) | GPT-5.4 Pro, single prompt | None during; experts decoded and rewrote; Math Inc. formalized | about 80 minutes of model time | Not stated | [snippet] |
| 9 of 353 open Erdős problems (May) | Google DeepMind | Lab agent; "a few hundred dollars" per problem | Autonomous | n/a | n/a | [FH-paper abstract] |
| Jacobian conjecture counterexample (Jul) | Alpöge | Claude Fable 5 (released) | Prompting not public | n/a | Formalized afterwards | [snippet] |

- Is there any documented case of a non-expert on a consumer subscription producing a Lean-verified research-level result? **I found none in which the non-expert carried it through to Lean.** Closest: Price (amateur discovery, but experts interpreted it and a company formalized it); Sumner (non-mathematician prompter, but unreleased model, lab compute, expert validation); Ilin (consumer plan and no hand-written code, but a professional mathematician, and an expert still found a statement bug afterwards); the relayed Hacker News comment (consumer plan, class (a), unverified, expertise unknown). [synthesis of the rows above]
- The "$100 OpenAI plan, 70% of the way to GRH/RH" story: per the coordinator it is a Hacker News self-report about **re-deriving a compressed version of an already machine-checked proof**, not progress on an open problem. My own searches found no trace of it; I can confirm only that the compression it parallels exists (`lsernau/grh-lean`, 76K-line compact edition). [relayed + FH-repo]

What the successful cases had in common:
- **Statement expressible over existing Mathlib definitions, in a few lines.** Zeta23's statements are "made against Mathlib's `riemannZeta` and `DirichletCharacter.LFunction`, with multiplicity via `analyticOrderAt`" and need only 15 auxiliary definitions; the RH proof targets Mathlib's own statement of RH. [FH-repo; snippet for RH] — [anthropics/zeta-23-lean](https://github.com/anthropics/zeta-23-lean)
- **A rich neighbouring library to adapt.** Zeta23 adapts files from PrimeNumberTheoremAnd. [FH-repo]
- **Where Mathlib lacked the area (PDE), the statement was still made self-contained in a short Mathlib-only file with elementary definitions,** and all heavy analysis was built privately: OpenAI's `VlasovMaxwell.lean`; Grad `Showcase.lean`. [FH-repo]
- **A trusted-statement / untrusted-proof split with comparator and an axioms audit** in every lab release. [FH-repo]
- **An expert owned the mathematics** in every case where the provenance is documented; the two partial exceptions (Price, Sumner) were followed by expert validation before anyone relied on the result.
- **Explicit-object problems fell first in this very field:** two of the three Grad counterexamples are closed-form fields found by a model and checked by hand; the third is an existence proof by Nash-Moser.

Does the stellarator problem share those properties?
- Shares: the statement can be made elementary and Mathlib-only (template exists); if a solution is an explicit formula, certification is cheap (Probes 1 and 2).
- Does not share: (1) there are several inequivalent definitions of quasisymmetry and of "non-trivial", so the one-line-trusted-statement property fails until T3/T4 are done and reviewed; (2) no neighbouring Lean library exists for quasisymmetry (zero hits everywhere I looked); (3) the expected answer is unknown: Garren-Boozer's expansion suggests generic non-existence, and a non-existence (rigidity) proof would be hard global analysis of the kind Constantin-Drivas-Ginsberg called out of reach, while Burby et al. hint the Killing-field question might be local; (4) the easiest half of the neighbouring problem was just harvested by experts with stronger model tiers than a standard subscription.

Which class each stellarator target is in, and the evidence for individuals succeeding in that class:

| Target | Class | Evidence that an individual on a consumer plan can do this class |
|---|---|---|
| Port/compress/re-derive parts of the Grad-conjecture Lean development; port its statement layer (part of T4) | (a) | RH compression by a small Anthropic group [FH-repo]; the relayed, unverified HN self-report. Weak but positive. |
| T1 Landreman certificate; T2 Grad-Shafranov; T3 local quasisymmetry equivalences; T5 near-axis algebra | (c) | Ilin's VML project (one person, $200 plan, 10 days) is the direct precedent, with the caveat that he is a mathematician and still shipped a statement bug. De Giorgi-Nash-Moser, elliptic PDE, Calderón: all expert teams. Positive for experts; **untested for non-experts**. |
| T6 exact quasisymmetry existence/non-existence; weakened Grad conjecture; "every quasisymmetry is Killing" | (b) | Experts with consumer-tier frontier models: Landreman, de Lacy-Sato (explicit constructions). Amateur: Price, on a short combinatorial-number-theory problem. Non-expert + deep analysis + Lean: **no documented case**. |

Realistic shape of a full-ambition attempt (my proposal; [estimate] throughout):

1. **Build first (weeks 1 to 4): the trusted statement layer (T4) and its test instances (T2, then T1).** Without these there is nothing to aim a model at and no way to know whether a "proof" proves the intended thing. Use `∞` not `⊤`; write non-vacuity lemmas; keep the statement file under about 300 lines; get one plasma theorist or the Lean Zulip to read it. Decide explicitly between the Burby et al. and Rodríguez et al. notions and state both.
2. **Calibrate the model before trusting a negative (week 2 to 3).** Give the subscription model the problem Landreman gave his (find explicit non-axisymmetric MHS equilibria with nested surfaces) without the September papers in context. If it cannot rediscover a known family, its failure to find a quasisymmetric one carries no information.
3. **Search the explicit side (weeks 3 to 8).** Ask for closed-form or finitely-parametrized candidates (deformations of axisymmetric seeds, as in Landreman's `A·B₀(A⁻¹x)` construction; the families of all three counterexample papers) and test `u·∇|B|` or the triple product symbolically, then numerically with the near-axis and equilibrium codes already cloned by other researchers in this project (pyQSC, DESC). A candidate whose quasisymmetry residual vanishes identically in a CAS is the only event that would justify full Lean effort; certifying it would then take days to weeks.
4. **In parallel, attack a local rigidity statement in Lean:** "a quasisymmetry of an MHS field on an open subset of R³ is a Killing field", under progressively weaker extra hypotheses, starting from Burby et al.'s Theorem 7.2. This is the one sub-question the literature flags as possibly local and therefore statable today. Use Aristotle-style disproof attempts on every intermediate lemma, as Ilin did (28 of 220 submissions were disproved, which he found more valuable than the proofs).
5. **What must exist before a global non-existence proof could even be formalized:** flux-surface/Boozer coordinates, flows as diffeomorphisms with smooth dependence, averaging over tori, probably elliptic theory on the torus. Some of it is in Tau Ceti; none is in Mathlib. This is the 100K-plus-line regime of the Grad certification and should not be started speculatively.

Early signals that it is or is not going anywhere:
- Positive: the statement file survives an outside expert's read with no semantic changes; T2 and T1 go through in days; the model rediscovers a known counterexample family unaided; a candidate's residual is identically zero symbolically; a restricted Killing-field lemma is proved with hypotheses that T1/T2 instances satisfy.
- Negative (stop or rescope): the main theorem's hypothesis list grows from week to week; the model "proves" non-existence quickly (almost certainly a vacuous or mis-stated theorem; test it against the axisymmetric instance, which must *not* be excluded by mistake); no rediscovery of known families; every candidate has a non-zero residual at third order in distance from the axis, mirroring Garren-Boozer; four weeks pass with only infrastructure and no mathematical statement about quasisymmetry proved.
- A hard budget rule borrowed from the documented projects: if phase 1 is not finished and externally reviewed in a month, the bottleneck is statement review, not compute, and a stronger plan will not fix it.

### Inferences
- The step-change of late 2026 came from model capability on (b), mostly inside labs with unreleased models. A subscription user gets the released models, which have demonstrably sufficed for expert-steered (b) in this exact area (Landreman) and for (c) (Ilin). What the subscription user lacks relative to those cases is the expert, not the model.
- The most defensible "ambitious" framing is therefore: build the statement layer and certificates (class (a)/(c), likely to succeed, useful regardless), and run the explicit-solution search and the local Killing-field question as cheap, clearly-bounded class (b) bets on top.
- Because this field moved from "conjecture" to "three counterexamples" in eight weeks, any class (b) target here should be assumed to be under attack by better-resourced groups already.

### Gaps
- The Hacker News comment, Quanta's and Scientific American's articles, Wikipedia's page and the Lean FRO/Mathlib Initiative statement were not accessible; their content above is [snippet] or [relayed].
- Which subscription tier the Grad-counterexample authors used, and how many model-hours, is not stated in their paper.
- I found one search hit (we0.ai, "GPT-6 Astra and Grad's Conjecture") that looked like secondary or generated commentary and did not use it.
- OpenAI's and Anthropic's repositories contain extraordinary claims beyond this topic (for example a refutation of the Exponential Time Hypothesis). I report only what the READMEs say and did not assess correctness; OpenAI's own README warns that unformalized results "could have issues".
- No search surfaced any 2026 attempt, by AI or otherwise, specifically on exact quasisymmetry existence; I cannot tell whether that is because none exists or because it is unpublished.
