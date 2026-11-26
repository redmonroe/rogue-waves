# Stellarator quasisymmetry for "Lean + Claude": tractability, objections, alternatives (state of things 26 Nov 2026)

Evidence labels used throughout:
- **[FH]** first-hand: I ran it on this 4-core / 15 GB machine, or read the cloned repo / the arXiv PDF or abstract page directly on 26 Nov 2026.
- **[WF]** a WebFetch summary of the page (a small model read the page for me; treat numbers as probably right, wording as paraphrase).
- **[SS]** search-result snippet only; I did not open the page. Lower confidence.

Access limits that shaped the notes: en.wikipedia.org and the GitHub REST/GraphQL API were not reachable, so issue trackers (DESC, SIMSOPT, pyQSC) and the Hugging Face leaderboard were NOT read. arXiv, `git clone` and PyPI worked.

Scope note: the coordinator changed the brief mid-task. The user wants to attempt a real open problem outright, modelled on 2026 AI + Lean results, not make a minor contribution. Sections 4 to 6 therefore assess each candidate both as a modest target and as a full-ambition target.

---

## 1. Open tools and data: what exists, licences, maturity, install size, runtimes

### Takeaway
Everything needed for near-axis and fixed-boundary stellarator work installs from PyPI in under a minute and runs in milliseconds to seconds on 4 cores; the tools are MIT/BSD, and the big ones (SIMSOPT, DESC, VMEC++, ConStellaration) had commits in the last five days. None of them uses or mentions Lean.

### Cited Findings
- **pyQSC** (PyPI name `qsc`, v0.1.3). [FH] `pip install qsc` in a fresh Python 3.13 venv took 0.6 s and installed a 73.6 kB pure-Python wheel. The wheel declares **no dependencies** (`install_requires` is commented out in `setup.py`), so `import qsc` fails with `ModuleNotFoundError: numpy` until numpy, scipy and matplotlib are installed by hand (31 s; venv then 261 MB). Licence file is BSD-style ("Redistribution and use in source and binary forms…", copyright 2020 Matt Landreman). Last commit on the default branch: 2024-10-19, so the repo looks dormant. — [pyQSC repo](https://github.com/landreman/pyQSC)
- [FH] pyQSC run results on this machine: all seven built-in paper configurations (`Qsc.from_paper("r1 section 5.1")` … `"r2 section 5.5"`) solved. First-order solves took 5–6 ms each; second/third-order solves 0.06–0.30 s. Example outputs: "r1 section 5.1" nfp=3, iota=0.418307; "r2 section 5.2" nfp=2, iota=−0.423724, r_singularity=0.4096. A resolution scan (nphi 31→481) held iota to 12 digits. A random survey of **1,000 first-order configurations took 3.97 s; 200 second-order configurations took 2.43 s** (single process). The Newton solver printed "did not get close to desired tolerance" for under 1% of random inputs (8 warnings across 1,200 solves) without raising an exception, so a survey has to filter on residuals itself. Total script wall time 11 s. — [pyQSC repo](https://github.com/landreman/pyQSC)
- **pyQIC** (PyPI name `qicna`): near-axis code that "does not assume quasisymmetry" (quasi-isodynamic). MIT licence, pure Python (numpy, scipy, matplotlib), last commit 2026-01-30. Not run here. [FH from clone] — [pyQIC repo](https://github.com/rogeriojorge/pyQIC)
- **SIMSOPT**: MIT, Python with C++/Fortran cores, `pip install simsopt`, depends on jax, scipy, sympy; MHD equilibrium needs separate VMEC or SPEC modules. Shallow clone 188 MB, last commit 2026-11-25. Ships the Landreman–Paul input files (`tests/test_files/input.LandremanPaul2021_QA`, `_QA_lowres`, `_QH_reactorScale_lowres`) and QA/QH example scripts, so a reproduction of the 2021 precise-quasisymmetry result is set up in the repo. Not installed or run here. [FH from clone] — [simsopt repo](https://github.com/hiddenSymmetries/simsopt)
- SIMSOPT now carries an `AGENTS.md` for AI coding agents: "A human has read every line… If an agent generated code that human cannot explain, it does not go in", and every change "must rest on the physics or numerical method… Cite the reference". PR template has an "I used AI to generate this PR" tick box. [FH] — [simsopt repo](https://github.com/hiddenSymmetries/simsopt)
- **DESC**: MIT, JAX-based equilibrium and optimisation code, `pip install desc-opt`, pinned to jax ≤ 0.12.1 and about 24 direct dependencies; tutorials include precise quasisymmetry, near-axis constraints from pyQSC/pyQIC, and coil optimisation. Shallow clone is 1.1 GB (includes publication data). Last commit 2026-11-24. `CONTRIBUTING.rst` points newcomers to "good first issue" and "help wanted" labels. Not installed or run here. [FH from clone] — [DESC repo](https://github.com/PlasmaControl/DESC)
- **VMEC++** (Proxima Fusion): MIT, C++ reimplementation of VMEC with Python API, "can run on a laptop". [FH] `pip install vmecpp` gave v0.9.4 in 8 s (venv 702 MB, pulls jax and netCDF4). Timings on this machine: Solov'ev tokamak <0.1 s; CTH-like fixed boundary 0.8 s; **W7-X fixed boundary (ns=99, mpol=ntor=12) 9.5 s wall, 42 s CPU**. — [vmecpp repo](https://github.com/proximafusion/vmecpp)
- **ConStellaration** (Proxima Fusion + Hugging Face): MIT code; dataset of QI-like plasma boundaries with VMEC++ equilibria; three benchmark problems (geometric; "simple-to-build" QI; multi-objective MHD-stable QI). The public leaderboard is frozen to `constellaration==0.2.6` and "evaluates every submission with exactly this version"; the README records that the leaderboard was mistakenly run on 0.3.0 from 2026-10-02 to 2026-10-08 and re-evaluated, which shows it is live. Install needs `libnetcdf-dev` and builds booz-xform; dev environment supports Python 3.10 only. Last commit 2026-11-21. Not installed here. [FH from clone] — [constellaration repo](https://github.com/proximafusion/constellaration); [paper, arXiv:2506.19583](https://arxiv.org/abs/2506.19583)
- ConStellaration compute, from the paper [FH]: SciPy baselines "ran for ∼40 hours on a machine with 4 vCPUs"; the strong ALM-NGOpt baseline "ran on a 96 vCPU machine for 18 hours (geometric problem) and 34 hours (simple-to-build)". — [arXiv:2506.19583](https://arxiv.org/abs/2506.19583)
- **coilstellaration** (Proxima, new in 2026): MIT, coil-set benchmark built on constellaration with its own leaderboard; needs gfortran, CMake, REGCOIL binary; Docker example asks for `--shm-size=16g`. Last commit 2026-11-16. [FH from clone] — [coilstellaration repo](https://github.com/proximafusion/coilstellaration)
- **QUASR**: "almost 370,000 quasisaxisymmetry and quasihelically symmetric devices along with coil sets", vacuum fields, archived on Zenodo and browsable at quasr.flatironinstitute.org. [FH abstract] — [arXiv:2409.04826](https://arxiv.org/abs/2409.04826). Formats are SIMSOPT and VMEC. [SS] — [Zenodo record](https://zenodo.org/records/13717741)
- **High-order near-axis codes**: StellaratorNearAxis.jl (Julia; Ruth 2024) [SS] and pyAQSC, "the first code for solving the near-axis expansion of anisotropic-pressure quasisymmetric equilibria to any order", auto-differentiable, with a 6th-order demonstration [FH abstract]. — [arXiv:2411.04352](https://arxiv.org/abs/2411.04352); [arXiv:2505.20475](https://arxiv.org/abs/2505.20475)
- A fast-ion benchmark dataset "BoSS" is mentioned in commentary only. [SS] — [blog commentary](https://blog.pebblous.ai/report/boss-stellarator-fast-ion-benchmark/en/)

### Inferences
- A near-axis parameter survey of millions of configurations is a few hours of laptop time; compute is not the constraint for near-axis work.
- ConStellaration is feasible on 4 cores only in the "slow baseline" sense (tens of hours per run), against reference baselines produced on 96 vCPUs.
- pyQSC's packaging defect (no declared dependencies) and its 2024 last commit suggest fixes there may not be reviewed quickly.

### Gaps
- QUASR's licence: not stated in any source I could read.
- DESC and SIMSOPT were not installed or timed; ConStellaration's forward model was not installed.
- Current ConStellaration leaderboard scores and number of entrants: not readable from here.

---

## 2. What the field has explicitly asked for, and which items need only a laptop

### Takeaway
The clearest explicit "open questions" list of 2026 comes from the pyQSC author himself (Landreman, 22 Sep 2026) and names exact quasisymmetry existence as open. Most other "future work" items are things the authors say they will do themselves; I could not read the issue trackers.

### Cited Findings
- Landreman, "Analytic toroidal 3D MHD equilibria…" (arXiv 22 Sep 2026), conclusion [FH, quoted]: "Since at least three families of exact solutions exist, it is likely that there are others. Similar analytic solutions could be sought with other numbers of field periods, or with ι having a constant value that is irrational… It remains unknown whether exact analytic solutions exist for non-axisymmetric toroidal magnetic fields with nested toroidal flux surfaces if we demand the fields be curl-free (vacuum) or force-free. Finally, a related open question is the existence of non-axisymmetric equilibria satisfying exact quasisymmetry." He also says the solutions "will be practically useful as tests of numerical equilibrium codes". — [arXiv:2609.26742](https://arxiv.org/abs/2609.26742)
- The same paper's acknowledgement [FH, quoted]: "These solutions were discovered using the artificial intelligence model GPT-6 Astra Pro, which was also used to draft parts of the manuscript. All equations were confirmed manually by the author." — [arXiv:2609.26742](https://arxiv.org/abs/2609.26742)
- Items on that list are being taken within weeks by specialists: Issan, Citrin, Landreman, Tracey broke stellarator symmetry (5 Oct 2026) [FH abstract]; de Lacy and Sato gave a real-analytic equilibrium with trivial Euclidean symmetry group and rotational transform 2 (12 Nov 2026), "discovered with the assistance of GPT-6 Astra Ultra" [FH]. — [arXiv:2610.07304](https://arxiv.org/abs/2610.07304); [arXiv:2611.16688](https://arxiv.org/abs/2611.16688)
- Ruth, Jorge, Bindel (near-axis divergence) [FH]: the expansion "diverges in vacuum due to ill-posedness" and can be regularised; "a moderate-order near-axis expansion (say, 4 ⩽ N ⩽ 6) could be used to explore the space of stellarators more effectively"; whether regularisation works with pressure "remains to be seen", with "small denominators for near-rational ι0" as the obstacle. — [arXiv:2411.04352](https://arxiv.org/abs/2411.04352)
- pyAQSC authors [FH]: "In the future, we will perform parameter space studies"; the expansion "is only convergent within a small volume around the axis". — [arXiv:2505.20475](https://arxiv.org/abs/2505.20475)
- QUASR paper [FH]: extending the workflow to knotted devices "is the subject of future work"; the paper already does the PCA / landscape analysis of the database itself. — [arXiv:2409.04826](https://arxiv.org/abs/2409.04826)
- ConStellaration limitations [FH]: degree of QI in the dataset is limited by the sampling budget; dataset covers plasma boundaries only, not coils. The benchmark exists "to lower the entry barrier for optimization and machine learning researchers". — [arXiv:2506.19583](https://arxiv.org/abs/2506.19583)
- Landreman–Paul 2021: quasisymmetry-breaking mode amplitudes "as small as the typical ∼50 μT geomagnetic field" for a 1 T field at aspect ratio 6. [FH abstract] — [arXiv:2108.03711](https://arxiv.org/abs/2108.03711)

### Inferences
- Laptop-feasible and requested: (a) searching for further exact solution families on Landreman's list; (b) leaderboard submissions to ConStellaration / coilstellaration; (c) moderate-order near-axis surveys. Of these, only (a) involves a statement that could be machine-checked.
- Laptop-feasible but already done by the authors or promised by them: QUASR dataset analysis, pyAQSC parameter studies, reproducing Landreman–Paul (inputs are shipped as test files, so a reproduction adds nothing new).
- A coil-complexity versus quasisymmetry-error study is what QUASR and coilstellaration already are, at a scale (370,000 devices) a 4-core machine cannot match.

### Gaps
- Open-issue lists for DESC, SIMSOPT, pyQSC, pyQIC, VMEC++ could not be read (GitHub API blocked). Whether maintainers have posted specific requests there is unknown.

---

## 3. Is there a computer-algebra or formal-verification angle the field would value?

### Takeaway
Until September 2026 I find no sign that plasma physics asked for formal methods. Since then, exact 3D MHD equilibria have become an active Lean + AI topic led by professional analysts, and one of those papers explicitly says computer-assisted (rigorous-numerics) approaches "should have immense value". The angle is real, and it is already occupied by specialists.

### Cited Findings
- Gómez-Serrano, Liehr, Taylor, "Counterexamples to Grad's conjecture" (arXiv 21 Sep 2026, 147 pages) [FH]: smooth MHS equilibria on solid tori with nested pressure tori and exactly C_N symmetry; "A Lean 4 certification of the above results is also provided." — [arXiv:2609.24739](https://arxiv.org/abs/2609.24739)
- Their Lean repository [FH from clone]: 5,847 `.lean` files, 566,654 lines, pinned to Lean 4.33.1, last commit 2026-09-22, no licence file at the root. `Showcase.lean` is 242 lines, "imports only Mathlib", and defines the MHS system directly (`IsMHS`: `cross (B x) (curl B x) + gradient P x = 0 ∧ divergence B x = 0` in the interior, B tangent on the boundary), nested pressure tori, and the symmetry group. — [Grad-Conjecture repo](https://github.com/lukasliehr/Grad-Conjecture)
- Their LLM disclosure [FH, quoted]: work "began in July 2026. In July, the authors set up the problem and provided a detailed roadmap… they used a combination of the LLMs GPT-5.6 Sol, Claude Fable 5, and Claude Opus 5 to fill in the technical details… The Lean code was produced by the same models with continuous guidance from the authors", later switching to "GPT-6 Astra and Claude Fable 5.1". — [arXiv:2609.24739](https://arxiv.org/abs/2609.24739)
- Same paper on rigorous numerics [FH, quoted]: "we also ran an independent computer-assisted approach. This latter attempt was abandoned once it became clear that the approach we present here would close, but it would nonetheless be interesting to systematically pursue such methods as they should have immense value." — [arXiv:2609.24739](https://arxiv.org/abs/2609.24739)
- Same paper on why this matters to the field [FH, quoted]: "The existence of three-dimensional equilibria with nested pressure surfaces is the central mathematical question of stellarator theory"; VMEC and DESC "assume nested flux surfaces from the outset". — [arXiv:2609.24739](https://arxiv.org/abs/2609.24739)
- Their equilibria are not stellarator-like: the field "vanishes exactly on a round magnetic axis". [FH abstract] Landreman's and de Lacy–Sato's have nowhere-vanishing fields. [FH] — [arXiv:2609.26742](https://arxiv.org/abs/2609.26742); [arXiv:2611.16688](https://arxiv.org/abs/2611.16688)
- Quasisymmetry is not formalised anywhere I looked: a text search for "quasisymmetr", "omnigen", "isodynamic" returns nothing in google-deepmind/formal-conjectures (24 Nov 2026 snapshot), and the Grad repo mentions quasisymmetry only in two bibliography entries. [FH] — [formal-conjectures repo](https://github.com/google-deepmind/formal-conjectures)
- Published errata in hand-derived near-axis papers exist: Landreman–Sengupta (Phys. Plasmas 2021) "wrong sign preceding the summation" in Eq. (11); Rodríguez (JPP 2024) several sign typos; an erratum to Plunk, Landreman, Helander (2019) fixing normalisations. In each case "numerical results and main conclusions are unaffected" or similar. [SS] — [JPP erratum](https://www.cambridge.org/core/product/7C7D1CCD9A8B580F7643C2FB30EF6383); [Princeton record of erratum](https://collaborate.princeton.edu/en/publications/erratum-solving-the-problem-of-overdetermination-of-quasisymmetri/)
- De Lacy–Sato [FH, quoted]: near-axis expansions "are formal and are truncated at low order in applications. When quasisymmetry is imposed they become overdetermined at third order [GB91a]"; Weitzner's all-order expansions come without "a proof of convergence". — [arXiv:2611.16688](https://arxiv.org/abs/2611.16688)
- Negative result: a search for Lean or formal verification of plasma control, power-grid or medical-device software in 2026 returned only generic Lean work (TorchLean, a Rust-to-Lean pipeline for cryptographic code, a self-described proof-of-concept hardware monitor). [SS] — [TorchLean, arXiv:2602.22631](https://arxiv.org/abs/2602.22631)

### Inferences
- "Verified derivation of the near-axis equations to high order" has weak demand: errata are sign/typo level and did not change results; codes already exist to arbitrary order; the series is known to be divergent/asymptotic, so high-order exactness has limited physical payoff.
- A Lean statement of exact quasisymmetry could be built on the existing `IsMHS` definitions. A coordinate-free form (a divergence-free vector field u with B × u = ∇ψ and u·∇|B| = 0, following the Burby–Kallinikos–MacKay line) uses only derivatives Mathlib already has. I did not write or check such a statement; this is my inference from reading `Showcase.lean` and the abstract of [arXiv:1912.06468](https://arxiv.org/abs/1912.06468).
- Checking a claimed explicit closed-form equilibrium is an elementary-function identity. That is a far smaller formalisation job than the 566,654-line Nash–Moser development.

### Gaps
- Whether anyone has started a Lean certification of Landreman's closed-form families: not found, not ruled out.
- I found no plasma-physics source before September 2026 asking for formal methods; absence of evidence from a handful of searches only.

---

## 4. Counterarguments, and the full-ambition question

### Takeaway
The 2026 record shows that AI + Lean can settle famous open statements, and that public subscription models were enough when an expert steered them in this exact subfield. It also shows that the big results came from lab staff with unreleased models, that autonomous runs at a few hundred dollars per problem solve only a few percent of hard open problems, and that the stellarator question is now being worked by its leading specialists with the same tools. Exact quasisymmetry existence is a legitimate open target with a clean verification path, low odds for an outsider, and small practical consequence.

### Cited Findings: what the 2026 results actually were
- Anthropic, 10 Aug 2026 [WF]: "an unreleased research version of Claude" raised the proven share of zeta zeros on the critical line from 41.6% to 67.2%. A non-mathematician staff member (Jarred Sumner) asked it to "take a real stab" at RH; about 60 subagents ran for roughly a day and a half; about 31 million output tokens across two Claude Code sessions; Anthropic mathematicians Alpöge and Furman validated it; Eric Easley worked with Claude on the Lean formalisation. Anthropic said it does not expect these techniques to prove RH. — [Anthropic research note](https://www.anthropic.com/research/riemann-zeta); paper [FH abstract]: "Proof discovered autonomously by Claude (Anthropic); verified and communicated by the listed authors" — [arXiv:2608.13637](https://arxiv.org/abs/2608.13637)
- Alpöge's RH proof, released 20 Nov 2026 [SS]: a Lean library of about 175,187 files / 17 million lines developed with Claude; the Lean FRO and Mathlib Initiative issued a joint statement on 21 Nov certifying statement and proof of the classical RH; the generalized claim is outside that statement; "no human can yet claim to understand the argument"; a full end-to-end comparator run was still listed as in progress on 22–23 Nov. — [Lean on X](https://x.com/leanprover/status/2124159121645801562); [Nature news](https://www.nature.com/articles/d41586-026-03812-7); [Scientific American](https://www.scientificamerican.com/article/anthropic-mathematician-claims-to-prove-riemann-hypothesis-most-famous-problem-in-all-of-math/); cautionary framing in [Unite.ai](https://www.unite.ai/levent-alpoges-unverified-riemann-claim-puts-claudes-lean-code-in-focus/)
- First-hand corroboration [FH]: in the formal-conjectures snapshot of 24 Nov 2026, `Millennium/RiemannHypothesis.lean` tags `riemannHypothesis` and the Dirichlet-L-function GRH as `research solved` with `formal_proof … at "https://github.com/leanprover-community/riemann-hypothesis"`; `git ls-remote` confirms that repository and `anthropics/zeta-23-lean` exist. I did not build or check the proofs. — [formal-conjectures repo](https://github.com/google-deepmind/formal-conjectures)
- OpenAI, 6 Oct 2026 [SS]: 722 manuscripts in 372 result families from "an internal model that OpenAI has not released"; average effort "the equivalent of just 3 hours of ChatGPT Pro 'thinking'" per result from a pool of about 4,000 problems; 300 of 719 remaining papers formalised in Lean; three to five withdrawn (sources conflict) and 14 revised. The Lean-verified quasi-RH version places zeros between 1/8 and 7/8. — [Quanta](https://www.quantamagazine.org/ai-proved-the-quasi-riemann-hypothesis-what-happens-now-20261019/); [Science](https://www.science.org/content/article/openai-just-changed-math-forever-exactly-how-no-one-yet-knows)
- An OpenAI claim was formalised independently afterwards: Cârstea and Oksanen announce a Lean 4 proof of a special case of the anisotropic Calderón uniqueness theorem (7 Nov 2026); "OpenAI's release did not include a Lean proof". [FH abstract] — [arXiv:2611.10537](https://arxiv.org/abs/2611.10537)
- Individuals with public tools:
  - Liam Price, described as a 23-year-old amateur, solved Erdős #1196 by prompting GPT-5.4 Pro in April 2026; validated by Tao, later formalised in Lean; sources disagree on whether his plan was Plus or Pro; Lichtman called the raw output "quite poor". [SS] — [Scientific American](https://www.scientificamerican.com/article/amateur-armed-with-chatgpt-vibe-maths-a-60-year-old-problem/)
  - K. M. Knausgård, single author, improved Claude's zeta constants (67.35% simple and on the line) using "OpenAI Codex and Anthropic Claude under the author's direction"; the Lean proofs "use three axioms beyond those of Mathlib… Each records that a search program returned true"; "an experiment in AI-assisted mathematical research". [FH] His background was not checked. — [arXiv:2610.08965](https://arxiv.org/abs/2610.08965)
- Systematic measurements:
  - Google DeepMind (21 May 2026): best agent "autonomously resolved 9 of 353 open Erdős problems at the per-problem cost of a few hundred dollars, proved 44/492 OEIS conjectures". [FH abstract] — [arXiv:2605.22763](https://arxiv.org/abs/2605.22763)
  - FrontierMath Erdős (6 Sep 2026): 68 hard open Erdős problems, five AIs, "$300 per problem. One (GPT-6 Astra) scored 3%, and all others scored 0%." [FH abstract] — [arXiv:2609.25050](https://arxiv.org/abs/2609.25050)
- The user's recollection of an individual getting "70% of the way" to GRH/RH on a $100/month OpenAI plan: **not found**. Closest matches are the Anthropic 67.2% result (lab, unreleased model), OpenAI's 7/8 zero-free half-plane (lab, internal model), and Knausgård's 67.35% (individual, Codex + Claude). [SS/FH as above]

### Cited Findings: objections specific to the stellarator target
- Design practice has moved toward quasi-isodynamic fields. Type One Energy's lead theorist: "we couldn't make quasi-symmetry work as well as we thought we could… it looks like QI is the prominent optimization choice". Proxima is QI. [SS] — [IEEE Spectrum](https://spectrum.ieee.org/stellarator). ConStellaration calls QI "a promising path to commercial fusion" and its benchmarks are all QI. [FH] — [arXiv:2506.19583](https://arxiv.org/abs/2506.19583)
- Numerical optimisation already reaches quasisymmetry error at geomagnetic-field level over a volume (Landreman–Paul). [FH abstract] — [arXiv:2108.03711](https://arxiv.org/abs/2108.03711)
- The standard heuristic points to non-existence: "The Garren-Boozer (GB) conundrum has been interpreted to mean that globally quasisymmetric magnetohydrostatic (MHS) equilibria, other than axisymmetric solutions, with isotropic pressure do not exist… the governing equations become overdetermined at the 3rd order." [FH abstract] — [arXiv:2505.20475](https://arxiv.org/abs/2505.20475)
- Existing rigorous partial results need relaxations: nearly quasisymmetric equilibria "with a small force" (Constantin–Drivas–Ginsberg) [FH abstract]; "weakly quasisymmetric" fields in asymmetric domains (Sato) [SS]. — [arXiv:2009.08860](https://arxiv.org/abs/2009.08860); [arXiv:2204.09241](https://arxiv.org/abs/2204.09241)
- No 2026 resolution of exact isotropic-pressure quasisymmetry was found. [SS, one search] Landreman listed it as open on 22 Sep 2026. [FH] — [arXiv:2609.26742](https://arxiv.org/abs/2609.26742)
- Crowding and speed: four exact-equilibrium papers from at least three groups (Gómez-Serrano et al.; Landreman; Issan/Citrin/Landreman/Tracey; de Lacy–Sato) appeared between 21 Sep and 12 Nov 2026. Three of the four state frontier-model use; I did not check the Issan et al. paper for an AI statement. [FH] — [arXiv:2609.24739](https://arxiv.org/abs/2609.24739); [arXiv:2610.07304](https://arxiv.org/abs/2610.07304); [arXiv:2611.16688](https://arxiv.org/abs/2611.16688)
- Reception of outsider, AI-assisted work: arXiv extended tightened endorsement (institutional email no longer sufficient) to all categories on 21 Jan 2026, and in May 2026 announced a one-year ban for unchecked AI content. [SS] — [The Next Web](https://thenextweb.com/news/arxiv-ai-slop-ban-researchers-preprint); [TechCrunch](https://techcrunch.com/2026/05/16/research-repository-arxiv-will-ban-authors-for-a-year-if-they-let-ai-do-all-the-work). erdosproblems.com stopped displaying problem statuses and solved counts on 6 Oct 2026, and Tao's AI-contributions wiki "is no longer updated as of June 30, 2026". [FH from repo README] — [teorth/erdosproblems](https://github.com/teorth/erdosproblems)
- Plausible-but-wrong output is documented even at lab scale: OpenAI withdrew papers over "a shared sign error", and "nobody has checked whether the Lean statements say what the paper claims" for parts of its release. [SS] — [Science](https://www.science.org/content/article/openai-just-changed-math-forever-exactly-how-no-one-yet-knows)

### Inferences
- **Does the evidence overturn "a non-specialist can only do small things"?** Partly. Named commercial models found exact stellarator-relevant equilibria (GPT-6 Astra Pro, Astra Ultra) and wrote half a million lines of Lean that compiles to a Mathlib-only statement (Claude Fable 5 / 5.1, GPT-5.6 / 6). These appear to be released products, so the capability is plausibly available on paid plans; I did not verify which plan tier or usage volume each paper needed. In every stellarator case an expert chose the ansatz, supplied the roadmap or checked every equation. The two non-specialist successes (Price; Sumner) involved a combinatorics problem and an unreleased model with 60 subagents respectively.
- **Does it transfer to one person on a consumer plan?** The FrontierMath Erdős figure (0–3% at $300 per problem, autonomous) is the best available base rate for hard open problems with public models as of September 2026. Models have moved since (the RH proof was released about eleven weeks later), so the base rate is probably higher now, by an amount I cannot source.
- **Exact quasisymmetry as a full-ambition target:** (i) explicitly open and named by the field's leading tool author; (ii) expressible on top of Mathlib definitions that already exist in a public repo; (iii) peacetime; (iv) a positive answer by explicit example would be checkable by any computer-algebra system and by Lean. Against: (v) the expected answer is non-existence, which needs a rigidity proof, not a search; (vi) the people best placed to win are already searching with the same or stronger models; (vii) the practical consequence is small, because designers use approximate quasisymmetry or QI either way. Landreman's own word for the value of exact existence results is "reassuring", plus benchmark utility.
- **Has numerical optimisation made the existence question practically moot?** For device design, largely yes. For code verification, no: exact solutions are wanted as tests of VMEC/DESC-class codes.
- The prerequisite-knowledge objection (Boozer coordinates, MHD equilibrium theory) weakens for a pure "find an explicit solution and verify the PDE" task and stays strong for anything involving physics judgement (which objective matters, what a designer would use).

### Gaps
- Which Claude model Alpöge used, the compute spent, and whether a comparable run is possible on a subscription plan: not found.
- Journal of Plasma Physics policy on AI-assisted or outsider submissions: not researched.
- arXiv policy claims are from news snippets; arXiv's own pages were not read.

---

## 5. Alternatives, including "open problem already expressible in Mathlib, peacetime, with practical consequences"

### Takeaway
I could not find a problem that satisfies all four wishes at once (open, already stated in Lean, peacetime, large practical consequence). In the largest public collection of Lean-stated open problems, applied areas are nearly empty: zero open statements in fluid mechanics, optimisation/operations research, control, numerical analysis, statistical mechanics or biology. The stellarator exact-solution question is unusually close to the user's criteria by comparison.

### Cited Findings
- formal-conjectures (Google DeepMind), snapshot 24 Nov 2026 [FH, my own tally of `@[category …]` attributes]: 1,952 `research open` statements, 3,142 `research solved`. Open statements by AMS class: number theory (11) 1,297; combinatorics (05) 528; convex/discrete geometry (52) 90; quantum theory (81) 44; information/coding (94) 22; PDE (35) 8 statements in 2 files; calculus of variations/optimal control (49) 1; optics/electromagnetism (78) 2. **Classes 65, 76, 82, 90, 92, 93 have none.** — [formal-conjectures repo](https://github.com/google-deepmind/formal-conjectures). The project's paper (May 2026) reported 1,642 open conjectures and says the benchmark "has already been leveraged for… the resolution of open research conjectures". [FH abstract] — [arXiv:2605.13171](https://arxiv.org/abs/2605.13171)
- The repo invites new statements: "consider adding (or even just opening an issue describing) your favourite conjecture", and "Follow the same AI usage conventions as Mathlib". [FH] — [formal-conjectures repo](https://github.com/google-deepmind/formal-conjectures)
- Open, Lean-stated, applied-leaning items in that repo [FH; file names and tags read directly, practical relevance is my gloss]:
  - `Millennium/NavierStokes.lean`: existence and smoothness on R³ and periodic, open.
  - `Paper/DeGiorgi.lean`: De Giorgi's conjecture for Δu + u − u³ = 0, "4 ≤ n ≤ 8 remains open" (each dimension listed separately "to enable partial solutions").
  - `Wikipedia/SquarePacking.lean`: four open finite packing cases (11 and 17 squares in a square; 21 circles in a square; 15 circles in a circle).
  - `Wikipedia/SnakeInTheBox.lean`: longest snake in dimension 9, open. `Wikipedia/Hadamard.lean`: Hadamard conjecture and order 668, open. `Wikipedia/SteinerSystem.lean`: one open. `Wikipedia/CernyConjecture.lean`: synchronising automata, open. `Arxiv/0805.3777/TensorRank.lean`: two open generic/maximal rank statements.
  - `OptimizationConstants/` (5 files): Tao's optimisation constants, e.g. the Sidon autocorrelation constant C_1a, with 2026 AI-driven bound improvements cited in the file header.
  - **Flag for exclusion (security/weapons overlap):** `Wikipedia/IntegerFactorization.lean` (cryptography), `Arxiv/2610.08796/BayesianSecurityGames.lean` (security games), `Arxiv/2104.00502/BarkerSequence.lean` and `Wikipedia/SparseRuler.lean` (pulse-compression / antenna-array uses include radar), `OpenQuantumProblems/` (quantum information).
  - `Arxiv/2607.27197/MaxwellConjecture.lean` records that the Maxwell point-charge conjecture was disproved in 2026, so it is no longer open.
- Erdős problems, 24 Nov 2026 [FH from repo]: 1,221 problems; 339 proved, 144 disproved, 104 otherwise solved; **570 "completely open"**; 323 have a solution formalised in Lean; 1,024 statements are formalised in formal-conjectures; 114 carry prizes. — [teorth/erdosproblems](https://github.com/teorth/erdosproblems). Secondary reporting gives 318 Lean-formalised solutions in October 2026 "up from 38 at the end of 2025". [SS] — [FourWeekMBA](https://fourweekmba.com/ai-tao-math-now-hides-its-problems-after-openais-372-proofs/)
- Mathlib missing-theorem list, 25 Nov 2026 [FH, my parse of `docs/1000.yaml`; approximate]: 1,179 entries, 248 carrying a declaration or URL, so roughly 930 unformalised. Unformalised entries with engineering relevance include Routh–Hurwitz, Kharitonov, Nyquist–Shannon sampling, Shannon source coding, Shannon–Hartley, Lax–Richtmyer, Lax–Wendroff, Perron–Frobenius, Gauss–Markov, Berry–Esséen, Glivenko–Cantelli, **Fisher–Tippett–Gnedenko** (extreme value theory), Cramér's large-deviations theorem, Nash–Moser. `Mathlib/Probability/Distributions/` currently has Bernoulli, Beta, Binomial, Cauchy, Exponential, Gamma, Gaussian, Geometric, Pareto, Poisson, Uniform; no Weibull, Gumbel, Fréchet or generalised extreme value file. — [mathlib4 repo](https://github.com/leanprover-community/mathlib4)
- Mathlib's attitude to AI-written code [SS]: Johan Commelin (Feb 2026, personal view): "there is no point in risking Mathlib's carefully curated high quality status by mass-adopting AI-generated code"; a disclosure policy for AI use in PR descriptions was under discussion in January 2026. — [Zulip archive](https://leanprover-community.github.io/archive/stream/287929-mathlib4/topic/policy.20proposal.3A.20mention.20AI.20use.20in.20PR.20description.3F.html)
- Energy benchmarks: ARPA-E GO Competition Challenge 3 selected teams in December 2022 and a Georgia Tech item dated 8 Apr 2026 reports a second-place finish; no open public leaderboard or L2RPN 2026 round was found. [SS] — [GO Competition news](https://gocompetition.energy.gov/news/arpa-e-announces-challenge-3-teams)

### Assessment table (inference, built from the findings above)

| Candidate | What is open | Stated in Lean today? | Who would use a result | How verified | Crowding | Full-ambition plausibility for AI + Lean on a plan | Weapons/crypto flag |
|---|---|---|---|---|---|---|---|
| Exact stellarator equilibria: new explicit families; vacuum/force-free with nested surfaces; exact quasisymmetry | Named open by Landreman, Sep 2026 | MHS yes (Grad repo); quasisymmetry no | Equilibrium-code developers (test cases); stellarator theorists | CAS identity + Lean; expert read | High since Sep 2026, specialists with frontier models | New family: possible. Vacuum/force-free: unknown. Exact QS: low, likely needs non-existence proof | None found |
| De Giorgi conjecture, n = 4…8 | Open | Yes | PDE analysts; indirect link to phase-transition models | Lean | Not assessed | Low; deep regularity theory, no search shortcut | None |
| Navier–Stokes regularity | Open | Yes | Mathematics; practical effect of a proof is indirect | Lean | Not assessed | Very low | None |
| Finite packing cases (squares/circles) | Open | Yes | Mostly mathematicians; packing practice does not wait on proofs | Lean + interval arithmetic | Not assessed | Moderate for a single case; compute-bound | None |
| Erdős problems | 570 open | 1,024 statements | Combinatorialists | Lean; site community | Very high; labs sweep them | Proven route for amateurs (Price), few % per attempt | None, but little practical use |
| Tao optimisation constants | Bounds improvable | 5 files | Analysts | Lean / certificates | High | Incremental bound gains realistic | None |
| Mathlib: Fisher–Tippett–Gnedenko, Weibull/Gumbel/Fréchet, Routh–Hurwitz, Nyquist–Shannon | Not open problems; missing library theorems | n/a | Mathlib users; later verified-engineering work | Lean + maintainer review | Low to moderate | Near-certain to be provable; merge is the hard part | None |
| ConStellaration / coilstellaration leaderboards | Better QI designs | n/a (not Lean) | Proxima Fusion and benchmark users | Leaderboard re-runs VMEC++ | Unknown | n/a; compute-limited on 4 cores | None |
| Verified control / dosing / structural software in Lean | No stated demand found | n/a | None identified | — | — | — | — |

### Inferences
- The nearest thing to "practical + open + Lean-ready" outside the stellarator question is finite packing and coding-design existence problems (packing cases, snake-in-the-box, Hadamard, Steiner systems). Their practical consequence is modest and several are compute-bound search problems.
- Formalising the *statements* of open problems from control, optimisation, fluids or numerical analysis would fill a measured gap (zero entries) in a benchmark that explicitly asks for contributions. That is a contribution with an identifiable user (the benchmark and every lab that runs on it) even before any proof.
- The extreme-value-theory items connect directly to the earlier rogue-wave work and have a definite user, but they are not open problems, so they do not meet the user's stated ambition.

### Gaps
- I did not find or assess a curated list of open problems in control theory or power systems with precise mathematical statements; candidates there are unassessed, not absent.
- Verified medical-dosing and structural-engineering calculations: no source found either way.
- Leaderboard participation numbers for ConStellaration and coilstellaration: unknown.

---

## 6. Best options overall and a first 4-to-6-week milestone for each

### Takeaway
Three options survive: (A) a time-boxed attempt on Landreman's list of open exact-equilibrium questions, with exact quasisymmetry as the stretch goal; (B) adding and attacking Lean statements of applied open problems via formal-conjectures; (C) extreme value theory in Mathlib as the low-risk fallback. Each has a check that someone other than the author performs.

### Cited Findings
- Landreman's open list and benchmark use of exact solutions. [FH] — [arXiv:2609.26742](https://arxiv.org/abs/2609.26742)
- A public Mathlib-only Lean statement of MHS equilibria exists to build on, and the Grad repo ships a `Comparator/` directory that checks statement correspondence and that only the three standard axioms are used. [FH] — [Grad-Conjecture repo](https://github.com/lukasliehr/Grad-Conjecture)
- "Comparator" checks are the accepted way to show a Lean proof matches a trusted statement: used for the Anthropic zeta result [WF], the RH certification [SS], and the Calderón formalisation [FH abstract]. — [Anthropic research note](https://www.anthropic.com/research/riemann-zeta); [arXiv:2611.10537](https://arxiv.org/abs/2611.10537)
- formal-conjectures accepts new statements by pull request with human review. [FH] — [formal-conjectures repo](https://github.com/google-deepmind/formal-conjectures)
- Extreme-value distributions and Fisher–Tippett–Gnedenko are absent from Mathlib as of 25 Nov 2026. [FH] — [mathlib4 repo](https://github.com/leanprover-community/mathlib4)

### Inferences: proposed milestones
**Option A. Exact 3D MHD equilibria (narrowed stellarator option).**
- Weeks 1–2: re-derive one Landreman family and confirm with a computer-algebra script that (∇×B)×B = ∇p and ∇·B = 0 hold identically; then state and prove the same in Lean on top of the public `IsMHS`-style definitions. This tests whether the person + Claude can do the verification half at all.
- Weeks 3–6: run a structured search for a family not yet published (constant irrational ι; other field-period numbers; then curl-free or force-free with nested surfaces). Write the quasisymmetry condition as a Lean definition and test every candidate against it.
- External check: the symbolic residual script and the Lean file run on someone else's machine with standard axioms only; the solution is run as a test case through VMEC++ (9.5 s for a W7-X-size case here); the result is sent to the authors who posed the question.
- Honest odds: the verification milestone is likely achievable; a genuinely new family is uncertain because specialists took three such steps in seven weeks; exact quasisymmetry is a long shot and may be false.
- Model caveat: the published discoveries used GPT-6 Astra Pro/Ultra; Claude Fable 5/5.1 is documented for proof detail and Lean code in the Gómez-Serrano paper, not for the discovery step.

**Option B. Lean statements of applied open problems, then attack one.**
- Weeks 1–3: choose three to five open problems from control, optimisation, coding/packing or fluid stability with no security or weapons use; write Lean statements; submit to formal-conjectures.
- Weeks 4–6: attempt the most tractable one (a finite packing or design-existence case is the realistic choice) with a Lean proof.
- External check: maintainers' review and merge of the statements; a proof checked by comparator against the merged statement.
- Honest odds: statement merges are likely; resolving one is uncertain and for packing cases depends on compute.

**Option C. Extreme value theory in Mathlib (fallback, not an open problem).**
- Weeks 1–6: Gumbel, Fréchet, Weibull distributions with basic API, then a first Fisher–Tippett–Gnedenko statement.
- External check: Mathlib maintainer review; merged PRs.
- Honest odds: proofs are near-certain; merge depends on meeting Mathlib style, where AI-written code meets documented scepticism.

**Not recommended as the main project:** a general "stellarator quasisymmetry design" effort with SIMSOPT/DESC. Lean has no role in those workflows, the dataset-scale studies are already done by funded groups, and commercial designs are moving to QI.

### Gaps
- No milestone above has been tried; feasibility of the Lean half of Option A (elementary-function identities in three variables) is untested here.
- Whether Landreman's supplemental materials already include machine-checkable verification scripts was not checked; if they do, weeks 1–2 of Option A reproduce existing work.
