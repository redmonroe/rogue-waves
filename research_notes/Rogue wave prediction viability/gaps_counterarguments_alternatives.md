# Rogue wave probability prediction: open gaps, counterarguments and alternative maritime targets (state of knowledge 26 Nov 2026)

Method note (applies to every section). All web material was collected on 2026-11-26 through web search (standard and extended). Direct page fetches failed for every publisher host tried (os.copernicus.org, nhess.copernicus.org, www.nature.com, ebi.ac.uk/Europe PMC: DNS or proxy refusal), so **no paper was read in full text in this session**. Web claims below rest on search-result summaries and snippets of the cited pages and should be re-checked against the full text before being quoted. Claims marked "(local)" were measured directly on the files in `/home/user/rogue-waves/code` (the paper's public code) and `/home/user/rogue-waves/data` (the 157 downloaded parquet files) and are first-hand.

## 1. Limitations that the authors and reviewers state for Häfner et al. 2021/2023 and FOWD

### Takeaway
The authors themselves flag buoy sensor bias, the 100-wave aggregation, causal-analysis assumptions, and unphysical behaviour of the fitted formula for narrow directional spread; all validation published so far is inside the CDIP network. I found no independent test of the 2023 closed-form model on another observing system as of Nov 2026.

### Cited Findings
Stated by the authors (PNAS 2023 / arXiv 2311.12579, published Nov 2023):
- Training was "solely on observations from buoys"; the paper says systematic sensor bias is common in buoys, can produce spurious causal relationships, and may hurt generalisation to other sensors. — [PNAS](https://www.pnas.org/doi/10.1073/pnas.2306275120); [arXiv](https://arxiv.org/pdf/2311.12579)
- Aggregating into 100-wave chunks "underestimate[s] the per-wave rogue wave probability in sea states in which rogue waves do not occur independently of each other" (e.g. strong group structure); the per-wave conversion assumes independence. The authors say these limitations could reduce the ability to detect causal pathways and underestimate true risk. — [PMC copy](https://pmc.ncbi.nlm.nih.gov/articles/PMC10691345/)
- The fitted directional-spread term diverges as spread goes to zero; the model only saw real seas with spread above about 0.2, and the authors propose a non-diverging replacement. — [arXiv](https://arxiv.org/pdf/2311.12579)
- The rogue criterion is crest-to-trough height H/Hs > 2.0, not crest height; the "Rayleigh is an upper bound" statement is made for crest-to-trough waves. In the most extreme conditions cited (spread 13°, steepness 0.008, r = 0.88) the model reaches the Rayleigh value p = 3.3e-4 per wave. — [PNAS](https://www.pnas.org/doi/full/10.1073/pnas.2306275120)
- Data filter: water depth 10 m to 4000 m and significant wave height at least 1 m. — [arXiv](https://arxiv.org/pdf/2311.12579)
- Crossing seas appear only as a possible pathway in the causal graph (high crossing angle influencing breaking), not as a dedicated analysis. — [arXiv](https://arxiv.org/pdf/2311.12579)
- Generalisation claims ("can generalize to unseen physical regimes, is robust to distributional shift", "outperforms existing wave theory on unseen stations from our dataset") are the authors' own, tested on held-out parts of the same dataset. — [arXiv](https://arxiv.org/html/2311.12579.pdf)
- The authors note crest-trough correlation r "is still lacking mainstream adoption" and "is not part of ECMWF's operational forecast" (as of 2023). — [arXiv](https://arxiv.org/pdf/2311.12579)

Stated in the 2021 papers:
- Sci Rep (12 May 2021): kurtosis, steepness and Benjamin-Feir index are "weak predictors"; Rayleigh is an upper bound in most conditions; the authors themselves doubt "whether the common rogue wave definition as any wave exceeding a certain height threshold is meaningful in practice". — [Sci Rep 2021](https://www.nature.com/articles/s41598-021-89359-1)
- Sci Rep 2021 applies a correction for 1.28 Hz discrete sampling described as conservative, maximum 4.2%; FOWD filtering already removes records with mean period below 5 s for 1.28 Hz CDIP data. — [Sci Rep 2021](https://www.nature.com/articles/s41598-021-89359-1)
- FOWD (JTECH 38(7), 2021): running-window QC with multiple flags; waves failing QC above 2 Hs and all waves above 2.5 Hs are logged to JSON rather than silently dropped; 38 deployments and 4 whole CDIP stations were blacklisted after visual inspection; about 67% of raw waves are filtered, mostly by the minimum-Hs requirement; filtered and unfiltered versions are both released. — [arXiv 2011.12071](https://arxiv.org/pdf/2011.12071); [JTECH](https://journals.ametsoc.org/view/journals/atot/38/7/JTECH-D-20-0185.1.xml)
- 80% of QC failures occur at only 5 CDIP locations (as reported in a search summary of a FOWD-using source). — [FOWD preprint copy](https://www.gfy.ku.dk/~mjochum/fowd.pdf)
- CDIP catalogue: about 750 years of surface elevation, almost all at 1.28 Hz. — [arXiv 2011.12071](https://arxiv.org/pdf/2011.12071)

First-hand facts from the public code and data (local, 2026-11-26):
- The 157 parquet files hold 14,985,244 aggregate rows and 106 columns, spanning 1993-09-02 to 2022-05-04. — [code repo](https://github.com/dionhaefner/rogue-wave-discovery) (local measurement)
- Geographic concentration by row count: US West Coast 71.0% (85 stations); a Southern California box (lat 32 to 38, lon -123.5 to -117, the paper's own "socal" subset definition) 52.2%; Hawaii 14.5%; western/central Pacific islands 6.2%; US East/Gulf coast 4.7% (40 stations but few rows); Alaska 0.6%; Caribbean 0.6%. Stations with median depth below 100 m hold 28% of rows; above 1000 m, 9%. — [code repo](https://github.com/dionhaefner/rogue-wave-discovery) (local measurement, before the paper's base constraints)
- 119,500 rows (0.80%) have a maximum H/Hs above 2.0. The dataset also carries `aggregate_100_max_rel_crest_height`; 21,848 rows (0.146%) exceed C/Hs > 1.25. The paper's label is `aggregate_100_max_rel_wave_height` with threshold 2.0, so the crest-height label is shipped but unused by the published model. — [constants.py in repo](https://github.com/dionhaefner/rogue-wave-discovery) (local)
- Hard-coded base constraints drop sea states with |skewness| > 0.5, |excess kurtosis| > 1 (code comment: "just no"), low-frequency drift energy fraction > 0.1, or valid-data ratio < 0.95. — [constants.py in repo](https://github.com/dionhaefner/rogue-wave-discovery) (local)
- Validation uses 26 held-out CDIP stations plus in-network subsets (socal, deep, shallow, summer, winter, Hs > 3 m, high/low frequency, low/high period, cnoidal, weakly nonlinear, narrow, wide). No non-CDIP data appear. — [constants.py in repo](https://github.com/dionhaefner/rogue-wave-discovery) (local)

Stated by later users and reviewers:
- Azevedo et al., Ocean Science 22, 367-386 (published 2 Feb 2026) used FOWD-CDIP for 2015-2021, found 81 rogue waves at four stations, and report that buoys record far more rogue waves than ERA5/ECMWF CY47R1 model estimates; they note FOWD rogue waves are at least 2 m high because of the Hs > 1 m filter, and that "residual differences occur where FOWD QC excludes short segments around spikes". — [Ocean Science](https://os.copernicus.org/articles/22/367/2026/)
- Same paper: buoys ride forward on crests and back in troughs, raising recorded mean level and underestimating the highest crests. — [Ocean Science](https://os.copernicus.org/articles/22/367/2026/)
- A referee of that paper objected that the model parameter used (expected maximum envelope height) "is not correlated with rogue wave height or occurrence", questioned comparing its 75% figure with the 0.72 predictive power of Cicon et al. 2024, and said nobody proposed a fixed r = 0.6 threshold. — [OS peer review](https://os.copernicus.org/articles/22/367/2026/os-22-367-2026-discussion.html)
- Häfner's PhD thesis ("An Ocean of Data", submitted 4 Apr 2022, defended 11 May 2022) closes with a discussion of implications for future rogue wave research and a "path to an improved rogue wave forecast". — [NBI thesis page](https://nbi.ku.dk/english/theses/phd-theses/dion-hfner/); [thesis PDF](https://nbi.ku.dk/english/theses/phd-theses/dion-hfner/PhD_Thesis_Dion_H_fner_2022.pdf)

### Inferences
- Candidate gaps a small effort could address with the public files alone: (a) calibration and reliability analysis of the closed-form model per station/region/depth rather than aggregate scores; (b) a crest-height (C/Hs > 1.25) version, since the label is already in the files; (c) sensitivity to the skewness/kurtosis/valid-ratio constraints and to the 100-wave independence assumption; (d) leave-one-region-out tests (e.g. train without the 52% Southern California block). None of these needs new data.
- An out-of-time test is possible in principle: the files end in May 2022, so about 4.5 years of newer CDIP displacement data exist that neither the authors nor (as far as I found) anyone else has scored against the model. It requires running the FOWD pipeline on raw CDIP files.
- The 40 East/Gulf-coast stations contribute under 5% of rows, so "157 stations" overstates geographic diversity; the effective sample is the California/Hawaii swell climate.

### Gaps
- I could not read the PNAS limitations section, the JTECH FOWD discussion, or the thesis final chapter in full; the specific "next steps" Häfner lists in the thesis are unverified.
- Whether the paper discusses crest height, shallow water below 10 m, or Lagrangian buoy motion as explicit limitations is unconfirmed.
- No published reviewer reports for the PNAS paper were found.
- No citing-paper sweep was possible; an independent validation may exist that the searches missed.

## 2. Other public datasets for independent validation: access terms and sizes

### Takeaway
Raw surface-elevation time series (needed to label individual rogue waves) are openly downloadable in bulk essentially only from CDIP; most other open archives publish bulk parameters or spectra, and fixed-platform laser/radar records are mostly industry-held. A truly independent sensor type is the hard part.

### Cited Findings
- CDIP: displacement records in netCDF4 via a THREDDS server; about 750 years total, almost all at 1.28 Hz; NOAA data catalogue lists per-deployment "sea surface displacement" datasets (e.g. a Datawell DWR-M3 near St. Augustine, FL, 2019, 1.280 Hz, about 22 m depth). — [arXiv 2011.12071](https://arxiv.org/pdf/2011.12071); [NOAA catalogue](https://data.noaa.gov/dataset/dataset/sea-surface-displacement-values-collected-in-situ-by-datawell-dwr-m3-directional-buoy-loc-30-00175)
- FOWD processing code is public and MIT-licensed (v0.5.2 archived on Zenodo); it accepts CDIP netCDF4 or a generic netCDF4 format and can run deployments in parallel; full FOWD has about 4.2 billion waves before filtering, 1.4 billion after. — [Zenodo](https://zenodo.org/records/4628203); [arXiv 2011.12071](https://arxiv.org/pdf/2011.12071)
- NDBC: its technical document states "transmission of raw time series data is not feasible with today's satellite limited message lengths"; the public archive lists spectral products (swden and directional coefficient files). — [NDBC Tech Doc 03-01](https://www.ndbc.noaa.gov/wavemeas.pdf); [NDBC data index](https://www.ndbc.noaa.gov/data/)
- Sofar Spotter Archive: free, public, hourly wave and inferred wind data, 2019 to March 2022, on AWS S3 (`s3://sofar-spotter-archive/`, no account needed), under the Sofar Data Access Agreement / Spotter Data Archive License (not CC-BY; research/non-commercial). — [AWS registry](https://registry.opendata.aws/sofar-spotter-archive/); [GitHub](https://github.com/sofarocean/spotter-data-archive)
- Copernicus Marine in-situ wave product (INSITU_GLO_WAV_DISCRETE_MY_013_045): integrated parameters ("TS" files) and spectra ("WS" files, since Dec 2020). — [QUID](https://documentation.marine.copernicus.eu/QUID/CMEMS-INS-QUID-013-045.pdf); [SEANOE](https://www.seanoe.org/data/00591/70345/)
- IMOS/AODN (Australia): wave buoy data are CC BY 4.0; delayed-mode data go to the National Wave Archive as NetCDF via THREDDS; the Dec 2025 QA/QC manual says raw displacements are sampled at 2.5 Hz and that the facility archives raw displacement data. — [IMOS record](https://researchdata.edu.au/imos-new-technology-wave-buoys/1714710); [IMOS QA/QC doc](https://content.aodn.org.au/Documents/IMOS/Facilities/coastal_wave_buoys/IMOS_Coastal_Wave_Buoys_QAQC_Documentation_v1.0.pdf)
- North Sea fixed sensors: the southern North Sea radar data in Teutsch et al. were provided by Shell (credited to Graham Feld and Shell); Ekofisk data are described as collected by Phillips Petroleum and "kindly provided"; Draupner is known from two published 20-minute series; FINO1 data are downloadable from the BSH portal (confirmed in a wind-validation paper; wave elevation series not confirmed). — [JPO 2022](https://journals.ametsoc.org/view/journals/phoc/53/1/JPO-D-22-0059.1.xml); [NORA3 paper](https://wes.copernicus.org/articles/6/1501/2021/); [Draupner](https://en.wikipedia.org/wiki/Draupner_wave)
- Christou and Ewans (JPO 44(9), 2014): 122 million quality-controlled waves from fixed platform sensors (22 North Sea sites, 5 Gulf of Mexico, 5 South China Sea, 1 NW Australia), 3649 rogue waves. — [VLIZ copy](https://www.vliz.be/imisdocs/publications/ocrd/381785.pdf)
- A fixed-tower lead: "45 years of directional wave recorded data at the Acqua Alta oceanographic tower" (Adriatic; data paper, 2025). — [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11802720/)
- USACE QCC Measurement Archive: quality-controlled NDBC-derived archive (data paper 2022). — [Sci Data](https://www.nature.com/articles/s41597-022-01344-z)

### Inferences
- Sofar's open archive and the Copernicus in-situ product give sea-state parameters and spectra but (per their descriptions) not wave-by-wave elevation, so they cannot supply rogue-wave labels; they could only test whether the model's input distribution shifts in the open ocean.
- IMOS (2.5 Hz Spotters and Waveriders, Southern Ocean swell climate, CC BY) is the most promising non-US buoy source if the raw displacement files are actually downloadable, but it is still a surface-following buoy and so does not resolve the sensor-bias question.
- A fixed-sensor test almost certainly requires a data request to an operator or institute (Shell, ConocoPhillips, Equinor, MET Norway, BSH, CNR-ISMAR), which is feasible for an unaffiliated person only if the holder agrees.

### Gaps
- Not verified: whether IMOS raw displacement is public; whether Acqua Alta includes elevation time series; FINO1 wave time-series terms; any open Ekofisk/Gorm archive.
- No usable search results on Canadian MEDS, Japanese NOWPHAS, or satellite products as rogue-wave validation sources. Satellites measure Hs and spectra (SAR), not single-point wave-by-wave elevation, so their role would be limited to inputs; this is background knowledge, not sourced here.
- Sizes of the non-CDIP archives in wave counts were not found.

## 3. Buoy versus fixed-sensor disagreement and likely transfer

### Takeaway
Published comparisons consistently find fewer rogue waves and lower crests from surface-following buoys than from radar/laser sensors, with the cause (sampling rate, Lagrangian motion, footprint) not settled. A model fitted to buoy wave heights is therefore unlikely to transfer in absolute level without recalibration; whether its parameter dependencies transfer is untested.

### Cited Findings
- Teutsch et al. (NHESS 20, 2665, 2020), southern North Sea, about 329 million waves in about 790,000 half-hour samples, 2011-2016: radar stations recorded more rogue waves than buoys; radar counts exceeded Forristall expectations at every station while buoys followed Forristall more closely or fell below second-order expectation; buoys were Waverider MkIII at 1.28 Hz, radars 2 to 4 Hz; subsampling radar data to buoy-like rates was tested; the authors say the gap may be due to measurement area or the systems themselves. — [NHESS](https://nhess.copernicus.org/articles/20/2665/2020/); [EGU poster](https://www.egu.eu/media/awards/ospp-award/2019/ina_teutsch.pdf)
- A later brief communication attributes the gap to technique: a high-resolution radar "is able to identify more rogue waves than a surface-following wave buoy with a lower sampling frequency". — [NHESS 2024](https://nhess.copernicus.org/articles/24/2065/2024/nhess-24-2065-2024.pdf)
- Ekofisk sensor comparisons: many Waverider records show negative skewness in higher sea states and lower crest heights than laser and radar; laser Hs is about 0.1 m or about 2.4% higher than Waverider. — [ResearchGate summary](https://www.researchgate.net/publication/363338759_On_the_Functionality_of_Radar_and_Laser_Ocean_Wave_Sensors); [Ocean Dynamics 2021](https://link.springer.com/article/10.1007/s10236-021-01468-7)
- Mixed evidence: a JTECH 2023 three-sensor comparison found radar gave the lowest Hs while laser and Waverider agreed well for wind sea and swell, with the laser higher at high frequency; a search summary of it says a Waverider "should not underestimate crest heights in some conditions". — [JTECH 2023](https://journals.ametsoc.org/view/journals/atot/40/5/JTECH-D-22-0108.1.xml)
- A Datawell technical note cites a report that the Waverider underestimates spectral energy above 0.3 Hz in high sea states. — [Datawell note](https://datawell.nl/wp-content/uploads/2022/08/datawell_publication_wr_oceanengineering-vol20.pdf)
- Baschek and Imai (2011), 16 Waverider buoys off California/Oregon: buoy statistics fall below the Gaussian curve; a citing paper adds that buoys tend to underestimate amplitude statistics (Allender et al. 1989). — [Baschek & Imai](https://tos.org/oceanography/assets/docs/24-2_baschek.pdf); [Cattrell et al. accepted MS](https://eprints.soton.ac.uk/422951/1/Paper2018JC013958_acceptall_final.pdf)
- Christou and Ewans 2014 (fixed sensors) concluded rogue waves "are not governed by sea state parameters" and are rare members of the normal population produced by dispersive focusing. — [VLIZ copy](https://www.vliz.be/imisdocs/publications/ocrd/381785.pdf)

### Inferences
- Rough arithmetic, not like-for-like: the local FOWD-CDIP files give about 0.8% of 100-wave blocks with H/Hs > 2, i.e. at least about 8e-5 per wave; Rayleigh gives exp(-8) = 3.35e-4; Christou and Ewans' 3649 in 122 million is about 3e-5. Definitions, QC and sea-state mixes differ, so these numbers show only that absolute rates are strongly dataset-dependent.
- The PNAS "Rayleigh is an upper bound" result could partly reflect buoy under-measurement rather than physics; the authors themselves list sensor bias as a risk. Testing this needs fixed-sensor data, which is the least accessible resource.
- The fixed-sensor literature and the buoy literature disagree on whether sea-state parameters matter at all (Christou and Ewans: no; Häfner: yes, via r and spread), which is itself an acknowledged open question that an outsider cannot resolve without the industry data.

### Gaps
- No source quantifies the Lagrangian or 1.28 Hz bias on rogue-wave rates against a co-located laser; the specific Magnusson/Donelan and Casas-Prat and Holthuijsen results were not retrieved.
- Exact radar versus buoy rogue counts in Teutsch et al. were truncated in snippets (1234 radar rogue samples quoted; buoy figure unclear).

## 4. Compute needed for meaningful contributions

### Takeaway
Evaluating the published closed-form model, calibration analysis and extreme-value statistics fit easily on 4 cores and 15 GB; the paper's full neural-network search as packaged does not (64 GB RAM, about 2 days on 16 cores), though a column-pruned re-implementation probably would.

### Cited Findings
- The repository README states that running all 24 experiments "will take a long time (about 2 days on a 16-core machine) and requires at least 64 GB of RAM", that reproducing figures "requires at least 32 GB of RAM", and that the data download is about 6.2 GB. — [rogue-wave-discovery README](https://github.com/dionhaefner/rogue-wave-discovery) (local copy read)
- The network is small: MLP layers (32, 16, 8), 50 epochs plus 50 SWAG epochs, batch size 1024, learning rate 1e-4, CPU JAX (`jax[cpu]==0.4.13`). — [repo constants/requirements](https://github.com/dionhaefner/rogue-wave-discovery) (local)
- Symbolic regression in the paper's notebook uses PySR 0.16.0 on a 100,000-row subsample with `niterations=40`, `populations=48`, `procs=48`, `maxsize=32`. — [repo notebook](https://github.com/dionhaefner/rogue-wave-discovery) (local)
- Reading five columns from all 157 files (15.0 million rows) with pyarrow took 1.5 s on the 4-core machine. — [repo data](https://github.com/dionhaefner/rogue-wave-discovery) (local measurement)
- PySR supports mini-batching for large datasets and a `procs` setting equal to cores used. — [PySR options](https://huggingface.co/spaces/MilesCranmer/PySR/blob/main/docs/options.md)
- External reference point: on a 1e5-point dataset one benchmark ran PySR for 10 h on 24 cores and found the same result after 5 h. — [Exhaustive Symbolic Regression](https://arxiv.org/pdf/2211.11461)

### Inferences
- The 64 GB figure probably comes from loading all 106 columns into pandas; ten float64 columns over 15 million rows is about 1.2 GB, so a pruned loader should fit in 15 GB. This is an estimate, not tested.
- Scaling the README figure linearly (32 core-days for 24 experiments) gives roughly 8 hours per experiment on 4 cores, so a handful of re-fits is feasible but a wide architecture search is not.
- PySR at 4 procs instead of 48 means roughly an order of magnitude less search per wall-clock hour; rediscovering or modestly varying the published formula is plausible over days, a broader search is not.
- Processing raw CDIP displacement through FOWD for an out-of-time test is I/O-bound and parallel by deployment; feasibility on 4 cores is unknown.

### Gaps
- No source gives FOWD processing time or memory. The paper's own compute (hardware, GPU or not) was not found in searches.

## 5. Role for formal verification / Lean 4

### Takeaway
Formal extreme-value theory is an active but tiny niche (one course project), ocean-wave statistics have not been formalised anywhere I could find, and Mathlib lacks most named continuous distributions; this is a real library gap but a mismatch with the actual bottleneck in rogue-wave prediction, which is empirical (sensor bias, data coverage), not deductive.

### Cited Findings
- Kalle Kytölä's ExtremeValueProject (Aalto course MS-EV0029, Feb to June 2025) aims to formalise the Fisher-Tippett-Gnedenko theorem in Lean; its blueprint says it "is still work in progress" and defines Gumbel, Weibull and Fréchet types via CDFs. — [blueprint](https://kkytola.github.io/ExtremeValueProject/blueprint); [Reservoir](https://reservoir.lean-lang.org/@kkytola/ExtremeValueProject)
- Mathlib has Gaussian (`Mathlib.Probability.Distributions.Gaussian.Real`), exponential, gamma, Poisson (PMF) and uniform; a Feb 2026 Zulip thread notes missing properties for the exponential; a maintainer remarked that only the uniform and Gaussian are used elsewhere in Mathlib. No Rayleigh, Weibull, Gumbel or Pareto file was found. — [Zulip: exponential](https://leanprover-community.github.io/archive/stream/287929-mathlib4/topic/Exponential.20Distribution.html); [Zulip: distributions](https://leanprover-community.github.io/archive/stream/217875-Is-there-code-for-X%3F/topic/Probability.20Distributions.html)
- Mathlib probability otherwise covers Markov kernels, sub-Gaussian variables, Markov/Chebyshev inequalities and Doob's martingale convergence. — [Degenne preprint](https://arxiv.org/pdf/2510.04070); [Doob formalisation](https://arxiv.org/pdf/2212.05578)
- Fluid mechanics in Lean exists at the PDE-theory level: a Lean 4 formalisation of the Galerkin construction for 2D Navier-Stokes (arXiv, Sept 2026), community Navier-Stokes projects, and a Nov 2026 arXiv paper using Lean as the deduction layer for certified computation of nonlinear water waves. — [2D Galerkin](https://arxiv.org/pdf/2609.33033); [water waves](https://arxiv.org/pdf/2611.07395)
- Related formal-methods work: Coq proofs of a wave-equation solver (Boldo et al.); formal methods for storm-surge finite-element code (Baugh and Altuntas, cited in a climate-modelling survey). — [Boldo et al.](https://arxiv.org/pdf/1001.4898); [survey](https://arxiv.org/pdf/2510.13425)
- Physics formalisation is growing (PhysLean; Lean4PHYS at ICLR 2026), but practitioners stress that Lean checks logic, not physics: physical "lore becomes named hypotheses", and proof assistants "don't say anything about the physics context of results". — [Lean4PHYS](https://proceedings.iclr.cc/paper_files/paper/2026/file/2821de2dc3501f55732e59d053e9e0cb-Paper-Conference.pdf); [Seiberg-Witten paper](https://arxiv.org/pdf/2607.06379); [PhysLean author blog](https://josephtoobysmith.com/blog.html)
- A 2026 position paper says formalisation is not yet practical for people without Lean expertise. — [arXiv 2610.14059](https://arxiv.org/pdf/2610.14059)
- No search result showed formalised ocean-wave statistics, Rayleigh/Tayfun/Forristall distributions, or peaks-over-threshold theory in Lean, Coq or Isabelle. — [statistical theory survey](https://arxiv.org/pdf/2609.04481)

### Inferences
- Feasible small Lean targets: define the Rayleigh distribution and prove the exceedance P(H > 2Hs) = exp(-8) under the narrow-band Gaussian hypotheses; the block-maximum identity 1 - (1 - p)^n behind the 100-wave aggregation; monotonicity and limit properties (e.g. the spread-to-zero divergence) of the published formula. These would be new to Mathlib-adjacent libraries but none bears on whether the model is right about the ocean.
- A formalised result would have no natural audience in oceanography; its audience is the Lean community. It should be judged as a mathematics-library contribution with a maritime theme, not as progress on rogue-wave prediction.

### Gaps
- I could not inspect the Mathlib source tree or the ExtremeValueProject dependency graph, so current completion status and the exact distribution list are unverified.

## 6. Counterarguments: why this may not be a good use of time

### Takeaway
Several large field studies conclude that sea-state parameters carry little predictive signal for individual rogue waves, the best-supported indicator (crest-trough correlation) is already published and being followed up by the original group, the public data are single-network buoys with known bias, and I found no evidence that any rogue-wave probability product changes operational decisions.

### Cited Findings
- Weak signal: Cattrell et al. (JGR Oceans 123, 2018) found that in bulk "only weak relationships can be seen" and that predictors will need to be location-specific. — [NOC record](https://nora.nerc.ac.uk/id/eprint/520662/)
- Christou and Ewans 2014 (122 million waves, fixed sensors): rogue waves "are not governed by sea state parameters". — [VLIZ copy](https://www.vliz.be/imisdocs/publications/ocrd/381785.pdf)
- Teutsch et al. 2020: predictors "can probably not be derived for an entire dataset"; about 55% probability of correctly detecting a rogue wave from background-field knowledge (search-summary wording). — [ResearchGate copy](https://www.researchgate.net/publication/339954829_A_Statistical_Analysis_of_Rogue_Waves_in_the_Southern_North_Sea)
- The signal that exists is modest in absolute terms: probability rises by about a factor of 10 from low to high r (as summarised for the Gemmrich group's work), and the model's worst case equals the Rayleigh value 3.3e-4 per wave. — [Weather and Forecasting 2023](https://journals.ametsoc.org/view/journals/wefo/38/11/WAF-D-23-0074.1.xml); [PNAS](https://www.pnas.org/doi/full/10.1073/pnas.2306275120)
- Base-rate problem: Breunung and Balachandran (Sci Rep, July 2024; 14 million 30-minute samples from 172 buoys) note an always-negative classifier scores 99.86% accuracy; their network catches about 75% of rogue waves 1 minute ahead and 73% at 5 minutes; generalisation was tested on two held-out buoys; they say they address prediction, not causality. — [Sci Rep 2024](https://www.nature.com/articles/s41598-024-66315-3)
- The niche is occupied by specialists: Cicon, Gemmrich et al. already built a probabilistic forecast from WAVEWATCH III using a semilog fit on r (2023), noting r is forecast with only "moderate accuracy"; Azevedo et al. (2026) and an Ocean Dynamics 2025 study of 440 eyewitness incidents extend the r line of work. — [WAF 2023](https://journals.ametsoc.org/view/journals/wefo/38/11/WAF-D-23-0074.1.xml); [Ocean Science 2026](https://os.copernicus.org/articles/22/367/2026/); [Ocean Dynamics 2025](https://link.springer.com/article/10.1007/s10236-025-01697-0)
- Models cannot yet supply the inputs well: ERA5 and ECMWF hindcasts "systematically underestimate extreme wave events due to spectral smoothing and spatial averaging". — [Ocean Science 2026](https://os.copernicus.org/articles/22/367/2026/)
- Operational context: ECMWF already outputs freak-wave parameters (BFI, kurtosis-based) but its documentation calls the kurtosis expression "far too involved" for operations; I found no data on shipping uptake or any documented case of a vessel avoiding a rogue wave due to a forecast. — [ECMWF wave model documentation](https://www.ecmwf.int/sites/default/files/elibrary/2016/79992-ifs-documentation-cy43r1-part-vii-ecmwf-wave-model_1.pdf); [ECMWF fact sheet](https://www.ecmwf.int/en/about/media-centre/focus/2020/fact-sheet-ocean-wave-forecasting)
- Even a FOWD-based paper by a research group drew substantive referee objections on indicator choice and comparability. — [OS peer review](https://os.copernicus.org/articles/22/367/2026/os-22-367-2026-discussion.html)
- Outsider publishing: a secondary source says IOP and Cambridge University Press accept "Independent Researcher" as an affiliation; fee waivers for unaffiliated authors are partial or case-by-case at the journals found. — [blog source](https://www.alivelearn.net/?p=5074); [ICES JMS](https://academic.oup.com/icesjms/pages/open-access)
- Counter-evidence in favour: the code and data are fully public and a university news release says authorities "can easily begin calculating" the probability; secondary analyses of FOWD do get published (Azevedo et al.). — [EurekAlert](https://www.eurekalert.org/news-releases/1008299); [Ocean Science 2026](https://os.copernicus.org/articles/22/367/2026/)

### Inferences
- The strongest reasons against: (1) the remaining predictable variance beyond r, spread, steepness and depth is plausibly small and may sit below the sensor-bias floor of the only open data; (2) a marginal skill gain in a statistical model has no identified user, because the limiting factor for forecasts is model-predicted r and spectra, not the mapping from parameters to probability; (3) the decisive experiment (fixed-sensor or open-ocean validation) needs data the person cannot download; (4) domain reviewers will weight measurement physics the non-specialist lacks.
- The most defensible rogue-wave contribution is therefore a replication and audit (calibration by region, crest-height variant, QC sensitivity, out-of-time CDIP test), which has honest scientific value as a negative or confirming result but low odds of uptake.

### Gaps
- No systematic evidence was found on acceptance rates or citation uptake for unaffiliated authors in physical oceanography; the outsider-disadvantage claim is plausible but unquantified here.
- Copernicus/Ocean Science fee levels and waiver rules for unaffiliated authors were not retrieved.
- The Dudley et al. 2019 review's position on predictability could not be checked.

## 7. Alternative open maritime problems for the same resources

### Takeaway
Container stowage planning (public benchmark suites with objective scores), reproducible ship weather routing, and evaluation of search-and-rescue drift-to-search pipelines look more tractable for a 4-core, no-funding effort than rogue-wave prediction; storm-surge and sneaker-wave work have clearer public-safety users but are respectively crowded and data-poor. Parametric roll, AIS collision risk and fuel prediction are poor fits.

### Cited Findings
Container stowage / berth planning:
- A 2024 EJOR literature survey publishes 162 benchmark instances (DTU repository, with a Julia reader on GitHub) and a research agenda; Larsen and Pacino (2021) give 3 vessel profiles and 27 instances; a 2025 paper introduces RCSPPSuite, described as the largest suite with industrial-sized instances. — [survey](https://arxiv.org/pdf/2307.07573); [benchmark repo](https://github.com/maritimelab/StowagePlanning_MasterPlanning_Benchmark); [RCSPPSuite](https://www.sciencedirect.com/science/article/abs/pii/S1366554525002650)

Weather routing / voyage optimisation:
- VISIR-2 (GMD 17, 4355, 2024) is an open Python routing model using waves, currents and wind, archived on Zenodo; its authors state the potential of weather routing "has not been reproducibly assessed, especially for short sea shipping". SIMROUTE (Copernicus Marine products) is described as a candidate for inter-comparison exercises, but no published comparison was found. — [GMD](https://gmd.copernicus.org/articles/17/4355/2024/); [SIMROUTE](https://marine.copernicus.eu/services/use-cases/simroute-comprehensive-ship-weather-routing-system-using-copernicus-marine)
- A multi-basin open benchmark (Aredah and Rakha 2026) and an Ocean Engineering 2026 result that routing on a failed deterministic forecast was worse than weather-blind routing are cited only second-hand in a GitHub issue; unverified. — [GitHub issue](https://github.com/Luminous-Dynamics/symthaea/issues/11587)

Search and rescue:
- A JMSE review (vol. 14, issue 20; 2026 by volume) on linking drift prediction to search planning argues that better drift accuracy does not necessarily produce better search outcomes and calls for calibration and sharpness evaluation of uncertainty-aware drift models under regional/seasonal shift. — [JMSE review](https://www.mdpi.com/2077-1312/14/20/1886)
- OpenDrift is used openly for leeway simulation (Frontiers 2026 Korean case study against the IAMSAR worksheet); leeway coefficients come from Allen and Plourde's 1999 USCG compilation (95 target types, 25 field studies, 63 classes); a methods paper recommends that raw leeway experiment data be published, implying much is not. — [Frontiers 2026](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2026.1852504/full); [Allen & Plourde](https://apps.dtic.mil/sti/tr/pdf/ADA366414.pdf); [leeway field method](https://arxiv.org/pdf/1111.0750)
- Search-path RL papers (2025) are simulation studies on OpenDrift probability maps with no shared benchmark found. — [maskable PPO](https://www.kjrs.org/journal/view.html?doi=10.7780%2Fkjrs.2025.41.6.17); [improved DQN](https://www.tandfonline.com/doi/full/10.1080/17538947.2025.2562555)
- Bellman's "lost at sea" problem: Isbell solved the minimax case (1957; full proof published by Joris 1980); Gluss (1961) gave approximations for the expected-distance criterion; the general lost-in-a-forest problem is listed as unsolved, with 2024-2026 preprints claiming partial or general solutions not yet settled. No formal verification of any of these was found. — [Gluss 1961](https://onlinelibrary.wiley.com/doi/abs/10.1002/nav.3800080108); [Wikipedia](https://en.wikipedia.com/wiki/Bellman%27s_lost-in-a-forest_problem); [arXiv 2412.10686](https://arxiv.org/pdf/2412.10686); [arXiv 2601.21867](https://arxiv.org/pdf/2601.21867)

Sneaker waves / coastal hazard:
- A 2018 AGU abstract counts 27 fatalities directly attributed to sneaker waves in Northern California and Oregon since 2005 and finds most occurred with peak period in the top quintile; an OSU-linked report cites about two drowning deaths a year across Oregon, Washington and northern California; OSU researchers name a localised predictive index as the goal; NWS Portland issues sneaker-wave statements. — [AGU abstract](https://ui.adsabs.harvard.edu/abs/2018AGUOSCD14B0033G/abstract); [OSU news](https://news.oregonstate.edu/news/far-storms-fuel-sneaker-waves-along-pacific-northwest-coast-new-research-suggests); [monitoring report](https://www.fondriest.com/news/monitoring-for-runup-signals-to-reduce-sneaker-wave-risk.htm)

Storm surge / tide:
- GESLA-3 hourly tide-gauge records with ERA5 predictors are the standard open setup; a 2020 global study covered 610 gauges, a 2025 European study only 9; known pitfalls are tidal-decomposition error as large as the signal and 30 to 50% underestimation of extremes with daily forcing. — [JGR ML 2025](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025JH000650); [EGUsphere 2025](https://egusphere.copernicus.org/preprints/2025/egusphere-2025-196/egusphere-2025-196.pdf); [ERL 2020](https://iopscience.iop.org/article/10.1088/1748-9326/ab89d6)

AIS collision / grounding risk:
- A comprehensive ML-on-AIS review exists; collision-risk labels are typically computed from kinematic formulas on simulated encounters rather than observed near-misses; grounding has far fewer papers than collision; a large Finnish AIS dataset was released in 2026. — [review](https://www.researchgate.net/publication/377750904_Harnessing_the_power_of_Machine_learning_for_AIS_Data-Driven_maritime_Research_A_comprehensive_review); [AIS-LLM](https://arxiv.org/pdf/2508.07668); [Finnish AIS dataset](https://arxiv.org/pdf/2609.12938)

Parametric roll; fuel and emissions:
- Parametric-roll ML studies rely on numerical simulation or tank tests (e.g. C11 and KCS hulls) with no public benchmark found. — [Ocean Engineering 2023](https://www.sciencedirect.com/science/article/abs/pii/S0029801823005218)
- EU MRV (THETIS) data are public but annual per-ship aggregates; a daily dataset (FuelCast) is CC BY-NC-ND 4.0; no standard benchmark split was found and reported fits on annual data (R² about 0.99) are easy. — [MRV benchmarking paper](https://www.mdpi.com/2071-1050/18/18/9545); [physics-informed ML paper](https://www.mdpi.com/2077-1312/14/19/1774)

### Inferences
Assessment per candidate (what is open / data / users / crowding), my judgement from the findings above:
- **Container stowage planning**: open in the sense that benchmark suites are recent and best-known solutions can be beaten; data are public instances; users are liner-shipping planners and OR researchers; moderately crowded but scored objectively, so an outsider's result is verifiable without credentials. Runs on 4 cores. Best fit for "tractable". Weakest on direct safety value.
- **Weather routing reproducibility**: open question is an independent, reproducible estimate of savings and of risk when forecasts fail; open code (VISIR-2) plus ERA5/Copernicus data; users are operators and decarbonisation analysts; academically crowded and commercially dominated, but reproducible comparisons are explicitly scarce.
- **SAR drift-to-search evaluation**: open question is whether drift-model improvements improve search outcomes, and calibration of drift uncertainty; open tools exist but open leeway ground truth is thin; users are coast guards, who are hard for an outsider to reach. Highest safety value, uncertain data.
- **Lean-suited maritime mathematics**: formalising Isbell's lost-at-sea minimax path is the one place the Lean setup clearly matches a maritime-named problem; it is a well-defined, apparently unformalised result, but its practical maritime value is essentially nil.
- **Storm surge**: fully open data and small-compute templates, real users (coastal agencies), but crowded with well-resourced groups; a small effort would most plausibly add a careful benchmark or negative result.
- **Sneaker waves**: clear public-safety user and few researchers, but I found no open event catalogue or runup dataset; likely blocked on data.
- **AIS collision risk, parametric roll, fuel prediction**: respectively crowded with modelled labels, lacking open data, and limited by coarse public data; weaker targets.

### Gaps
- Not verified: MarineCadastre AIS terms and size; Global Drifter Program suitability for leeway validation; USCG SAROPS evaluation literature; berth-allocation benchmark sets; any sneaker-wave open dataset; the Aredah and Rakha benchmark.
- Crowding is judged qualitatively from search results; no paper counts were obtained.
