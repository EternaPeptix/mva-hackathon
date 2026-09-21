> Package copy: local links/navigation are adapted for this finite review package. Labels identify underlying files that are not included. Scientific text is retained; approved status/source metadata changes are recorded separately; original and packaged hashes are recorded in PROVENANCE.

> **Retained dated method detail.** This appendix preserves the original analysis or proposal at its stated date. The current [Track 1 report](../jeweled_track1_report_submission.md), [Track 2 report](../jeweled_track2_report_submission.md) and [methods form](../jeweled_methods_form_answers.md) govern interpretation; see the [revision and evidence scope](../PROVENANCE.md#evidence-scope). Historical delivery or ranking statements are not current status.

# MVA pilot protocol and feasibility packet

20 September 2026 UTC. **Prospective design for collaborator review; no laboratory experiment
has been performed.** Material access, a laboratory partner, exposure ranges and confirmatory
sample size are unconfirmed. The quantities below are planning scenarios, not observed effects,
selected patient doses or a promise of completion within a month.

This packet makes the [assay decision plan](assay-decision-plan.md) operational.
The question is whether a verified BUB1B defect can be improved without merely preventing
division, losing damaged cells or selecting a different abnormal clone. Arimoclomol remains
an upstream mechanistic probe and acetylcysteine a downstream comparator. The
[Track 2 report — archival document not included](../PROVENANCE.md#archival-material) retains their exposure and evidence limits.

## 1. Model qualification before drug screening

The preferred case-relevant comparison is a verified patient-derived model with a precisely
corrected derivative and appropriate editing controls. Cell availability is not established.
An engineered model can test an allele mechanism, but cannot establish the child's phase,
tissue response or diagnosis. Do not replace a missing model with an assumed genotype.

While phase is unknown, the proposed engineered comparison has six model classes:

| Class | Purpose |
|---|---|
| Sham-edited wild type | Background and editing-process control |
| Leu737Ter / wild type | Effect of the nonsense allele alone |
| Asn1002Lys / wild type | Effect of the missense allele alone |
| The two alleles in trans | Test the proposed two-copy mechanism |
| Both alleles in cis, with a wild-type partner | Test the unresolved alternative configuration |
| Precisely corrected derivative of the affected engineered model | Test genetic reversibility in a related background |

These classes are proposed, not available lines. Verify endogenous genotype, local phase,
cell identity, culture quality, passage history and baseline copy-number state before assigning
labels. Preserve parental and sham controls, and document each independently derived editing
family. Do not disable p53 or select only fast-growing clones to make the primary model convenient.
Record substrate, cell density, oxygen and cell-cycle distribution consistently across arms.

First establish RNA/protein effects, a reproducible functional defect and interpretable genetic
rescue. Quantify protein in a matched cell-cycle context rather than interpreting a changed bulk
population composition as restored protein regulation. A cDNA expression system may examine
protein behavior, but cannot establish endogenous splicing or nonsense-mediated RNA decay.
If model construction fails or the phenotype is not reproducible, report that limitation;
it does not disprove the patient's allele or justify replacing the causal model silently.

## 2. Observe the first division of defined starting cells

Select fields before treatment under a masked, prespecified sampling rule. Assign a founder ID
to each eligible cell present at time zero. Record baseline cell-cycle information when available.
Follow that founder's **first** mitosis and its descendants. A daughter's later division is a
secondary lineage event, not another independent founder or primary observation.

At the prespecified acute endpoint, every starting founder must have exactly one first-event
status. The proposed first feasibility window is 48 hours; the pilot must establish whether it
captures an interpretable cell cycle without excessive censoring in the chosen model.

| Status | Definition and treatment in analysis |
|---|---|
| A: scorable anaphase in the first mitosis | Visible-error status is scored using the frozen image criteria |
| U: the first mitosis reached anaphase but error status was unscorable | Retain the observation; never recode as no error |
| B: tracked to the endpoint without mitotic entry | Keep as a competing outcome |
| S: adjudicated exit from the first mitosis without an observed anaphase | Record slippage separately from successful division |
| D: death before a first-mitosis outcome could be assigned | Record whether interphase or mitotic, without dropping the founder |
| O: still in the first mitosis at the window boundary | Administrative endpoint, not a confirmed permanent arrest |
| L: tracking lost before classification | Record time/reason; do not infer death, arrest or normal division |

A, U or S fixes the founder's first-mitosis category. Later death, a second division or
tracking loss does not overwrite it; those remain secondary outcomes. D and L apply only
before that first category is determined. B and O require follow-up through the endpoint.

The census must satisfy `N0 = A + U + B + S + D + O + L`. Within A, let E be the number
with at least one prespecified visible segregation error. Keep error types and cytokinesis
outcomes separately; a cell with several error types counts once in E. Imaging criteria must
distinguish lagging material, bridges, multipolar segregation and uncertain events. Reviewers
remain masked to treatment and genotype until scoring and adjudication are frozen.

The proposed **linked primary readouts** are:

- `E / A`: visible-error frequency among scorable first anaphases; undefined when A is zero.
- `(A − E) / N0`: observed first-anaphase yield without a scored visible error per starting
  founder; this is not proof of a normal karyotype or successful daughter-cell function.

Also report A/N0, every competing-status proportion, mitotic timing, cytokinesis failure,
daughter death and tracking completeness. A drug that changes entry or survival changes the
population represented by E/A. It cannot advance on a smaller E/A alone. If U is nonzero,
the error fraction among observed first anaphases lies between `E/(A+U)` and `(E+U)/(A+U)`
without assuming how those unscorable events behaved. This does not solve missingness in L;
retain separate censoring and sensitivity analysis.

For illustration only, 100 founders with 80 scorable anaphases and 16 errors yield 20% errors
and 64% observed yield without a scored error. A constructed cytostatic scenario with 20
anaphases and 2 errors has a lower 10% error rate but only 18% yield. These are arithmetic
counterexamples, not simulated or observed results for either candidate drug.

## 3. Separate acute correction from later selection

Use spontaneous divisions for the primary comparison. A checkpoint or oxidative-stress
challenge can serve as a separately labeled assay-control arm; it cannot replace the baseline
genotype-dependent phenotype. In all arms, retain live-cell counts and the full founder census.

Link daughter fate and clone abundance back to the first scored division. Proposed sampling
includes baseline, the end of acute exposure and a further 72 hours after washout, with the
actual calendar window frozen after assay calibration. Sample karyotypes under a condition-blind
rule and record the number selected, successfully assayed and unresolved. Report absolute
abnormal-cell/lineage counts as well as their fraction. No measured chromosome state should be
assigned to unassayed cells solely from appearance or a relative growth signal. A limited FISH
panel cannot establish genome-wide normality, and bulk copy number can miss small subclones.

Longer adaptation surveillance is a later requirement if a signal advances; this short pilot
cannot exclude long-term selection. In human cell-line experiments, selected chromosome changes
under sustained checkpoint inhibition altered mitotic timing and reduced instability. That is
a reason to examine lineage and karyotype, not evidence that our drugs do this in the child.
[Primary adaptation study](https://pubmed.ncbi.nlm.nih.gov/36859339/).

Other primary models link missegregation to arrest/clearance and nuclear-state changes.
Their relevance here is the need to distinguish restored function from escape from an arrest
response; they do not establish BUB1B-specific drug efficacy.
[Arrest and clearance](https://pubmed.ncbi.nlm.nih.gov/28633018/),
[nuclear checkpoint](https://pubmed.ncbi.nlm.nih.gov/39779939/).

## 4. Exposure and measurement qualification

The three nonzero concentrations in the capacity example below are **unfilled planning slots**.
Before drug work, the pharmacologist and assay laboratory must complete a concentration/time
table specifying formulation, analyte, measured medium concentration, stability, free/total
interpretation, sampling times, calibration range and the source supporting relevance. Preserve
the completed table before unmasking response data.

For arimoclomol, the previously reviewed high-concentration NPC experiments do not set a
BUB1B assay target. For NAC, total NAC measured after reduction is not interchangeable with
free reduced NAC, cysteine or glutathione. The [arimoclomol](mechanistic-bridge.md) and
[NAC](nac-exposure.md) reviews remain controlling evidence boundaries.
No case-specific achievable exposure range has been demonstrated.

Measure exposure/recovery separately from biological response. Keep intentionally high
mechanistic controls separate from the translational comparison. Include medium/pH and
cell-free fluorescence-interference controls, and an orthogonal functional readout for the
stress branch. ROS dye intensity alone cannot establish mitochondrial rescue. Vehicle and
single-drug conditions precede combinations. If exposure cannot be justified or measured,
report exploratory engagement rather than practical repurposing efficacy.

## 5. Pilot allocation and inferential units

The following is an **illustrative engineered-line resource plan**, contingent on the required
models. It is not a power calculation. Three editing founder families, each contributing the
specified model classes, are observed across three independent culture preparations, with two
technical wells per condition. Culture preparation, editing family, clone and well identifiers
must remain distinct. Thousands of cells or frames do not turn one engineered background or one
patient into a patient cohort.

| Stage | Proposed factors | Wells per readout |
|---|---|---:|
| Baseline/model qualification | 6 model classes × 3 founder families × 3 preparations × 2 wells | 108 |
| One drug after model qualification | Affected model, corrected derivative and sham wild type × 3 families × 3 preparations × (vehicle + 3 exposure slots) × 2 wells | 216 |
| Two drugs run separately | Two complete single-drug stages | 432 |

Additional cell-free controls, destructive RNA/protein assays, karyotyping, failed acquisitions,
high-exposure controls and confirmatory experiments are not included in these counts. If only
patient material is available, design a separate patient/corrected-control allocation rather
than relabeling patient cultures as independent editing families. Reducing scope because of
availability must be explicit; omitted configurations remain untested.

Randomize conditions within appropriate preparation/founder blocks and balance plate position,
acquisition order and time. Keep a masked assignment key separate from the scoring file. The
pilot estimates variation and missingness; it does not make confirmatory efficacy claims.
Analyze matched treatment contrasts with the actual hierarchy and repeated lineages retained.
With only three founder families, uncertainty about between-family variation will be substantial.

Before fresh confirmation, freeze the primary contrasts, minimal meaningful functional change,
permitted division/survival tradeoffs, exposure bounds, censoring/exclusion rules, multiplicity
policy and sample-size rationale using pilot estimates. A failed significance test is not an
equivalence or safety result. There is no invented numerical efficacy threshold in this packet.

## 6. Imaging, storage and calendar feasibility

The planning calculator assumes 48 hours at five-minute intervals, including both endpoints
(577 frames), four fields per well, three z planes, two channels, 1024 × 1024 pixels and
two stored bytes per pixel. These are **unvalidated acquisition assumptions**. Pilot microscopy
must demonstrate adequate focus, temporal resolution, tracking and acceptable phototoxicity.
Do not reduce acquisition quality solely to make the capacity estimate smaller.

Raw pixel bytes are `wells × 4 × 3 × 2 × 577 × 1024 × 1024 × 2`.
The table uses decimal TB, assumes no compression gain, and excludes metadata, analysis outputs,
temporary files and storage headroom.

| Stage | One raw copy | Two raw copies | Acute instrument days at 24 simultaneous wells | At 72 simultaneous wells |
|---|---:|---:|---:|---:|
| Baseline | 3.136 TB | 6.273 TB | At least 12 | At least 6 |
| One drug | 6.273 TB | 12.546 TB | At least 18 | At least 6 |
| Baseline plus two separate drugs | 15.682 TB | 31.365 TB | At least 48 | At least 18 |

These lower bounds preserve three separate culture preparations. They exclude model construction,
calibration, washout imaging, setup and failed/repeated runs; they are not a promised schedule.
At 24 wells, one acquisition timepoint requires 576 image planes, allowing only 0.521 seconds
per plane including stage movement, readout, focus and other overhead within five minutes.
At 72 wells the corresponding limit is 0.174 seconds across 1,728 planes. Actual throughput
must be measured. A single full 24-well session produces about 697 GB of raw pixels under
these assumptions, so the current genome-analysis machine is not assumed to provide the
required image archive or instrument capacity.

Secure dedicated acquisition storage and a separately verified retained copy before starting.
Confirm bandwidth and data-integrity checks. Obtain written laboratory/core estimates for
model creation/QC, imaging occupancy, analytical chemistry, karyotyping and analysis labor;
no vendor prices or available staffing have been verified here. Stage one drug at a time
after qualification, while retaining the full stated scope and documenting deferred work.

## 7. Minimum handoff records and decisions

| Record | Minimum fields |
|---|---|
| Model registry | Material origin, model class, verified genotype/phase, editing family, clone, passage, identity/culture checks, baseline copy-number assay |
| Culture/well registry | Preparation, day/time, well, starting cells, substrate/oxygen, masked condition, randomization block, exclusions with reasons |
| Exposure table | Compound/formulation, analyte and units, nominal and measured concentrations, sampling times, stability/recovery, calibration, free/total interpretation |
| Founder/lineage ledger | Founder and parent/daughter IDs, time zero, first-event category, mitotic times, scored/uncertain error types, tracking loss, later fate |
| Orthogonal outcomes | Viable counts, functional endpoint, karyotype sampling/assay denominator and unresolved calls, absolute clone counts, treatment/washout timing |
| Analysis provenance | Frozen scoring rules, masked-code release, raw-image identity, analysis version, all planned conditions and failed runs |

Advance an upstream hypothesis only with a reproducible, genetically interpretable defect,
measured relevant exposure and joint functional improvement not explained by altered entry,
death, tracking loss or clone replacement. Better downstream function with persistent
segregation errors remains protection, not checkpoint repair. An apparently favorable result
with poorly measured exposure or wide/confounded endpoint uncertainty is inconclusive.
None of these laboratory outcomes alone would establish clinical benefit or cancer prevention.

The immediate deliverable is this reviewable plan, its
[capacity calculation — underlying file not included](../PROVENANCE.md#evidence-scope) and
[bounded source audit — underlying file not included](../PROVENANCE.md#evidence-scope). The capacity calculation
ran on 512S2 and was checked with a separate closed-form calculation. No laboratory partner
has been contacted, materials ordered or experiment represented as completed.
