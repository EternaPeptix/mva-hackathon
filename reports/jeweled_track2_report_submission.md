# jeweled — Track 2: test everolimus for BUB1B-associated muscle dysfunction

**Research submission draft, revised 21 September 2026; not submitted.** We select
**everolimus** as the single lead medication for testing whether normalization of
excessive **mTORC1 signaling** improves muscle function in a qualified BUB1B model
while preserving growth and regeneration. Recorded muscle atrophy and failure to
thrive make this a relevant, bounded objective. No patient pathway measurement,
drug-response experiment or clinical benefit has been demonstrated. Model and
laboratory access remain unestablished.

The proposed chain is **BUB1B-associated dysfunction → potentially harmful muscle
mTORC1 excess → everolimus engagement → improved absolute muscle function without
lost replacement capacity**. The middle disease-specific links are hypotheses to
test. This proposal does not claim chromosome correction, developmental reversal
or cancer prevention.

## Why this is our lead

We choose on four grounds: disease-relevant evidence, an established action of an
approved drug, a discriminating functional experiment and feasibility of further
translation. This is a qualitative research-priority decision, not a computed
probability that the drug works. The [official Track 2 criteria](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py)
reward rigor, impact, innovation and scalability.

Our preference is **narrow and conditional**. A completed secondary analysis of
the public mouse counts did not corroborate an upregulated mTORC1 transcriptional
program. We retain everolimus because the separate muscle-signaling observation
and established drug action still define a more direct approved-drug experiment
than the available alternatives. We do not claim a convergent molecular signature
or a validated best therapy.

The disease link comes from increased muscle mTORC1 signaling in heterozygous
BubR1+/X753 mice before overt muscle aging. **X753 is not L737Ter**, and mouse
L1002P models human L1012P, not this case's N1002K. No rapalog rescue was tested in
that study; abundance and signaling were measured at different ages. It nominates
a pathway, not a demonstrated therapy. [Sieben et al., Figure 5](https://www.jci.org/articles/view/126863/version/2).

Everolimus has established FKBP12/mTORC1 pharmacology and pediatric experience in
other indications. We select it pragmatically for its labeled TSC experience in
younger children and dispersible formulation. Sirolimus also offers an oral
solution and therapeutic monitoring; **neither has demonstrated comparative MVA
superiority**. [Sirolimus label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=400da0a5-2a85-44ab-a2c5-bb9161750db3).
The [everolimus prescribing information](https://www.novartis.com/us-en/product/pi/pdf/afinitor.pdf)
supports drug feasibility; it does not supply an MVA dose or safety conclusion.
Published rapamycin results are class-level support, not everolimus experiments.

The strongest argument against our selection is that human MVA expression-rescue
experiments place **usable BUBR1 restoration closer to the primary defect**.
However, choosing arimoclomol would still require both N1002K abundance-rescuability
and a demonstrated drug-to-functional-BUBR1 link. We give the established
everolimus-to-target connection greater weight for this approved-drug submission,
and concentrate the experiment on whether the nominated target is harmful and
actionable. [BUBR1 rescue precedent](https://pubmed.ncbi.nlm.nih.gov/20516114/).

## Evidence pipeline: what is established and what must be earned

| Link | Evidence available now | Required next result / failure criterion |
|---|---|---|
| **1. Candidate genotype** | Track 1 retains BUB1B Leu737Ter + Asn1002Lys, subjective exact-pair EPCR 0.8, phase unresolved and missense VUS. | Establish phase where feasible; retain alternatives. EPCR is not a probability of pathway activation or drug benefit. [Track 1](jeweled_track1_report_submission.md). |
| **2. Sequence consequence** | Completed reconstruction: L737* predicts 736 amino acids; N1002K predicts a full-length altered protein. In cis, the stop masks the downstream missense under canonical translation. | Measure endogenous RNA/protein. Predicted NMD is not measured decay; neither phase branch proves mTOR excess. [Analysis](../workflow/track2-model/README.md). |
| **3. Disease-to-pathway rationale** | Same-gene mouse muscle signaling association; our public-count reanalysis adds no positive mTORC1 transcriptional corroboration. | Demonstrate excess signaling and impaired function in matched myogenic cells with genetic controls. Normal activity in a qualified relevant model removes this lead's central rationale. |
| **4. Pathway dependence** | No case-specific causal result exists. Different-genotype muscle studies show both benefit from reducing excess and harm from excessive inhibition. | Graded partial RPTOR attenuation toward the corrected-cell range must improve function without lost differentiation/regeneration. A lower phosphoprotein signal alone fails advancement. |
| **5. Drug implementation** | Everolimus engages the nominated pathway; no BUB1B response has been measured. | Drug versus vehicle must reproduce useful function at justified, measured exposure. Genetic rescue with drug failure rejects this compound implementation. |
| **6. Durable net benefit** | Population sensitivity analysis demonstrates why favorable fractions can hide lost output. | Preserve absolute muscle function, viable differentiated output and progenitor replacement; follow after washout and through independent confirmation. No simulation score substitutes for these observations. |
| **7. Oncology and translation** | RMS history is recorded; active tumor status and regimen are unconfirmed. | Apply the actual tumor/regimen counter-screen and exposure/safety review. Credible tumor protection, treatment interference or unacceptable loss of tissue capacity blocks advancement. |

The [evidence-pipeline appendix](appendices/track2-evidence-pipeline.md) connects
sources to claims, specifies the two muscle experiments and records decision rules.
Completed computation is distinguished from proposed laboratory work.

## What the public mouse data actually add

We analysed deposited **GSE134781** counts from three mice per genotype, using
all 50 mouse Hallmark programs and three planned contrasts. The analysis retained
12,989 genes and accounted for within-program gene correlation. No program met
FDR <0.05 across the 150 tests. The mTORC1 set in X753 versus WT was weakly
downward relative to background (P=0.219; adjusted P=0.983), with median gene
log2 fold change −0.085. It did not become significant in any sensitivity analysis.
[Full analysis, figures and replay](../workflow/track2-rna/README.md).

This result **does not support the proposed pathway through this RNA signature**.
It does not negate the separate phosphorylation observations: RNA and kinase
activity are different measurements, and our analysis differs from the original
DESeq2/STRING method. Nine bulk-tissue samples with unavailable batch information
also cannot exclude a pathway effect. Neither this null result nor the generic
population simulation ranks drugs. We therefore make pathway qualification the
first decisive experiment, retain competing mechanisms and avoid presenting
literature association plus computation as independent therapeutic validation.

## The experiment that decides the proposal

Use consented patient-derived myogenic material with isogenic correction where
feasible. If unavailable, use explicitly labeled engineered models, including
the unresolved trans/cis alternatives rather than assuming phase. Fibroblasts can
help qualify allele biology but cannot settle a muscle-specific hypothesis.
Verify endogenous products, chromosome state and a reproducible defect before
treatment; retain independent editing families and preparation identities.

**Separate two stages of muscle biology.** In already mature muscle constructs,
compare absolute contractile force per construct from a fixed pretreatment cohort,
alongside viable tissue and repeated-stimulation performance. In a separate
progenitor cohort, measure expansion, differentiation/fusion and subsequent
functional muscle yield per starting cohort. Normalized force and percentages
are secondary: removing weak cells or shrinking a construct must not inflate the
primary claim. If only a surrogate contraction assay is feasible, qualify it and
state that limitation instead of claiming measured force.

Measure complementary mTORC1 substrates within matched cell states, such as
phospho/total p70S6K and 4EBP1, with nutrient, density and stimulation controls.
Bulk differences can reflect different myoblast/myotube proportions. Use graded,
inducible **partial** RPTOR attenuation as an orthogonal pathway test, with
appropriate targeting controls; aim toward the corrected-cell signaling range,
not complete ablation. Apply genetic and drug interventions to both the mature
and progenitor cohorts and to affected and corrected models.

The central alternative explanation is that increased signaling is compensatory.
Normalization with restored function supports a harmful component of excess
activity; worsening supports caution about suppressing compensation. A negative
result with failed exposure or model qualification remains inconclusive. Genetic
attenuation and everolimus are not assumed to be interchangeable.

This distinction has experimental precedent. Rapamycin improved established
myopathy in muscle-specific Tsc1-deficient mice, a different genotype.
[Castets et al.](https://pubmed.ncbi.nlm.nih.gov/23602450/). Conversely, rapamycin
impaired regenerating myofiber formation/maturation, and raptor ablation caused
muscle dystrophy. [Ge et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2793064/),
[Bentzinger et al.](https://doi.org/10.1016/j.cmet.2008.10.002).
These results require a useful signaling window and preservation of regeneration.

The next question is whether a useful effect persists after drug withdrawal.
Compare continuous exposure with a pulse followed by verified washout, recording
both target recovery and absolute functional output. Include sirolimus as the
published-compound reference and compare measured exposure, not equal nominal
concentrations. Do not presume a pulse protects growth: a randomized 2026 study
in 40 older adults found no primary functional advantage from timed weekly
sirolimus plus exercise; its sensitivity analyses favored placebo. This is
counterevidence to a generic timing claim, not a test of MVA.
[RAPA-EX-01](https://doi.org/10.1002/jcsm.70274). The
[pipeline appendix](appendices/track2-evidence-pipeline.md#compare-treatment-windows-without-presuming-a-safer-schedule)
also considers context-dependent regeneration results.

For everolimus, a pharmacologist must justify formulation, exposure duration,
medium binding, uptake and measured drug concentrations. **Whole-blood trough
values cannot be copied into culture medium as equivalent free exposure.** Keep
high mechanistic controls separate from the translational comparison. Measure
engagement and function together, including follow-up after washout. No patient
dose or useful case-specific exposure range is proposed.

Track spontaneous progenitor divisions and daughter outcomes to detect selection.
The existing [founder-ledger endpoints](appendices/founder-ledger.md) help identify
arrest and selective loss; anaphase counts do not prove successful cytokinesis or
normal karyotype. Muscle function and regenerative output are the primary benefit
criteria; reduced segregation-error frequency is not required or assumed.

Randomize within preparation/editing-family blocks and blind outcome scoring.
Use a feasibility pilot to establish maturation windows and estimate variation;
freeze meaningful benefit, allowable loss of regenerative capacity, exposure and
adverse-effect margins before independent confirmation. Cells and wells do not
create independent patients. A nonsignificant result is not a growth-preservation
or tumor-safety pass.

## What the completed simulation contributes

The [retained model](../workflow/track2-model/README.md) combines actual candidate
sequence consequences with explicit population assumptions. At one illustrative
setting, selective abnormal-cell killing improves euploid fraction by **16.8
percentage points** while differentiated-cell time falls **11.0%**. With assumed
35% slower division alone, additional per-cell function needed to break even
reaches **3.7%** and **12.5%** at the two arbitrary observation horizons in this
grid. These are counterexamples and conditional thresholds, **not measured
everolimus effects, patient timescales or tissue-force predictions**.

The model supports choosing absolute function and replacement endpoints. It does
not select a drug or establish mTOR activity from a genome. Both model layers ran
on a Mac Studio; source and results are retained for replay. No new laboratory
evidence is implied by numerical checks.

## Other routes considered and not selected as the lead

These weaknesses explain the present choice; “not selected” does not mean
biologically disproved. Positive evidence is retained beside each limitation.

| Unelected route | Strongest reason to consider it | Why not selected / what could change the decision |
|---|---|---|
| **Arimoclomol / usable BUBR1 restoration** | Selected mutant-expression rescue in human MVA biology makes this the strongest proximal alternative. [Study](https://pubmed.ncbi.nlm.nih.gov/20516114/) | N1002K instability/rescuability and drug-mediated functional rescue are both unknown. HSP90 dependence does not establish HSP70 rescue. Genetic abundance rescue followed by exposure-qualified drug rescue could displace our lead. |
| **Senolysis / dasatinib ± quercetin** | Genetic senescent-cell clearance improved selected tissue outcomes in BubR1 hypomorphic mice. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3468323/) | INK-ATTAC does not establish drug selectivity. Burden, cell-type susceptibility and replacement capacity are unmeasured; quercetin does not supply an approved-drug indication. Require selective clearance with increased net function and preserved reserve. |
| **PP2A-B56 recruitment** | Engineered kinetochore recruitment rescued alignment in human MVA cells. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3654266/) | Strong mechanistic comparison, but no demonstrated repurposed drug supplies the same spatial recruitment. Alignment does not establish sustained fidelity. Retain as a genetic rescue reference. |
| **Elamipretide / mitochondrial support** | Defined cardiolipin pharmacology and accelerated FDA approval for Barth-associated muscle strength at ≥30 kg, based on open-label knee-extensor findings. [Label](https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/215244s000lbl.pdf) | No demonstrated BUB1B mitochondrial phenotype or rescue; [randomized primary endpoints](https://www.nature.com/articles/s41436-020-01006-8) in Barth were not met. A measured mitochondrial defect with sustained functional rescue would elevate it. |
| **NAC / redox protection** | Plausible protection of already stressed cells. [Fibroblast study](https://doi.org/10.1002/cam4.101) | Genotype-dependent oxidative injury and useful exposure are unestablished. Probe interference and overlapping metabolic effects limit interpretation. A lower ROS signal does not establish tissue benefit. |
| **Readthrough / NMD manipulation** | Could recover a longer product from the stop allele. [Study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5886253/) | RNA substrate, inserted residue, protein function and exposure are unresolved. Longer protein need not work; increasing an interfering product could worsen function. First measure L737*. The claim that UGA can never yield leucine is not used. |
| **Proteasome blockade or HSP90 inhibition** | Turnover studies identify a protein-stability mechanism. [Study](https://pubmed.ncbi.nlm.nih.gov/20516114/) | Blocked degradation need not produce usable protein; HSP90 inhibition destabilized tested mutants. Broad effects and an untested N1002K response make these poor lead rescue strategies, not proven allele-specific contraindications. |
| **IL-6/IL-6R blockade** | Particular CIN-high cancer models showed an actionable dependency. [Study](https://doi.org/10.1038/s41586-022-04847-2) | Requires an actual tumor dependency; CIN alone and RMS history are insufficient. It is a separate oncology question, not the selected muscle-rescue proposal. |
| **Growth hormone** | Relevant to the growth phenotype if true deficiency exists. [Consensus](https://academic.oup.com/ejendo/article/186/6/P35/6853681) | Deficiency is unestablished; cancer predisposition generally argues against routine use, with individualized exceptions. It does not address the nominated cellular mechanism. |
| **Splice switching** | Could help if an actionable transcript defect were demonstrated. | None is established. Coding location neither proves abnormal splicing nor excludes it; an ASO needs RNA and frame/function evidence before nomination. |
| **Sirolimus, within the elected mTOR route** | The actual rapamycin compound has preclinical muscle/CIN evidence. [Muscle study](https://pubmed.ncbi.nlm.nih.gov/23602450/) | Not mechanistically rejected. Everolimus is the pragmatic lead; no head-to-head superiority is claimed. Exposure or stage-specific results could favor sirolimus instead. |

## Scope, safety and contribution

Failure to thrive, muscle atrophy, nephrocalcinosis and RMS history constrain
interpretation. Nephrocalcinosis is not a measurement of filtration. Everolimus'
infection, renal, metabolic, mucosal, pulmonary and wound-healing risks and drug
interactions require clinical pharmacology review. Pediatric experience elsewhere
does not settle this child's benefit–risk balance. The oncology stage requires
confirmed tumor/regimen information; anticancer uses do not demonstrate safety
in this RMS context. The dated [counter-screen](appendices/oncology-counterscreen.md)
supplies general design principles, not everolimus-specific evidence.

Our submission's competitive claim is a specific, falsifiable muscle-rescue
experiment connected to the candidate genotype, with computation that can expose
its weaknesses. The [Track 2 rubric](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py)
maps to concrete deliverables:

| Criterion | What this submission contributes | What remains unearned |
|---|---|---|
| Scientific rigor — 35% | Allele reconstruction, unfavorable public-data analysis, source-scoped drug comparison and genetic-versus-drug functional tests. | Case-relevant pathway causality and drug rescue. |
| Potential impact — 25% | A test of preserved muscle function and regenerative reserve, tied to recorded muscle atrophy and failure to thrive. | Magnitude or durability of benefit; neither genome correction nor clinical response is predicted. |
| Innovation — 25% | A falsifiable signaling-window experiment that separates function from selection and growth suppression, informed by genotype branches and population counterexamples. | A demonstrated useful window; using an mTOR inhibitor itself is not novel. |
| Scalability — 15% | Public-count replay and reusable genotype/pathway qualification plus absolute-output endpoints. | Transfer requires relevant models and empirical dependence in each new genotype; myogenic assays are a practical bottleneck. |

A decisive negative experiment may improve understanding; an inconclusive assay
does not automatically close a mechanism or establish impact.

The resource need is a qualified myogenic model and a partner able to measure
contraction, differentiation and exposure. Access is a prerequisite, not a
completed milestone. This is a research proposal, not a treatment recommendation.
The current [evidence pipeline](appendices/track2-evidence-pipeline.md) governs this
lead; dated abundance/NAC appendices remain historical methods.

**AI use:** Anthropic Claude Code and OpenAI Codex, personal subscriptions,
training on the user's content disabled (user-confirmed). Exact historical
model/plan names are unverified. Computation uses owned machines; candidate-level
evidence summaries are reviewed with hosted assistants. Public release and
competition submission have not occurred.

**Synthetic narration:** The existing Cartesia Sonic 3.6/Lauren movie presents the
earlier arimoclomol/NAC proposal and does not represent this elected lead. It is
retained in the private GitHub release; Vimeo is on hold. The [proposed pitch
text](jeweled_track2_pitch_reconciled.md) is unrecorded. Cartesia plan/tier and
training/retention settings remain unverified separately from Claude/Codex.
