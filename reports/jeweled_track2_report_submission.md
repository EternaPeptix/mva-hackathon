# jeweled — Track 2: mTORC1 normalization versus usable BUBR1 restoration for muscle renewal

**Research submission draft, revised 21 September 2026; not submitted.**

The case's [recorded muscle atrophy and failure to thrive](../workflow/configs/lirical-subject.json)
motivate a concrete target: preserve muscle function and its capacity for renewal.
Our [genomic hypothesis](jeweled_track1_report_submission.md) is BUB1B
p.(Leu737Ter) plus p.(Asn1002Lys); phase is unresolved and the missense remains VUS.
The proposed mechanism tests therefore establish how the surviving allele products
affect muscle before interpreting a drug response.

We nominate
**everolimus as the narrow lead** and **niacin (nicotinic acid) as the primary
orthogonal challenger**. Test them separately, not in combination. Arimoclomol is
a restoration reserve; dasatinib is a reserve requiring harmful senescent burden
and selective susceptibility. These are approved-medication hypotheses, not
established treatments or a claim that one drug is universally best.

The comparison distinguishes two useful outcomes: **downstream muscle gain despite
persistent segregation defects**, or **restored BUBR1 function that reduces new
errors and sustains productive progenitors**. Both face the same standard:
**greater absolute muscle output with preserved progenitor reserve and no worsened
chromosome fidelity**. This makes the proposed decision concrete: advance a useful
drug implementation, identify which mechanism it supports, or reject a response
that improves a marker while losing productive tissue. Neither route is assumed
to repair existing aneuploid genomes. Laboratory access and patient-derived
experiments remain unestablished; the work below specifies future tests.

## Why this focused set is the submission choice

Everolimus leads because it implements its proposed target through established
pharmacology and has pediatric formulation, exposure and monitoring experience.
Niacin has human muscle NAD engagement and the stronger proximal causal rationale,
but niacin-to-SIRT2-to-functional-mutant-BUBR1 remains unproved. **Prioritizing causal
restoration instead would put niacin first.** Including these two distinct routes
is more robust than their ordering. Genetic intervention evidence is not credited
to either drug, and no weighted score or therapeutic probability is assigned.

| Same criterion | Everolimus | Niacin |
|---|---|---|
| Disease connection | Same-gene muscle mTORC1 signaling association in a different truncation model; causality unproved. | Human MVA expression rescue supports usable BUBR1 restoration for selected mutants; no niacin or N1002K rescue. |
| Named-drug implementation | Established mTORC1 inhibition; useful normalization in affected muscle remains unproved. | Human muscle NAD response; utilization, SIRT2 engagement and functional BUBR1 rescue must be distinguished. |
| Phenotype relevance | No demonstrated BUB1B rapalog functional rescue; weak negative and BUBR1-loss evidence remain relevant. | Small, uncontrolled other-disease strength observations; genetic SIRT2 did not establish skeletal-muscle rescue. |
| Pediatric feasibility | Approved pediatric uses and characterized monitoring, with material toxicity and interaction constraints. | Approved prescription product, but pediatric efficacy, safety and PK are unestablished. |
| Result that advances it | Exposure-qualified improvement in absolute function with preserved renewal and no worse fidelity; causal pathway evidence additionally supports the mTOR interpretation. | The same functional and preservation standard; BUBR1 dependence and improved fidelity additionally validate the restoration claim. |

Recorded failure to thrive and muscle atrophy make loss of growth or renewal a
decisive failure for everolimus. Pediatric pharmacology establishes feasibility,
not MVA safety. Nephrocalcinosis is not a measured filtration deficit; current renal
function, tumor status and actual oncology regimen remain unconfirmed. Those facts
constrain translation without preventing the present hypothesis nomination.

The [official criteria](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py)
permit candidate-medication hypotheses. Completed biological proof is not a
prerequisite for this submission. Selective senescent-cell clearance has the
strongest demonstrated same-gene tissue/function mechanism, but its genetic result
does not establish an eligible drug implementation. Burden qualification can
proceed alongside the two main routes; it need not wait for their failure.

## Lead: everolimus to test harmful mTORC1 excess

The [BubR1+/X753 study](https://www.jci.org/articles/view/126863/version/2)
found increased phosphorylation of p70 S6 kinase and 4EBP1 in gastrocnemius from
three-month-old mice (two mice per genotype in the signaling blot). This is direct
pathway-activity evidence, but small and not a drug rescue. X753 is not L737Ter;
mouse L1002P models human L1012P, not N1002K. Abundance and signaling were measured
at different ages. Our [public mouse RNA reanalysis](../workflow/track2-rna/README.md)
found no positive mTORC1 Hallmark corroboration (X753–WT P=0.219, global FDR=0.983).
This adds no transcriptional support; it does not rebut the phosphorylation result.

**A specific upstream hypothesis is altered insulin-receptor (IR) trafficking.**
[Choi 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5560052/) established a
nonmitotic BUBR1-dependent endocytosis mechanism and increased insulin sensitivity
in BubR1-hypomorphic mice. [Choi 2019](https://www.nature.com/articles/s41467-019-09318-3)
extended the trafficking mechanism; SHP2 inhibition prolonged AKT signaling in
liver, but not skeletal muscle. These findings do not establish harmful mTORC1
activation, an N1002K effect or a muscle-disease mechanism. If myogenic mTORC1 excess
is found, compare endogenous IR internalization and signaling duration in affected
and corrected cells. The [focused follow-up](appendices/track2-evidence-pipeline.md#conditional-follow-up-insulin-receptor-trafficking)
tests the proposed chain **BUBR1 dysfunction → delayed IR uptake → prolonged
signaling → harmful mTORC1 activity**. Each arrow requires evidence; improved
insulin sensitivity alone is not pathology. The Sieben study already proposed
this connection; our contribution is a proposed allele-aware muscle test, not
discovery of the link. This refines the existing hypothesis without changing the medication
shortlist or making the experiment a submission prerequisite.

The [everolimus label](https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=2150f73a-179b-4afc-b8ce-67c85cc72f04)
describes FKBP12-mediated mTORC1 inhibition and pediatric TSC indications from age
one for SEGA and age two for associated seizures. Whole-blood concentration
monitoring and an oral-suspension formulation support implementation feasibility.
Their indication-specific exposure ranges are neither MVA doses nor free culture
concentrations. Established target inhibition does not establish beneficial
normalization, preserved myogenesis or acceptable systemic effects in this child.

**Why it leads over the earlier arimoclomol idea:** everolimus has an established
drug-to-target link and a muscle signaling observation in a BUB1B model. The
arimoclomol restoration idea is closer to the primary protein defect, but
[HSP90 dependence and expression rescue of selected mutants](https://pubmed.ncbi.nlm.nih.gov/20516114/)
do not demonstrate arimoclomol-mediated rescue. Both drugs have pediatric approved
uses elsewhere; approval alone does not distinguish them. The [arimoclomol exposure
review](appendices/mechanistic-bridge.md) adds an unresolved translation gap.
This is a narrow advantage in the justification for testing a named drug, not
demonstrated MVA efficacy or greater potential impact. Arimoclomol would move ahead
with endogenous mutant-function and muscle rescue at justified exposure, especially
if everolimus cannot preserve productive renewal.

Counterevidence remains explicit. A [2013 BubR1 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3785294/)
reports preliminary unpublished observations of no apparent age-phenotype benefit
with rapamycin, without treatment details sufficient for a definitive negative.
[Goutas et al., Fig. 6](https://pmc.ncbi.nlm.nih.gov/articles/PMC10149375/) found lower
BUBR1 after 50 nM everolimus for 24 hours in late-passage human Wharton-jelly MSCs,
but not early-passage cells. These were not BUB1B-mutant cells. This warns of
state-dependent effects; it does not establish MVA harm or negate the hypothesis.

First test genotype-associated excess signaling at matched cell state and whether
partial RPTOR attenuation toward the corrected-cell range improves function.
Then determine whether everolimus implements that benefit at justified exposure.
Measure BUBR1, productive divisions and chromosome fidelity alongside signaling
and mature-muscle output. The [pipeline](appendices/track2-evidence-pipeline.md)
retains sirolimus as the within-class reference and compares continuous treatment
with pulses followed by verified washout; intermittent inhibition is not presumed
to preserve growth. Absence of harmful excess or functional failure despite
engagement would move this hypothesis down.

## Orthogonal challenger: niacin to restore usable BUBR1

[Human MVA expression-rescue experiments](https://pubmed.ncbi.nlm.nih.gov/20516114/)
support residual-function rescue for selected mutants, not niacin or N1002K.
[North et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/) connected NAD/SIRT2
to BUBR1 stability. NMN increased BUBR1 in cells and aged wild-type mice; the
hypomorphic-mouse functional/lifespan intervention was **genetic SIRT2 elevation**,
not niacin or NMN treatment. Genetic SIRT2 did not correct small size, cataracts or
lordokyphosis; this is not demonstrated skeletal-muscle rescue.

[Pirinen et al.](https://pubmed.ncbi.nlm.nih.gov/32386566/) studied niacin in five
adults with mitochondrial myopathy and ten controls. Muscle NAD restoration and
strength changes in this small, open-label, nonrandomized study support human muscle
target engagement, not BUB1B efficacy. In a [separate pediatric myopathy case](https://pubmed.ncbi.nlm.nih.gov/39173541/),
NAD correction accompanied little clinical response and continued progression.
**Niacin → functional BUBR1** and **this allele → rescuable function** remain unproved.
NMN is a mechanistic reference; niacin, nicotinamide/niacinamide and NMN are distinct
interventions and cannot inherit one another's drug effects or approval status.

[Human-myotube experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC10963276/)
show modest, context-dependent NAD responses to nicotinic acid, with different
statistical support across panels. The [pipeline](appendices/track2-evidence-pipeline.md)
retains exact results and requires NAPRT/downstream precursor-utilization tests in
each intended cell state. A parent-drug culture failure does not alone exclude a
systemic metabolite mechanism, but hypothesized conversion cannot count as success.
Measure baseline NAD without manufacturing deficiency to obtain a rescue signal.

[Prescription niacin labeling](https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=dd60fe81-3d54-4a29-9502-601be89c575e)
establishes an approved product, with pediatric efficacy, safety and PK unestablished.
Preparations cannot be silently substituted. No useful culture exposure or child
dose is inferred from the adult study. Niacin moves ahead if qualified exposure
produces useful endogenous BUBR1-dependent function; NAD elevation alone cannot pass.

## The future experiment that distinguishes the hypotheses

| Stage | Existing evidence | Required discriminator |
|---|---|---|
| **Genotype and model** | BUB1B Leu737Ter + Asn1002Lys; subjective exact-pair EPCR 0.8, phase unresolved, missense VUS. [Track 1](jeweled_track1_report_submission.md). | Matched affected/corrected preparations with reproducible defects; establish phase where feasible. EPCR is not a drug-response probability. |
| **Endogenous products** | Reconstruction predicts a 736-amino-acid truncation and full-length N1002K. In cis, the stop masks the downstream substitution under canonical translation. | Allele-specific RNA, full-length/truncated protein, localization and turnover. Predicted NMD is not measured decay; cis does not establish half-normal usable protein. |
| **Mechanism qualification** | mTOR association and selected-mutant expression rescue supply different precedents. | Partial RPTOR normalization; separate physiological WT/N1002K/actual-truncation titrations to identify helpful, inert or interfering products. Measure function at matched abundance. |
| **Separate drug arms** | Everolimus has established target pharmacology; niacin has human muscle NAD engagement. | Vehicle versus everolimus versus niacin at separately justified concentrations/durations, with corrected controls and exposure/engagement measurements. No combination assumption. |
| **Function and durability** | No biological drug result for this genotype. | Absolute output, preserved progenitor reserve and no worsened fidelity; independent confirmation, withdrawal follow-up and actual tumor/regimen assessment. |

Engineered models must retain unresolved phase alternatives and be labeled as
engineered. Fibroblasts can qualify allele biology but cannot settle muscle benefit.
Increased protein alone is insufficient: additional N1002K may remain defective,
while chaperone-mediated improvement per molecule remains a distinct reserve
hypothesis. N1002K-specific SIRT2 binding is unmeasured; sequence position establishes
neither preserved nor disrupted interaction.

Use two functional cohorts. Follow progenitors through productive divisions,
differentiation and replacement, measuring absolute contractile output per starting
cohort. Separately measure established mature-muscle force and repeated-stimulation
performance from its pretreatment cohort. Area-normalized force, fusion percentages
and survivor fractions are secondary. **Meaningful improvement in either cohort
can support advancement, provided output in the other cohort, progenitor renewal
and chromosome fidelity meet the prespecified preservation margins.** Improvement
in both cohorts is not required; benefit in one cannot offset unacceptable loss
in the other. This rule applies equally to both drugs.

Retain death, arrest, non-entry, division time, clone expansion and tracking loss.
The [founder ledger](appendices/founder-ledger.md) detects selection; it is not an
independent efficacy result. Fewer observed segregation errors in a slower culture
cannot establish rescue, and apparently normal anaphases do not prove normal
karyotypes. Verify chromosome outcomes independently.

| Mechanism-specific result | Interpretation |
|---|---|
| Everolimus normalizes qualified excess and improves absolute function with preserved reserve/fidelity | Supports the mTOR hypothesis; correcting the primary segregation defect is not required. |
| Niacin improves usable BUBR1, spontaneous fidelity and productive output with interpretable SIRT2/BUBR1 dependence | Supports the restoration hypothesis and can reverse the lead order. |
| Useful function improves without the claimed pathway/BUBR1 dependence | Retain the observation; investigate a different mechanism rather than claiming the proposed chain. |
| Genetic rescue works, but a drug fails despite adequate exposure and engagement | Target success, implementation failure within the tested scope. |
| Markers improve without useful function, or benefit requires unacceptable reserve/fidelity loss or unjustified exposure | Do not advance that implementation. |
| Model, exposure, precursor utilization or functional assay is inadequate | Uninterpretable nonresponse; repair the limitation or report infeasibility. |

Randomize and mask within independent preparation/editing families. After a pilot,
freeze useful-effect, allowable-loss, missing-data and multiplicity rules before
confirmation; cells and images are nested observations, not independent replicates.
The four possible drug-by-benefit-cohort claims form one confirmatory family;
each claim requires both benefit and preservation. The pipeline specifies the
multiple-comparison rule and freezes one regimen per drug before confirmation.
A resolved negative requires adequate precision and applies only to the tested
implementation, tissue, genotype, exposure and duration. The [governing pipeline](appendices/track2-evidence-pipeline.md)
specifies these future stages; missing laboratory evidence does not defer the
present two-medication hypothesis choice.

## What the completed computation contributes

The [genome-informed model](../workflow/track2-model/README.md) reconstructs candidate
sequence consequences and explores declared population assumptions. One illustrative
selective-killing scenario raises euploid fraction by 16.8 percentage points while
differentiated-cell time falls 11.0%. This is a completed assay-design result:
a favorable cell percentage can coexist with less modeled tissue output, so the
proposed comparison retains each starting cohort and measures absolute function.
The values are synthetic counterexamples, not drug effects or patient timescales.

The [GSE134781 analysis](../workflow/track2-rna/README.md) includes nine mice,
12,989 retained genes and 50 Hallmark programs across three contrasts. No program
met global FDR <0.05. This secondary analysis is not independent biological
replication or a drug ranking. Both completed analyses ran on a Mac Studio.

**AI use and data handling:** Claude Code/Codex personal subscriptions, with
training disabled (user-confirmed); exact historical model/plan names are
unverified. Controlled raw case reads and genome-scale payloads remain excluded
from this private repository. [Full methods disclosure](jeweled_methods_form_answers.md).

## Other routes considered and not selected as the lead

| Reserve / alternative | Strongest reason to consider it | Weakness and condition for reconsideration |
|---|---|---|
| **Dasatinib reserve / senolysis** | [INK-ATTAC](https://www.nature.com/articles/nature10600) improves selected BubR1-hypomorph tissue outcomes. [Patent Example 11](https://patents.google.com/patent/WO2015116735A1/en) reports D+Q adipose-marker reductions in that model. | Engineered clearance differs from a drug. Patent evidence supplies no functional rescue, sample size or doses; p16 reduction was not conventionally significant. Peer-reviewed D+Q functional evidence concerns [Ercc1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4531078/). [Dasatinib activity in human fat-cell progenitors](https://pubmed.ncbi.nlm.nih.gov/25754370/) does not establish muscle-supporting-cell susceptibility, and monotherapy cannot inherit combination efficacy. D+Q eligibility remains unresolved because quercetin is not an approved drug; the rules do not expressly require an approved regimen. Qualify harmful burden, selective susceptibility and increased absolute output. |
| **Sirolimus — within-class alternative** | Actual rapamycin reference; [rapamycin improved Tsc1-deficient muscle](https://pubmed.ncbi.nlm.nih.gov/23602450/). | No demonstrated MVA superiority over everolimus. Retain as a pharmacologic reference; require an exposure-qualified functional window rather than transferring rescue from Tsc1 myopathy. |
| **Arimoclomol — restoration reserve** | [Pediatric NPC approval with miglustat](https://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=5feffc0e-453d-47fa-91dd-38d4952309bc&version=5) and [other-protein processing results](https://doi.org/10.1016/j.ymgme.2025.109103). | No demonstrated BUBR1 rescue. HSP90 dependence of tested MVA mutants does not establish HSP70 rescue of N1002K. The [exposure appraisal](appendices/mechanistic-bridge.md) records substantial concentration and compartment gaps. Promote if it improves endogenous mutant function at justified exposure; failed quantity-only rescue does not exclude improved folding/function per molecule. |
| **Acipimox / alternative NAD-support implementation** | [Human mitochondrial effects](https://pmc.ncbi.nlm.nih.gov/articles/PMC4375076/) and an [approved lipid indication](https://www.medicines.org.uk/emc/product/1357/smpc). | Human muscle NAD/SIRT2/BUBR1 engagement is unestablished; direct NAD evidence used 10 mM mouse myotubes. A [two-week older-adult trial](https://doi.org/10.1515/teb-2024-0030) showed no functional or muscle-protein-synthesis benefit; the separate diabetes trial showed rebound fatty acids and worse insulin sensitivity. These are limited, different-disease studies, not an MVA exclusion. |
| **Sodium phenylbutyrate** | Approved products and [functional rescue of mutant ABCA1](https://pmc.ncbi.nlm.nih.gov/articles/PMC3520540/). | Different protein; no located BUBR1-specific rescue. That experiment used 10 mM for 24 hours, with substantial exposure translation unresolved. Broad chaperone/HDAC activity does not select a BUBR1 mechanism. |
| **Teprenone / geranylgeranylacetone** | Japan-market-approved drug; a [13-versus-13 human atrial study](https://pubmed.ncbi.nlm.nih.gov/31302249/) found HSPA1/HSPB1 induction. [PMDA product record](https://www.pmda.go.jp/PmdaSearch/rdDetail/iyaku/2329012C1026_1?user=1). | HSP90 was not measured, nor BUBR1 rescue; atrial target engagement is not myogenic function. The Japanese label reports no pediatric clinical trials. It does not displace arimoclomol without a functional BUBR1 result. |
| **PP2A-B56 recruitment** | [Engineered kinetochore targeting](https://pmc.ncbi.nlm.nih.gov/articles/PMC3654266/) rescued alignment in MVA cells. | Strong genetic comparator; no demonstrated approved drug supplies the same spatial recruitment. Alignment alone is not sustained chromosome fidelity. |
| **Elamipretide** | [Accelerated FDA approval](https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/215244s000lbl.pdf) for Barth-associated muscle strength at ≥30 kg. | [Randomized primary endpoints](https://www.nature.com/articles/s41436-020-01006-8) were not met; no BUB1B-specific mitochondrial defect or rescue. Requires a measured defect and useful sustained function. |
| **NAC** | Plausible protection of stressed cells; [fibroblast precedent](https://doi.org/10.1002/cam4.101). | No demonstrated genotype-dependent redox mechanism or useful exposure. Lower ROS, with possible probe interference, cannot establish rescue. |
| **NMD inhibition / transcript stabilization** | Could expose a useful truncated product if endogenous RNA decays. [Amlexanox experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC3562214/) provide a different-gene precedent. | Measure Leu737Ter RNA and actual 1–736 product function. Historical "731X" is a different frameshift product; preserved arrest does not establish normal alignment/silencing. Amlexanox also affects readthrough, so separate RNA stabilization from the resulting protein species. No BUB1B rescue or useful systemic exposure established. |
| **Stop-codon readthrough** | Could recover a longer product; [readthrough study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5886253/). | Requires available transcript, sufficient readthrough and a functional resulting sequence. NMD reduces substrate but does not establish zero substrate. This differs from retaining a useful truncation; neither strategy is currently qualified. We do not use the absolute claim that UGA can never yield leucine. |
| **Proteasome blockade / HSP90 inhibition / checkpoint inhibition** | Turnover and attachment mechanisms identify possible intervention points. [MVA study](https://pubmed.ncbi.nlm.nih.gov/20516114/). | Accumulating protein need not restore function; HSP90 inhibition destabilized tested mutants. Checkpoint suppression conflicts with the proposed nonmalignant fidelity objective. These are poor rescue leads, not proven allele-specific clinical contraindications. |
| **Chloroquine / hydroxychloroquine / autophagy blockade** | [Tang 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3532042/) found chloroquine-sensitive trisomic fibroblasts; [Goutas 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10149375/) found increased BUBR1 with hydroxychloroquine in late-passage WJ-MSCs. | Tang's BUBR1-low experiment tested AICAR/17-AAG, not chloroquine; chloroquine did not suppress the tested human CIN/MIN lines. Neither transfers to hydroxychloroquine muscle rescue. The high-concentration Goutas abundance observation also lacks functional rescue. Selective suppression must improve absolute tissue output before reconsideration. |
| **Metformin / selective energy stress** | [Tang 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3532042/) reported weaker suppression of trisomy 13/16 cell accumulation than with AICAR. | An antiproliferation result, not BUB1B restoration or muscle-function rescue. Requires relevant selectivity and net tissue benefit; pathway membership cannot supply them. |
| **IL-6/IL-6R blockade** | [CIN-high cancer experiments](https://doi.org/10.1038/s41586-022-04847-2) identify a dependency in particular models. | Requires an actual tumor dependency; RMS history and CIN alone are insufficient. Separate oncology question. |
| **Growth hormone** | Potential relevance if true hormone deficiency exists. | Deficiency is unestablished; [consensus](https://academic.oup.com/ejendo/article/186/6/P35/6853681) generally discourages use in cancer-predisposition syndromes, with individualized exceptions. No nominated cellular rescue mechanism. |
| **Splice switching / gene-specific experimental agents** | Could address an observed actionable transcript defect. | No such defect is established; an ASO or novel peptide is not an existing approved medication for this proposal. Coding location neither proves nor excludes abnormal splicing. |

## Scope, translation and competitive contribution

This submission nominates hypotheses with a future ex-vivo validation plan.
Everolimus's immune, metabolic, growth/renewal and interaction constraints, niacin's
hepatic/metabolic/muscle risks and possible support of abnormal clones require
candidate-specific assessment. The [current oncology evidence and counter-screen](appendices/track2-evidence-pipeline.md#rhabdomyosarcoma-evidence-and-treatment-compatibility)
distinguish mixed temsirolimus trial results from limited direct everolimus data.
Neither supports assumed cancer protection or compatibility with the child's
unconfirmed regimen. No patient dose, developmental reversal or clinical efficacy
is claimed.

| Official criterion | Submission contribution | Evidence still needed |
|---|---|---|
| Scientific rigor — 35% | Two distinct mechanisms, actual drug-implementation gaps, counterevidence and scoped stopping rules. | Genotype-relevant mechanism and useful drug response. |
| Potential impact — 25% | Preserve absolute muscle function and productive progenitor reserve. | Magnitude, durability and translation of benefit. |
| Innovation — 25% | Mechanism-specific perturbations and common functional accounting distinguish normalization, restoration and false rescue. | Demonstrated mechanism; familiar drug repurposing alone is not claimed novel. |
| Scalability — 15% | Replayable sequence/public-data analyses and reusable mechanism/drug decision rules. | Qualification in each new allele, tissue and disease. |

The focused choice is settled for submission; biological efficacy remains open.
A consented model, myogenic functional-assay laboratory and pharmacology support
are needed for the partner-ready two-drug study. Informative negative results can
reject scoped hypotheses; failed assays cannot. Media remains on hold during review.
