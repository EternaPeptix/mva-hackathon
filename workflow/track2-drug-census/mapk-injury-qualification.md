# MAPK-inhibitor protection from fungal-filtrate injury

**FR180204 and SP600125 produced genuine, short-lived protection in injured A549 cells.** This earns an injury-response lead, not a BUBR1-restoration claim. SP600125 also has separate evidence of spindle-checkpoint inhibition, making it a particularly poor candidate for a faithful-division rescue strategy. SB203580 supplies a negative comparator under the tested conditions.

## Experiment and measured benefit

[Sharon et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3055868/)
tested A549 lung carcinoma cells exposed to *Aspergillus fumigatus* culture
filtrate. Cells were serum-starved for 16–18 hours, pretreated with inhibitor
for two hours, washed, and then exposed to filtrate. This is prevention of an
experimental insult, not treatment of established injury or a BUB1B-deficient
model. The filtrate contains several secreted products; heat treatment and
PrtT deletion do not identify a single injurious molecule.

| Compound | Actual exposure | What the source supports |
|---|---|---|
| **FR180204**, ERK1/2 inhibitor | 50 µM pretreatment; Calbiochem, no historical SKU/lot specified | Preserved adhesion, morphology and actin; higher XTT signal through early injury. Representative 12-hour PI-positive/annexin-V-negative events were 10.7%, versus 88.7% with wild-type fungal filtrate. |
| **SP600125**, nominal JNK inhibitor | 25 µM pretreatment; Sigma, no historical SKU/lot specified | Similar early protection, with 18.3% PI-positive/annexin-V-negative events in the representative 12-hour plot. Its pharmacology cannot be assigned exclusively to JNK. |
| **SB203580**, p38 inhibitor | 25 µM pretreatment; Calbiochem, no historical SKU/lot specified | Did not preserve the tested morphology/actin outcomes. A matched demonstration of p38-substrate inhibition was not supplied, so this is not a general disproof of p38-dependent protection. |

Figure 3B measures phospho-ERK alongside total ERK and phospho-c-Jun alongside
total c-Jun at three hours after challenge. It supports pathway modulation,
without quantitative densitometry or a genetic/drug-resistant rescue assigning
the survival response to one target.

Figure 3C shows the limitation clearly. At 12 hours, the two inhibitor arms
retain approximately 95–100% XTT signal versus approximately 20% with carrier;
by 24–36 hours, they retain only approximately 23–29%. These are visual estimates
from the plotted bars, not extracted individual measurements. Both inhibitor
bars are marked significant at 24 hours; only SP600125 is marked significant
at 36 hours. Thus substantial protection fades, but it is inaccurate to say
all benefit disappears at 24 hours. XTT measures metabolic signal, not a direct
count of surviving or faithfully dividing cells.

The paper reports three independent XTT experiments with technical triplicates.
The Methods describe error bars for replicates of a representative experiment,
while the legend describes three independent experiments. That ambiguity in
the variance denominator remains unresolved. Stars indicate P<0.05 without
exact P values or a specified multiplicity adjustment.

Figure 4B contains representative flow plots from three repeats described as
similar, not aggregate effect estimates. Methods remove filtrate and wash
before collecting adherent cells; they do not clearly include the detached
fraction. Differential recovery therefore remains a possible source of bias.
The reported fractions cannot be treated as percentages of all starting cells.
There is no long-term colony recovery, completed-daughter count, animal
drug-treatment result, normal-cell therapeutic window or muscle-function
measurement.

## What BUB1B means in this paper

BUB1B is a **fungal-challenge RNA observation**, not a drug-response protein
measurement. Supplementary Table S7 compares wild-type with PrtT-deficient
fungal filtrate after eight hours. In worksheet `Supp.Table 1G`, cells A79:E79
identify probe 7982663, BUB1B, Entrez 701, a displayed T-test P value of 0.004
and fold change 3.067484662576687. The positive change agrees with its listing
among upregulated genes in main Table 4. The table does not provide a separate
BUB1B-adjusted P value; the broader selection uses its stated FDR/filtering
procedure. It is a bulk A549 transcript comparison, with one PrtT sample
removed as an outlier, leaving three wild-type and two PrtT samples.

Neither inhibitor was administered in that transcript contrast. BUBR1 protein,
protein half-life, checkpoint competence and segregation fidelity were not
measured. The RT-PCR validation concerns IL8/MCP1, not BUB1B. This cannot be
converted into evidence that FR180204 or SP600125 restores BUBR1.

## Why SP600125 needs a separate checkpoint warning

[Schmidt et al., 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1369161/)
found that SP600125 overrides spindle-poison-induced arrest independently of
JNK. In U2OS cells, effective override began around 2.5 µM. Mps1 inhibition and
rescue by a drug-resistant Mps1 mutant support the mechanism; kinetochore
BubR1 recruitment also fell. BJ-Tert fibroblasts resisted full checkpoint
override at 10 µM, so the result is context-dependent. This is not proof of
patient-cell harm, but it is direct contrary evidence for selecting this
compound to restore checkpoint function. This secondary source was checked
through its abstract and indexed primary Results/legends; original figures
and complete Methods were not obtained. It receives no complete-paper review
or additional corpus-reading credit.

## Decision and practical availability

Keep FR180204 as an injury-response research reserve and SP600125 as a
mechanistic counter-probe with a fidelity liability. Keep SB203580's negative
result at its measured scope. None displaces the provisional restoration and
muscle-compensation candidates: the insult, tissue and outcomes do not match
the proposed MVA benefit, and no relevant exposure window is established.
Niacin and apitegromab likewise remain unproved in MVA; this is not a head-to-head
efficacy comparison or a numerical ranking.

Current manufacturer listings exist for
[FR180204/SML0320](https://www.sigmaaldrich.com/US/en/product/sigma/sml0320),
[SP600125/S5567](https://www.sigmaaldrich.com/US/en/product/sigma/s5567), and
[SB203580/S8307](https://www.sigmaaldrich.com/US/en/product/sigma/s8307).
They establish research-catalog availability, not an authenticated historical
lot, clinical formulation, safety or permission to administer them.

## Reading depth

The complete substantive PLOS main body, all five main tables and four figure
legends were read; Figures 3 and 4 were visually inspected. Supplementary Text
S1 was read as extracted text, but its embedded images were not inspected.
Only Table S7's header and BUB1B row were independently checked; no complete
supplementary-gene-table review or raw expression reanalysis is claimed.
A second reader checked the doses/timing, Figures 3–4, flow-recovery limits and
Table S7. Raw article text, images and supplementary files remain outside Git.
