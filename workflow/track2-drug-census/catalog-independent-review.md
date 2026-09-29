# Catalog-independent gene-literature intervention review

Initial review: **21 September 2026**; Amla cross-study follow-up: **29 September 2026**. This adds a bounded screen of the saved gene corpus to [the existing intervention review](literature-review.md). It is research triage, not a treatment recommendation. No drug experiment or clinical efficacy estimation was performed.

## Result and comparison boundary

**The clinical catalog did not exhaust the restorative literature.** Seven additional primary intervention routes receive explicit appraisals, alongside a deeper assay-level review of the already cited Sod1/ascorbate/NAC study. A direct chromosome-count result comes from **7,8-dihydroxyflavone (7,8-DHF) in aged mouse oocytes**. The **Phyllanthus emblica** comparison contains one BubR1 mRNA/chromosome-endpoint study and two same-group cancer/chemotherapy-context papers. **Ascorbate/NAC** improved chromosome outcomes in Sod1-deficient oocytes. These are experimental effects beyond changing a protein marker, although none is a BUB1B-MVA drug-rescue experiment.

The other studies add melatonin-associated oocyte maturation, EUK-207-associated hearing preservation, folate-associated protection from induced apoptosis, and context-specific LiCl/pifithrin observations. They deserve explicit disposition rather than disappearance at an alias gate. This is **not** evidence that they outrank niacin as a restoration research nomination or apitegromab as a muscle-compensation comparison. Their clinical exposure, population, tissue and actual-allele transfer remain unestablished. Everolimus remains a conditional downstream probe; finding more hypotheses does not validate its signaling premise.

## What was actually screened

The input was the existing `literature-gene-title-abstract.json` in `/Users/jeweled/mva-hackathon/track2-drug-census-work/`, with query `(TITLE_ABS:BUBR1 OR TITLE_ABS:BUB1B OR TITLE_ABS:"BUB-R1" OR TITLE_ABS:"mosaic variegated aneuploidy")`. This pass did not expand retrieval to additional records or follow reference chains into a new corpus.

| Step | Count | Meaning |
|---|---:|---|
| Input inventory | 1,460 unique source/ID pairs | Every record receives an inventory row; clinical names/catalog membership play no role. |
| Available abstracts | 1,432 | The remaining **28** have no saved abstract; these are unknown, not negative. |
| Broad benefit-word title/abstract filter | 422 | Titles inspected; selected sentences aid triage, not full abstract adjudication. |
| Within the above, generic intervention-word match | 236 | A priority flag only; the other **186** benefit matches were retained for title review. |
| Independent non-benefit title/context sensitivity sweep | 172 | Captures age, meiosis, muscle, compounds and related topics even without the benefit vocabulary. This found the Sod1/ascorbate/NAC study. |
| Selected titles inspected | 594 | Union of 422 and 172; **not 594 fully reviewed abstracts**. |
| All missing-abstract titles also inspected | 28 | Fifteen are outside the 594; **609 distinct titles** were seen in total. |
| Complete saved abstracts actually read | **21** | Fixed record-level ledger below and in the reproducible inventory. |
| New primary-source verification | **7 methods/results passage sets + 1 primary abstract** | Not seven complete systematic full-text reviews; figure captions and relevant methods/results were inspected. |

The exact case-insensitive filters are:

```text
benefit in title + abstract:
rescu|restor|protect|amelior|improv|alleviat|revers|attenuat|prevent|supplement

generic intervention in title + abstract (priority only):
treat|administr|supplement|compound|chemical|pharmacolog|inhibitor|activator|agonist|antagonist|extract|peptide|antioxidant|diet|vitamin|hormon

sensitivity in title, only for records without a benefit match:
oocyt|ageing|aging|senescen|muscle|proger|aneuploid|pharmacolog|drug|treat|antioxidant|diet|vitamin|hormon|supplement|compound|extract
```

The executable selection/read ledger is [gene_record_screen.py](gene_record_screen.py), with all 1,460 records in [gene-record-screen.tsv](gene-record-screen.tsv) and counts in [gene-record-screen-summary.json](gene-record-screen-summary.json). Input SHA-256: `1df717bce779851e8f8f02654cc476de320cf368ad3cbc10dad1499b5a809e90`. The script records historical reviewer actions and rejects another snapshot; rerunning filters does not manufacture a manual review. Its machine-priority counts are 236 benefit/intervention, 186 benefit-only, 172 sensitivity-title and **866 not prioritized**. Title-reading flags record the reviewer's statement that this exact selected set was inspected; the 21 complete abstract reads have an explicit PMID ledger.

The corpus includes 113 preprints, 14 patents, 13 theses, 9 PMC-only and 6 AGR records alongside 1,305 MED records; **72 records carry a review-related publication-type tag**. These categories do not all represent independent experiments. Review conclusions, docking/drug-database suggestions and title-only matches were not promoted into drug efficacy. Missing abstracts, full-text-only interventions, filter vocabulary misses, non-indexed work and unreviewed records remain coverage gaps. The 866 non-prioritized records and other unread abstracts remain **unadjudicated**, not rejected candidates.

## Added appraisals

### BDNF, MSC secretome and 7,8-DHF — a distinct chromosome-fidelity research route

In aged, non-BUB1B-mutant mouse oocytes, MSC secretome improved meiotic timing, chromosome alignment and kinetochore BubR1/ZW10 measures. Recombinant **BDNF (250 pg/mL)** reproduced timing/alignment effects; BDNF neutralization and PI3K inhibition challenged the proposed mechanism. **7,8-DHF (5 nM)** decreased in-vitro aneuploidy from **35.60 ± 2.62% to 13.16 ± 3.53%** (mean ± SD). In-vivo treatment (**5 mg/kg intraperitoneally for six days**) improved spindle/chromosome arrangement and subsequent embryo development, without significantly changing fertilization. [PMID 42286606, Figures 1–6 and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC13501573/).

**Disposition:** retain an experimental restoration comparator. The drug has a measured chromosome endpoint missing from abundance-only hypotheses; its meiotic aging model supplies no MVA muscle result or rescue of Leu737Ter/Asn1002Lys. Do not transfer the secretome/BDNF BubR1-localization assay to 7,8-DHF as though every intervention had every endpoint. In-vivo chromosome counts were not the demonstrated endpoint. BUBR1-specific causal mediation was not established by knockout/addback, and the report's TrkB-agonist designation alone does not establish target specificity. LY294002 is a mechanism-blocking control, not the beneficial intervention. The work remains a single preclinical report.

### Amla / *Phyllanthus emblica* aqueous fruit extract — chromosome endpoints and cell-state risk

The 2016 study tested a prepared aqueous fruit extract (**20–160 µg/mL**) in **NCM460 colon epithelial cells**, a spontaneously immortalized line with high background genomic instability. Several 72-hour experiments reported fewer micronuclei, nucleoplasmic bridges/nuclear buds, chromosome misalignment, lagging chromosomes and other abnormal divisions. Bulk **BubR1 mRNA** rose after 72 hours. The anaphase/metaphase ratio fell in a separate assay after four hours of exposure and one hour of release; this proxy also permits slowed mitotic progression and does not establish BUBR1 mediation. The study reported no significant apoptosis/necrosis increase within this concentration range. [Guo and Wang, PMID 27598149, Figures 1–6 and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5037716/).

Two same-group papers supply context rather than independent replication of the BubR1 finding. In 2013, water extract in **COLO320 colon adenocarcinoma cells** (20–160 µg/mL, 24–96 h) increased chromosomal-instability markers and apoptosis while delaying mitotic progression; no BubR1 assay was reported. The highest-dose 96-hour condition could not be scored for instability because too few binucleated cells remained. [Guo et al., PMID 25305222](https://doi.org/10.1024/0300-9831/a000169). In 2017, a related aqueous preparation increased chemotherapy-associated growth suppression in **Colo205 cancer cells**. In **NCM460**, it reduced measured chemotherapy-induced micronuclei and multinucleation, suppressed colony formation after reseeding, lowered mitotic index and increased apoptosis. BubR1 was not measured. The normal-cell reduction in damage burden therefore includes a possible arrest/selection effect, with potential cost to productive progenitors. [Guo et al., PMID 29204983, Figures 1–5 and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC5742286/).

These preparations, cell lines and exposures differ. The 2017 abstract/normal-cell Methods specify 0.05 µg/mL mitomycin C and 1 µg/mL cisplatin, while its tumor Figure 1 caption/body use 5 ng/mL and 10 ng/mL respectively; do not claim a matched therapeutic window. Extract alone at 320 µg/mL inhibited NCM460 growth by about 28% in the 2017 study, so the 2016 concentration-range finding is not a general safety result. No active component was isolated, and no matched/isogenic cancer comparison, MVA-specific cancer surveillance or muscle efficacy was demonstrated. MVA cells could fall on an adverse-response side of this cell-state-dependent hypothesis.

**Competitive disposition:** admit a **high-priority laboratory comparator among reviewed consumer preparations**. One study measures a BubR1 transcript response and chromosome-related endpoints in the same cell line; this adds a different observation to niacin's cross-intervention NAD/SIRT2 rationale without proving functional protein restoration. Niacin separately has human muscle NAD engagement and a defined medicine; apitegromab has pediatric motor-function evidence outside MVA. Test authenticated, compositionally characterized extract lots beside niacin and NMN in affected/corrected myogenic cells, measuring functional full-length BUBR1, new segregation errors, division rate, absolute muscle output, viability and progenitor renewal. Use a matched tumor counterscreen. Oral exposure and equivalence to retail products remain unknown, and individual constituents cannot inherit the mixture's findings.

**Source scope:** this 29 September follow-up checked the 2013 publisher PDF, 2016 full primary article and 2017 full primary article/publisher PDF. It expands context for the existing Phyllanthus decision; it adds no new ledger entry or original-queue reading count. Nocodazole and cytochalasin B remain challenge/assay reagents.

### Ascorbic acid and NAC — chromosome protection in a different genetic lesion

In **Sod1-knockout mouse oocytes**, in-vitro maturation under 21% oxygen was associated with decreased kinetochore BUBR1 and chromosome-segregation abnormalities. Adding **500 µM L-ascorbic acid or 500 µM NAC** increased maturation and normal 20-chromosome outcomes, and reduced premature sister-chromatid separation. The ≤19/≥21-chromosome categories were not each significantly different among treatment groups. Figure 7 lists 50 control, 39 ascorbate and 52 NAC oocytes. [PMID 39845481, Figure 7](https://pmc.ncbi.nlm.nih.gov/articles/PMC11751902/).

**Disposition:** strengthen the contextual antioxidant comparator beyond the previously cited NAC protein-abundance observation. Sod1 deficiency and oxygen-rich culture impose a different causal lesion from BUB1B-MVA. The treatment experiment does not establish that restoration of kinetochore BUBR1 mediates its benefit or that the endogenous MVA variant remains functional. No muscle or clinical exposure bridge is demonstrated. NAC is an existing decision route with an additional source, not a new unique drug merely because this source was missed.

### Melatonin — functional cryopreservation protection, stage-specific marker evidence

**10⁻⁷ M melatonin** across vitrification, warming and maturation increased mouse germinal-vesicle oocyte maturation to MII (**72.67% versus 61.30%** among morphologically normal oocytes, not all vitrified inputs) and improved spindle morphology, redox and mitochondrial measures. Other tested concentrations (10⁻⁹, 10⁻⁵ and 10⁻³ M) had lower maturation than vitrification alone. Bubr1 mRNA recovery was demonstrated at the GV stage; it was not a uniform all-stage correction. [PMID 31480299, Table 2 and Figures 1–6](https://pmc.ncbi.nlm.nih.gov/articles/PMC6770451/).

**Disposition:** a contextual cytoprotection reserve with a narrow concentration-dependent experimental effect. The relevant paper did not demonstrate BUB1B-mutant rescue, chromosome-number correction, muscle function or BUBR1-dependent mediation. Embryo/blastocyst benefits discussed from earlier publications are not outcomes of this same GV-oocyte experiment. This is stronger than marker-only evidence for oocyte maturation, but weaker than a chromosome-fidelity rescue claim.

### EUK-207 — organismal function in a non-BUB1B aging model

The synthetic SOD/catalase mimetic **EUK-207** attenuated hydrogen-peroxide-induced DNA-damage/senescence measures in cochlear explants. In **SAMP8** mice, **0.2 mg/kg/day** by subcutaneous pump over three months reduced hearing-threshold deterioration and hair-cell loss: mean final auditory brainstem response thresholds were **55.7 ± 3.5 versus 75.1 ± 4.7 dB** (mean ± SEM) for treatment and vehicle. The methods allocated ten mice per group. [PMID 30693443, Figures 5 and 9](https://pmc.ncbi.nlm.nih.gov/articles/PMC6614136/).

**Disposition:** a functional antioxidant/hearing-protection reserve. The paper's BubR1 decreases belong to the oxidant-exposure and aging characterization; do not automatically attribute restored BubR1 to the drug. No BubR1-mutant animals, BUBR1-dependence test, muscle rescue or MVA exposure bridge is shown. The organismal functional effect deserves acknowledgment, but it does not establish general antioxidant equivalence or compete directly with human muscle evidence.

### Folic acid — apoptosis protection with a weak BUB1B bridge

In bovine **Mac-T mammary epithelial cells**, **20 µM folic acid for 18 hours** reduced apoptosis induced by 400 µg/mL S. aureus lipoteichoic acid (**p=0.034**). BUB1B emerged among cell-cycle/network targets in transcriptomic analyses. Cow blood-transcriptome trends were directionally similar, but the paper explicitly states those individual-level gene differences were **not significant**, and the subclinical mastitis cows were not confirmed to have S. aureus infection. [PMID 39506684, Figures 1 and 7 and discussion](https://pmc.ncbi.nlm.nih.gov/articles/PMC11542246/).

**Disposition:** retain an evidence-record entry, with low priority for the BUB1B restoration comparison. The cellular anti-apoptotic effect is real evidence in its model; BUB1B is a proposed associated target, not an established mediator. Less apoptosis without division-fidelity or productive-renewal measurements cannot establish MVA progenitor rescue. These results are neither a demonstrated cow disease-cure experiment nor a human MVA treatment effect.

### LiCl — pharmacological rescue of Axin1 depletion, not the FRZB result

After **Axin1 siRNA** in mouse oocytes, **LiCl** improved spindle morphology and pronuclear formation. The methods report 2.5 µM LiCl in culture. BubR1 retention at kinetochores characterized the defective, checkpoint-activated oocytes; it was not a BUB1B mutation being repaired. [PMID 27284927, Figures 2 and 4](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0157197).

**Disposition:** a contextual Wnt/cytoskeleton mechanism observation, not a pharmacological implementation of the existing BubR1-H/H **sFRP3/FRZB genetic** rescue. It cannot supply generic lithium efficacy for MVA or an exposure assumption for a clinical lithium formulation. Keep the LiCl experiment distinct from that genetic opportunity.

### Pifithrin-alpha and NAC — arrest relief is not proven genomic repair

The accessible primary abstract reports that NAC rescued chromium-induced ATP depletion/cell-cycle arrest in L-02 hepatocytes. **Pifithrin-alpha** countered reduction of Mrc1 and BubR1 and inhibited cell-cycle arrest. Full-text treatment details were not reviewed; verification here is limited to the publisher's abstract. [PMID 22886373](https://link.springer.com/article/10.1007/s11010-012-1423-7).

**Disposition:** record contextual pharmacology without promoting p53 inhibition to a restoration lead. Releasing a damage response is not evidence of correct division or durable competent progenitors. This qualifies any broad wording that there is “no pifithrin drug rescue”: an imposed-toxicant cellular rescue is reported, while **BUB1B-MVA functional rescue remains unshown**. It does not overturn the same-gene p53 genetic counterevidence already appraised.

## Reconciliation and read ledger

The eight PMIDs above were absent from `decisions.json` and `literature-review.md` at the start of this pass. However, the Sod1/ascorbate/NAC source was already cited and summarized in `decision-report.md`, so it is a deeper review of existing evidence, not a newly discovered study. Melatonin occurred as a search term, and pifithrin/lithium routes were discussed without these experiments. Existing NAC, p53 and Wnt routes are refined rather than counted again. The eight sources correspond to different preparations/interventions and are **not an eight-drug unique-entity census**.

The **21** completely read saved abstracts are PMIDs: `42286606`, `27598149`, `31480299`, `39506684`, `30693443`, `27284927`, `22886373`, `40470792`, `37201634`, `36539842`, `32078082`, `19846757`, `40750410`, `37782908`, `39845481`, `21550059`, `26883501`, `28229533`, `41579174`, `42000973`, `29067790`.

The other thirteen full-abstract readings were used to distinguish the intervention:

| Records | Abstract-level disposition |
|---|---|
| 40470792, 36539842, 19846757, 29067790 | BUBR1, HtrA2, NAMPT or SIRT2/BUBR1-variant expression manipulations; do not assign genetic rescue to an unadministered drug. |
| 37201634 | Aristolochic acid I is the damaging exposure. |
| 32078082 | Polyamine/putrescine withdrawal in arginase-deficient CHO culture; not evidence of MVA rescue by supplementation. |
| 40750410, 37782908 | Propargylglycine blocks engineered tumor-cell killing; PAF(C-16) restores malignant behaviors after BUB1B knockdown. These are not therapeutic restoration of normal MVA tissues. |
| 21550059 | ROS signaling and BubR1 in aged vascular smooth muscle; kinase/oxidase inhibitors blocked induced expression rather than demonstrating restoration. |
| 26883501, 28229533 | BubR1-H/H barrier/MDSC phenotypes; the NF-κB genetic reduction result belongs to ERCC1 mice. |
| 41579174, 42000973 | Gentamicin-injury transcriptomic/drug-database proposals and MLH3-loss pathology, respectively; no protective drug experiment established by these abstracts. |

A fair competition experiment should admit 7,8-DHF and defined antioxidant/preparation hypotheses to the same actual-allele, chromosome-fidelity, productive-renewal and absolute-muscle-output qualification used for niacin. Different dimensions can favor different candidates: an oocyte chromosome assay is informative for restoration, while human muscle-function evidence informs compensation. Neither dimension alone settles the medication choice. This review establishes additional testable routes within the saved corpus; it does not establish clinical benefit, exhaustive mechanism coverage or a uniquely best medication pair.
