# What reading the gene literature adds

**21 September 2026. All 1,432 available saved abstracts were read.** The source-order review reads complete saved abstracts
without a catalog-name or machine-priority filter. Its current coverage is in
[`gene-abstract-review-summary.json`](gene-abstract-review-summary.json), with an
individual assessment for every source record in
[`gene-abstract-review.tsv`](gene-abstract-review.tsv). Authored notes are retained
separately from the compiler; source changes cannot inherit historical reading.
There are 226 records primarily classified as administered-drug studies and 231
follow-up flags across all 1,460 records; these are not counts of independent
experiments, restorative drugs or new leads. An abstract lacking a drug mention
is not proof that the full paper lacks one.
Saved thesis abstracts ending in “continues” remain incomplete source texts.

## New leads and their competitive position

| Route | Positive evidence actually reported | What separates it from the current nominations | Decision |
|---|---|---|---|
| **Geraniin** | In noncancerous CCD841 colon cells, treatment raised checkpoint-gene transcripts including BubR1, strengthened checkpoint behavior and reduced mitotic aberrations/CIN. HCT116 cancer cells showed the opposite checkpoint direction and increased damage/death. | Actual chromosome-fidelity pharmacology is closer to our safety endpoint than a protein-expression association. However, the abstract supplies no validated human exposure, approved-product bridge, BUBR1 mediation or muscle result. Noncancerous-cell proliferation was slightly inhibited. | **Prioritize a new fidelity comparator for source and exposure qualification.** Neither reject for lacking the actual alleles nor promote directly to an existing-medicine nomination. [Primary abstract](https://pubmed.ncbi.nlm.nih.gov/30403302/). |
| **TG101209** | A renal-cancer study reports BUB1B thermal stabilization, with an opposing destabilizing derivative. | Thermal resistance is not cellular protein supply or functional checkpoint rescue; quantitative CETSA annotations have not been checked against the original figure. | **Protein-stability probe**, with abundance, turnover and fidelity tests before promotion. See [assay review](target-assay-review.md). |
| **Bletilla striata extract** | Reduced ovariectomy-associated rat bone loss and inhibited osteoclast formation; BUB1B appears in a thermal-proteomics/virtual-screening target network. | Genuine organismal benefit, but no isolated active ingredient, direction of BUBR1 action, causal mediation, muscle result or qualified drug exposure established in the accessible abstract. | **Retain a mixture-level reserve.** Do not assign extract efficacy to individual docked constituents. [Primary abstract](https://pubmed.ncbi.nlm.nih.gov/40220420/). |
| **PLK1 reduction** | Knockdown suppressed centrosome amplification in PCS/MVA patient cells with elevated Plk1 activity. | Strong same-disease cellular evidence, but genetic intervention. In mitosis PLK1 also supports BUBR1/PP2A-dependent attachment; generalized inhibition need not reproduce an interphase benefit. | **New target opportunity with a timing constraint**, not an elected drug. [PCS study](https://pubmed.ncbi.nlm.nih.gov/19503101/), [attachment mechanism](https://pubmed.ncbi.nlm.nih.gov/23079597/). |
| **GW5074** | RAF inhibition blocked pan-JAK-inhibitor-induced endoreduplication in HL-60 leukemia cells. | Pharmacological prevention of an imposed polyploidy phenotype; reduced BubR1 phosphorylation is not restoration of deficient BUBR1. No nonmalignant muscle or drug-eligibility bridge. | **Contextual pharmacology reserve**, below geraniin for the proposed fidelity comparison. [Primary abstract](https://pubmed.ncbi.nlm.nih.gov/20703093/). |

Geraniin's manuscript, exposure and original figure-level controls remain to be
qualified; its abstract is stronger than a docking hit but cannot supply missing
effect sizes or a therapeutic window. The accessible Bletilla and PLK1 abstracts
likewise do not establish all methods or drug implementation. These are concrete
follow-ups, not five newly proven medicines. Niacin and apitegromab also lack
actual-allele rescue; the comparison uses mechanism, administered-agent outcome,
model transfer, exposure/eligibility and harms consistently.

**Current nomination remains niacin restoration plus separate apitegromab muscle
compensation.** The new findings add research alternatives rather than establish
superior expected clinical efficacy. In particular, geraniin warrants checking
before declaring that the broad search has produced no further leads.

## Existing routes sharpened by the review

- **Aurora-B compensation:** the KARD study reports pharmacological restoration
  of kinetochore attachments after specified BUBR1 defects. This complements the
  already appraised ZM447439 experiment and conference AZD1152 mouse result;
  it is not a new patient-drug rescue or proof that all Aurora inhibition helps.
  PLK1 inhibition in the same abstract reduces BUBR1-associated PP2A recruitment.
  A separate abstract reports partial alignment rescue after **BUBR1 depletion**,
  but worsening after BUB1/BUB3 depletion; its exact inhibitor is unnamed in the
  abstract. This strengthens a context-specific compensation hypothesis without
  establishing accurate anaphase, muscle function or an eligible medicine.
  [KARD study](https://pubmed.ncbi.nlm.nih.gov/23079597/),
  [depletion study](https://pubmed.ncbi.nlm.nih.gov/18199686/).
- **Rapamycin:** the full published fly paper explicitly reports improved
  neuroblast counts after larval feeding, 40 µM for approximately 48 hours.
  Its genetic autophagy/TOR interventions did not restore progeny numbers or
  brain size. Do not transfer those unrescued endpoints indiscriminately between
  genetic and drug arms, or treat the journal abstract's omission of rapamycin
  as withdrawal of the result. The model concerns spindle-checkpoint-deficient
  fly neural tissue, not BUB1B muscle. [Results and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC13121444/).
- **NAD biology:** NMN-mediated protein abundance and SIRT2-transgenic survival
  are different interventions. Neither result directly establishes niacin rescue.
  The repeated abstracts do not multiply this evidence. [Primary study](https://pubmed.ncbi.nlm.nih.gov/24825348/).
- **Lithium/GSK3:** arrest release and lower checkpoint recruitment in one study
  counter any assumption that Axin1-oocyte morphology rescue generalizes to
  chromosome fidelity. [Primary study](https://pubmed.ncbi.nlm.nih.gov/30250048/).
- **Suan-Zao-Ren-Tang B7E2:** checkpoint amplification with vinorelbine is an
  administered chemosensitization finding. The outcome is cancer killing; no
  healthy-cell restoration follows from more BUBR1 phosphorylation or a longer
  arrest. The primary results use 30 µg/ml B7E2 plus 10 nM vinorelbine for the
  16-hour checkpoint assay. MRC-5 fetal lung fibroblasts were less affected in a
  separate cytotoxicity comparison, but normal-cell chromosome rescue was not
  established. A short acute B7E2-alone tolerance assay is not chronic safety.
  Individual ingredients and the mixture retain separate identities.
  [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12570107/).

## Source integrity and unresolved evidence

PMID 25262777, an Arabidopsis checkpoint paper, is explicitly marked as retracted
in its assessment through [the retraction notice, PMID 27874987](https://pubmed.ncbi.nlm.nih.gov/27874987/).
It must not support a therapeutic inference. The available 2020 JCI correction
changes Figure 2E tumor images and methods; it does not supply quantification
or total-protein denominators for the separately discussed mTOR blot.
[Correction](https://www.jci.org/articles/view/144781).

For the 28 records without saved abstracts, nine accessible short bodies were
read (corrections, commentaries and news); 19 remain at title/metadata depth.
The individual rows identify this boundary. Follow-up flags also retain drug
identity, dose, model, genotype and original-figure questions from read abstracts.
They are unresolved questions, not additional positive findings.

The earlier alias triage and 21-abstract sensitivity pass retain their historical
denominators. This deeper review does not make all 7,410 retrieved source records
fully reviewed, collapse versions into independent studies, or evaluate the
efficacy of every drug in the 15,357-entity catalog. World-wide completeness and
unpublished efficacy remain unestablished. The publication/submission condition
is therefore not asserted to have been met.
