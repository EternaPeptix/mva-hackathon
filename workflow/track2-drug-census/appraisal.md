# Target-validation appraisal: BUB1 versus BUB1B, and BUBR1 rescue candidates

Primary-source check completed 2026-09-21. This is a bounded literature appraisal,
not a systematic review or treatment recommendation. Concentrations below identify
published experiments; they are not proposed clinical doses.

**Decision:** do not classify BAY-320, BAY 1816032 or 2OH-BNPP1 as validated direct
BUB1B/BUBR1 inhibitors on the basis of the 2026 hepatocellular-carcinoma paper.
Their original target-validation literature concerns **BUB1**, a separate checkpoint
protein. A genuine direct-binding, rescue-related BUBR1 research candidate was
identified: **UPCP**, an experimental peptide studied in progerin-driven HGPS.
It has not been shown to restore a pathogenic BUB1B allele or MVA tissue function.
No selective direct small-molecule BUBR1 stabilizer/activator with MVA rescue was
verified in this search beyond the separately considered NAD/SIRT2 and chaperone
hypotheses. That is a search outcome, not proof that no such compound exists.

## 1. What the 2026 HCC paper actually establishes

Sun et al. call the three compounds selective BUB1B inhibitors. In HepG2 cells,
they report RAS/MAPK **mRNA** changes after 24 hours of BAY-320 or 2OH-BNPP1
at 10 µM, or BAY 1816032 at 1.4 µM. These transcriptional responses do not
themselves demonstrate BUBR1 binding, selective inhibition, restored checkpoint
function or pathway phosphorylation. Their clinical-expression associations
are a separate result; their drug-sensitivity analysis uses transcriptomic
pRRophetic/GDSC predictions, not a prospective treatment-response trial.
The publisher's June 2026 update concerns equal author contributions and does
not resolve the target assignment. The inhibitor-based attribution to BUB1B
should therefore be withheld pending independent target engagement and genetic
specificity experiments. This criticism does not, by itself, refute the paper's
expression/prognosis associations. [Sun et al., Scientific Reports 2026,
doi:10.1038/s41598-026-36364-x](https://www.nature.com/articles/s41598-026-36364-x).

## 2. Original chemical validation and the wild-type/engineered distinction

| Compound | Measured primary evidence | Appropriate census classification |
|---|---|---|
| **BAY-320** | Baron et al. assayed human wild-type LAP-Bub1 purified from mitotic HEK293T cells against histone H2A, with kinase-dead Bub1 as a control. The biochemical IC50 was **680 nM at 2 mM ATP** (Supplementary file 1). Cellular H2A-T120 phosphorylation and Shugoshin/CPC localization were affected in HeLa/RPE1 cells. Paclitaxel combinations impaired chromosome segregation/proliferation. | **BUB1 catalytic inhibitor**, with endogenous cellular pathway evidence. The engineered tag and kinase-dead controls do not make it a BUBR1 inhibitor. [Baron et al., eLife 2016](https://doi.org/10.7554/eLife.12187); [original supplementary potency table](https://doi.org/10.7554/eLife.12187.018). |
| **BAY 1816032** | Direct BUB1 binding and a BUB1–compound crystal structure support the target. Reported BUB1 binding Kd **3.3 nM** contrasted with **57 nM** for the leading kinase off-target LOK/STK10. In nocodazole-arrested HeLa cells, H2A-T120 phosphorylation was inhibited after one hour with IC50 **29 ± 23 nM**. Across 43 tumor cell lines, median growth-inhibition IC50 was **1.4 µM**. | **BUB1 inhibitor**. The low-nanomolar target readout and micromolar growth readout are different endpoints. Neither validates inhibition of BUB1B. [Siemeister et al., Clinical Cancer Research 2019, PMID 30429199](https://aacrjournals.org/clincancerres/article/25/4/1404/10115/Inhibition-of-BUB1-Kinase-by-BAY-1816032). |
| **2OH-BNPP1** | Kang et al. tested full-length **wild-type human Bub1** and its C-terminal kinase fragment; biochemical IC50 was approximately **250 nM**. The compound came from a class designed for analog-sensitive engineered kinases, but Bub1 naturally has the small glycine gatekeeper that permits binding. The original cellular phenotype statement was explicitly unsupported by displayed data. | **Biochemical WT BUB1 inhibitor with limited cellular specificity**, not an engineered-BUBR1-only inhibitor and not validated WT BUBR1 inhibition. Do not confuse the origin of the chemical class with the genotype of the target assayed. [Kang et al., Molecular Cell 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2644263/). |

Two later experiments constrain use of those annotations. Amalina et al. found
that **10 µM BAY-320**, but not **10 µM 2OH-BNPP1**, efficiently inhibited an
engineered H2B–Bub1-kinase-domain cellular reporter; only BAY-320 effectively
displaced endogenous Bub1-dependent Shugoshin. This engineered reporter is
**Bub1**, not BubR1; the paper also warns about BAY-320 off-target effects at
the required concentrations. It does not prove that 2OH-BNPP1 can never inhibit
Bub1 in any cell context. [Amalina et al., Royal Society Open Science 2021,
PMID 34925867](https://doi.org/10.1098/rsos.210854).

The BAY 1816032 study additionally found 2OH-BNPP1 activity against PDGFRβ,
CSF1R, VEGFR2 and VEGFR3 at **30–100 nM**. Its TGFβ/SMAD effects were not
reproduced by the more selective BUB1 compounds. This creates an additional
off-target explanation for downstream transcriptional changes, independently
of the BUB1/BUB1B name error. [Siemeister et al., supplementary Table S3 and
cellular characterization](https://aacrjournals.org/clincancerres/article/25/4/1404/10115/Inhibition-of-BUB1-Kinase-by-BAY-1816032).

**Control against overcorrection:** a separate compound, **bubristatin-1
(BRT-1)**, does have published direct BubR1-antagonist evidence. Huang et al.
reported biochemical IC50 **1.2 µM**, affinity capture, binding-site mutants
that retained activity but resisted inhibition, and inhibition of CENP-E
phosphorylation with **3 µM** compound. It acts on wild-type BubR1; engineered
resistance mutants test specificity rather than create sensitivity. This
research antagonist impaired attachment/spindle behavior; it is not a
stabilizer, activator or MVA rescue agent. [Huang et al., Cell Research 2019](https://www.nature.com/articles/s41422-019-0178-z).

## 3. Rescue-related findings beyond the existing NAD/SIRT2 and chaperone hypotheses

### UPCP: direct binding, disease-specific rescue evidence

Zhang et al. designed a cell-permeable progerin-derived peptide, **UPCP**.
Pulldown and microscale thermophoresis showed direct BUBR1 binding. In HGPS
and progerin-expression models it disrupted progerin–BUBR1 sequestration,
improved localization/checkpoint response and increased BUBR1 abundance.
The abundance mechanism included **PTBP1-dependent BUBR1 mRNA stabilization**;
it should not be described simply as stabilization of mutant BUBR1 protein.
The concentration series for BUBR1 abundance was **0–10 µM**; interaction
competition experiments used **6, 12 and 18 µM**. In LmnaG609G/G609G mice,
intraperitoneal **20 mg/kg every other day**, beginning at 14 days of age,
was associated with improved progeroid phenotypes and longer survival
(15 treated versus 17 vehicle mutant mice in the survival analysis).
These are HGPS, not BUB1B-mutant MVA, experiments. Transfer requires evidence
that the relevant sequestration/transcript-loss mechanism exists in the
target cells; the peptide cannot be assumed to restore a truncated protein.
Classify as **direct-binding research comparator with HGPS rescue**, not a
validated MVA drug or general catalytic activator. [Zhang et al., Nature Aging
2023, PMID 37118121](https://www.nature.com/articles/s43587-023-00361-w).

**Correction checked:** Figure 7 originally mislabeled esophagus and lamina
propria as aorta and aortic adventitia. The corrected paper must be used;
those images do not establish aortic rescue. The correction does not retract
the separately reported binding or survival experiments.
[Author correction, 2 May 2023](https://www.nature.com/articles/s43587-023-00427-9).

### Trichostatin A: real abundance effect through HDAC inhibition

In synchronized HeLa cells, **100 nM trichostatin A (TSA)** maintained higher
BubR1 protein levels during release from nocodazole arrest, including with
cycloheximide to suppress new protein synthesis. The study linked K250
acetylation to reduced ubiquitination/proteolysis. This is a measured
pharmacological stabilization effect, but TSA acts on deacetylases; it is
not demonstrated to bind and selectively stabilize BUBR1. The experimental
endpoint was mitotic regulation, not MVA tissue rescue. [Choi et al., EMBO
Journal 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2684026/).

A later primary study showed that sustained acetylation can keep cells in
metaphase even after chromosome congression, and that HDAC2/3-mediated
deacetylation participates in checkpoint silencing. Thus more persistent
BubR1 is not automatically better division or better tissue output. This
K250 mitotic mechanism must also remain distinct from SIRT2-linked
stabilization at other sites. Classify TSA as an **indirect mechanistic
stabilization probe**, not an established substitute for niacin or the
everolimus probe. [Park et al., FEBS Journal 2017, PMID 28985013](https://pubmed.ncbi.nlm.nih.gov/28985013/).

### Indibulin: checkpoint activation secondary to spindle damage

In MCF-7 cells treated for 24 hours, the fraction of mitotic cells positive
for kinetochore BubR1 increased from **5% without compound** to **38%, 74%
and 95% at 150, 300 and 600 nM indibulin**, respectively. The experiment
also documented aberrant spindles, reduced kinetochore tension, mitotic
arrest and subsequent apoptosis. These results explain commercial
“BUBR1 activator” labels without establishing direct BUBR1 binding or repair
of deficient BUBR1. Classify as a **microtubule-directed antimitotic with
secondary SAC activation**, not a rescue/stabilization lead. [Kapoor et al.,
Scientific Reports 2018, PMID 30120268](https://www.nature.com/articles/s41598-018-30376-y).

### CENP-E and Aurora B: mechanistic comparators, not direct stabilizing drugs

Purified CENP-E stimulated recombinant Xenopus BubR1 autophosphorylation;
human-cell work connected the corresponding phosphorylation site to
chromosome alignment/checkpoint behavior. This is a protein-partner
activation experiment, not evidence for a deliverable small-molecule
BUBR1 agonist. [Guo et al., Journal of Cell Biology 2012](https://rupress.org/jcb/article/198/2/205/36875/CENP-E-dependent-BubR1-autophosphorylation).

Another study showed that disrupting BUBR1–B56/PP2A binding impairs chromosome
congression, including mechanistic findings in MVA-derived cells; Aurora B
inhibition suppressed the attachment/congression defect in the experimental
system. Suppressing a downstream phenotype does not establish BUBR1 protein
stabilization, recovery of every checkpoint function or durable developmental
benefit. Retain such compounds as **pathway-bypass controls**, with the tested
cell/genetic context explicit. [Xu et al., Biology Open 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3654266/).

## 4. Consequences for the drug census

1. Separate **BUB1** and **BUB1B/BUBR1** target fields. Record changes in BUBR1
   abundance, localization, phosphorylation and cellular dependence as distinct
   observations, not interchangeable evidence of direct binding.
2. Preserve UPCP as a substantive experimental comparator, TSA as indirect
   stabilization evidence, and BRT-1 as a direct antagonist/control. Exclude
   generic SAC activation by spindle poisons from a list of demonstrated
   BUBR1-restoring agents.
3. None of these results independently establishes niacin efficacy or justifies
   switching the existing lead/probe ordering. A proposed upstream rescue must
   demonstrate functional allele product, appropriate checkpoint timing,
   segregation fidelity and tissue output in relevant cells. A downstream
   intervention must demonstrate its own output benefit without calling a
   selected surviving fraction genetic correction.
4. The search covered direct BUBR1 binding/stabilization/activation, the three
   named inhibitors' original studies, chemical-genetic controls, acetylation,
   UPCP/progerin and commercial activator claims. No docking-only result or
   expression/drug-sensitivity correlation was upgraded to target validation.
   These classifications are experimentally bounded; they make no individual
   treatment or pediatric-safety claim.

## 5. Follow-up: dasatinib and ceritinib biochemical activity records

The retained [clinical-drug activity rows](bub1b-clinical-drug-activities.tsv)
contain two clinical drugs. An activity record means that a compound was
assayed; it does not establish a useful interaction or a rescue indication.

**Dasatinib:** CHEMBL4806967 is **patent WO2020083926A1**, not a journal
article. Example 3 describes duplicate, single-point **1 µM** testing across
**320 wild-type kinases** with the ProQinase Kinase Profiler. Figure 3 defines
its weakest inhibition category as **greater than 75% residual activity**.
The retained BUB1B row, assay CHEMBL4807760/activity 23127857, falls in that
category: **less than 25% inhibition**, not greater than 75% inhibition and
not evidence of activation. It supplies neither a BUB1B IC50 nor evidence of
cellular checkpoint rescue. Classify as **low-inhibition kinome-screen
observation at the tested concentration**, preserving the patent provenance.
The assay scope and category definition were checked in the primary patent;
the BUB1B-specific category is retained from ChEMBL's curation rather than
independently re-read from the figure image.
[Primary patent, Figure 3 and paragraphs 513–514](https://patents.google.com/patent/WO2020083926A1/en).

**Ceritinib:** CHEMBL6121040 maps to El Hafi et al., *European Journal of
Medicinal Chemistry* **285:117247 (2025)**, PMID 39818011. The retained
assay CHEMBL6125423/activity 29222147 reports recombinant BUB1B inhibition
with **IC50 13.08 µM (13,080 nM)**, using a **60-minute [γ-33P]-ATP
reaction** and scintillation readout. These numerical/method details are
**ChEMBL-curated and were not independently checked against the full article's
figure or methods**, which were unavailable through the accessible publisher
and author-repository routes. The primary abstract was checked: the study
concerns 40 TG-101209 derivatives, Caki-1 cancer-cell viability/apoptosis,
enzyme inhibition and thermal-shift experiments. It is not an MVA rescue
study. Neither this single-protein result nor the abstract establishes
ceritinib's selectivity across the kinome or a BUBR1-mediated cellular benefit.
Classify as **reported weak micromolar biochemical inhibition; no demonstrated
rescue**, not a direct stabilizer/activator or a clinically validated BUB1B
target. Thermal stabilization reported for TG-101209 in that paper is a
binding-assay observation, not proof of restored checkpoint function, and
must not be transferred to ceritinib.
[Primary paper, doi:10.1016/j.ejmech.2025.117247](https://doi.org/10.1016/j.ejmech.2025.117247);
[indexed primary abstract](https://pubmed.ncbi.nlm.nih.gov/39818011/);
[ChEMBL assay record](https://www.ebi.ac.uk/chembl/explore/assay/CHEMBL6125423).

Neither record upgrades cancer-directed inhibition to BUB1B loss-of-function
rescue. The census should retain them with the measured direction and evidence
limits, without changing the lead/probe ordering or inferring clinical safety
from the compounds' clinical-drug status.

## 6. Newly reported direct binding: griseofulvin is not a rescue result

Pham et al., *Bioorganic Chemistry* **180:110130 (2026)**, PMID 42330849,
report that **griseofulvin binds directly to Mad1 and BubR1**, disrupts the
Mad1–Mad2 interaction and promotes multipolar spindles and mitotic cell death.
The distinct biosynthetic shunt product **norlichexanthone** is reported to
bind **TPR**, impairing checkpoint assembly and nucleocytoplasmic transport;
its target and effects must not be assigned to griseofulvin. The primary
abstract and publisher highlights support this direction. Full experimental
methods/figures were not accessible here: binding constants, compound
concentrations, cell models and orthogonal specificity controls remain
unverified in this appraisal. Preserve the result as an **author-reported
direct-binding mechanism with antimitotic effects**, not an independently
validated selective BUBR1 agonist or mutant rescue.
[Primary paper, doi:10.1016/j.bioorg.2026.110130](https://doi.org/10.1016/j.bioorg.2026.110130);
[primary abstract](https://pubmed.ncbi.nlm.nih.gov/42330849/).

Griseofulvin has an existing approved oral antifungal product, as shown by
the label's ANDA091592 record. That status does not establish a clinically
useful BUBR1 interaction or MVA safety/efficacy.
[Current product label](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=98829d90-b289-4bec-a625-7af412a5cb2d).
This is a substantive literature finding missing from the retrieved curated
BUB1B mechanism rows, so their zero count must remain explicitly
database-specific. It supplies no demonstrated muscle preservation or
functional restoration and does not itself justify lead promotion.
