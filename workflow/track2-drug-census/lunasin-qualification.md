# Lunasin: a BUB1B transcript lead with an unclosed protein/function bridge

**28 September 2026. Retain as a food-derived restoration reserve.** Lunasin
has an actual peptide-treatment BUB1B RNA signal, so it merits more than
dismissal as an incidental food mention. A processed GEO matrix reanalysis
reproduces the reported RWPE-1 fold change and adds the RWPE-2 comparison, but
uses one array per condition and is not independent biological replication.
It does not displace niacin's restoration hypothesis or apitegromab's
muscle-compensation rationale. None has demonstrated MVA rescue; this is a
research-priority judgment, not measured clinical superiority.

## The direct BUB1B result

[Galvez et al., 2011](https://doi.org/10.1080/01635581.2011.539312)
(PMID 21526452, PMC3210036), Table 2, reports **2.23-fold BUB1B RNA
upregulation** after synthetic lunasin, 2 µM for 24 hours, versus PBS in
non-tumorigenic RWPE-1 prostate epithelial cells. TTK and MAD2L1 also appear in
the increased-transcript list. The corresponding tumorigenic RWPE-2 cells did
not meet the paper's greater-than-twofold broad-response criterion.

The [GEO GSE2992 series matrix](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2992)
reports processed values as `2^(RMA-normalized data)`. Its four samples contain
one control and one treated array for each line. The GPL96 annotation maps probe
203755_at to BUB1B. The reanalysis gives an RWPE-1 treated/control ratio of
**2.229** (rounds to the paper's 2.23) and an RWPE-2 ratio of **1.252**.
This is consistent with the published table, not a replication: the GEO page
explicitly says **“Citation missing.”** Its title/topic, cell lines, treatment,
authors and RWPE-1 value support an accession-to-paper match, but that linkage
is inferred rather than cited by GEO. The values are from the same likely
source arrays, with no replicate-level variance; no inferential test is
possible. The [reanalysis output](lunasin-gse2992-reanalysis.json),
[source manifest](lunasin-gse2992-source-manifest.json) and
[replay script](reanalyze_lunasin_gse2992.py) preserve the checksums and mapping.

This was a microarray finding. A BUB1B-specific confidence interval, adjusted
P value and independent PCR confirmation are not supplied. The four genes
confirmed by PCR and protein experiments were HIF1A, PRKAR1A, TOB1 and THBS1;
**BUB1B was not among them**. Histone-acetylation and promoter experiments
support a THBS1/chromatin mechanism, not a demonstrated BUB1B-specific one.
No BUBR1 protein, mutant allele, chromosome-fidelity or muscle-function rescue
was measured in this study. RWPE-1 is a cultured non-tumorigenic line, not an
MVA patient-derived normal tissue.

The full main text, both tables and all seven figure captions were read;
Figures 1 and 6 were visually inspected. The deposited processed series matrix
was reanalyzed; raw CEL files were not independently reprocessed, and no new
biological experiment was performed. A separate source check agreed with the
narrow transcript-level conclusion; it was not a blinded review of the whole
proposal.

## Human exposure and counterevidence

| Source and reading depth | What it contributes | Limit for MVA |
|---|---|---|
| [2020 plasma conference abstract](https://pmc.ncbi.nlm.nih.gov/articles/PMC7258671/), read completely | In a 31-person crossover study, an enriched soy extract supplied reported lunasin exposure of 50 mg/day for eight weeks. Fasting plasma ELISA values were 0–10 ng/mL, with higher mean immunoreactivity during supplementation. | Conference-level report, not full pharmacokinetic qualification. It does not establish intact active peptide or intracellular muscle exposure comparable with the culture experiment. |
| [2019 cardiometabolic trial](https://pubmed.ncbi.nlm.nih.gov/30955400/), complete primary abstract only | The 31-person crossover trial reported no significant improvements in its measured cardiometabolic endpoints. | The matching design/demographics may represent the same cohort as the plasma report; do not count independent replication. Different endpoints from MVA. |
| [2019 ALS trial](https://pubmed.ncbi.nlm.nih.gov/30663902/), complete primary abstract only | Fifty participants, open-label with historical controls: no significant histone-acetylation or progression benefit. | Weakens a generic neuromuscular-benefit claim, but does not establish MVA inefficacy. Full trial methods and product details were not qualified here. |
| [1999 transfection study](https://doi.org/10.1038/8676), publisher abstract only | Intracellular expression of GFP-tagged lunasin disrupted spindle/chromosome behavior and caused cell death. | A genetic-expression experiment, not proof that oral or externally supplied peptide causes those effects. It supplies a reason to test faithful division directly. |

As a scale check only, 10 ng/mL expressed as intact 5.5-kDa peptide would be
about 1.8 nM, roughly 1,100-fold below the 2 µM culture concentration
([purified lunasin mass characterization](https://pmc.ncbi.nlm.nih.gov/articles/PMC3326064/)).
This is not a pharmacokinetic comparison or an exposure threshold: the plasma
ELISA does not establish that the measured immunoreactivity is intact active
lunasin, and plasma and culture-medium concentrations are different matrices.

## Why it stays a reserve

The positive RNA direction makes lunasin a candidate for testing whether more
**usable BUBR1** can be produced. The missing steps are measured protein
abundance/function, an exposure bridge, and faithful productive divisions in a
deficient model. These are not waived for niacin either. Apitegromab has a
different rationale supported by functional evidence in another muscle disease;
lunasin's source does not supply comparable muscle output.

Synthetic peptide, enriched extract and dietary soy are different interventions.
Neither the culture result nor ELISA detection establishes that eating soy
reproduces the experiment, or that a consumer supplement has a useful MVA
exposure. Product composition, intact peptide delivery and batch identity remain
unqualified. A future result showing only higher RNA or faster growth would
not settle the restoration question.

This source was recovered from the SGC screen's false **BSP** match: BSP means
bisulfite PCR in this paper, not bromosporine. The false match is preserved;
lunasin receives its own appraisal. [Structured evidence](lunasin-qualification.json)
records hashes and reading limits. No extra outer-corpus completion credit is added.
