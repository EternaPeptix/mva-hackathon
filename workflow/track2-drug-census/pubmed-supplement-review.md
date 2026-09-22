# Supplemental PubMed version and mechanism appraisal

Reviewed **21 September 2026**. This is a research evidence appraisal, not a treatment recommendation. It examines three mechanistically relevant preprint records in the supplemental PubMed retrieval. It does not classify all supplemental records or add experimental results.

## Version reconciliation

All three have peer-reviewed journal versions. Their separate PubMed identifiers must not be counted as independent studies or independent rescue evidence.

| Supplemental preprint | Journal version used for appraisal | Existing census status |
|---|---|---|
| [PMID 42619900](https://pubmed.ncbi.nlm.nih.gov/42619900/), DOI `10.64898/2026.07.30.741836` | [PMID 42647139](https://pubmed.ncbi.nlm.nih.gov/42647139/), *Cancer Research Communications* 2026;6:2110–2125; DOI `10.1158/2767-9764.CRC-26-0560`; [PMC13555252](https://pmc.ncbi.nlm.nih.gov/articles/PMC13555252/) | Journal record already present in the broad corpus; an abstract lacking a catalog drug-name match still warrants mechanism review. |
| [PMID 38410481](https://pubmed.ncbi.nlm.nih.gov/38410481/), DOI `10.1101/2024.02.15.580539` | [PMID 38987356](https://pubmed.ncbi.nlm.nih.gov/38987356/), *Scientific Reports* 2024;14:15912; DOI `10.1038/s41598-024-66244-1`; [PMC11236977](https://pmc.ncbi.nlm.nih.gov/articles/PMC11236977/) | Journal record already represented by paclitaxel and BAY-1217389 rows in `literature-mentions.tsv` and `literature-triage.tsv`. |
| [PMID 42539170](https://pubmed.ncbi.nlm.nih.gov/42539170/), DOI `10.64898/2026.07.23.740390` | [PMID 42667124](https://pubmed.ncbi.nlm.nih.gov/42667124/), *FASEB Journal* 2026;40:e72246; DOI `10.1096/fj.202600627RR`; [PMC13525491](https://pmc.ncbi.nlm.nih.gov/articles/PMC13525491/) | Journal record already represented by the chorionic-gonadotropin background row. That hormone mention is not the experimental intervention. |

PubMed explicitly links the KIF18A and trophoblast preprints to their journal updates. The pancreatic study has the same authors and study, with the journal title ending “in a pancreatic cancer cell line.” The journal versions control the appraisal below. Full journal text was inspected through Europe PMC XML for KIF18A and pancreatic cancer and through the publisher for trophoblast development. No raw abstracts or full texts were added to the repository.

## KIF18A inhibition: checkpoint-dependent tumor killing

Tabe et al. tested the preclinical KIF18A inhibitor **AM-9022** in SCLC cells. Resistant models had impaired MAD1/BUBR1 kinetochore recruitment. Concurrent **reversine** treatment weakened checkpoint arrest and allowed abnormal chromosome segregation. A different sequence—reversine pretreatment, washout/recovery, then AM-9022—partially restored killing in a resistant model. The authors explicitly distinguish this schedule-dependent escalation of chromosomal instability from generic sensitization by checkpoint inhibition. The journal paper reports 25 screened SCLC lines and describes the mechanistic findings as in vitro; AM-9022 pharmacology and off-target effects remain incompletely defined. [Journal full text, Methods, Figures 4–6 and Discussion](https://pmc.ncbi.nlm.nih.gov/articles/PMC13555252/).

**Classification:** direct pharmacology in cancer, with BUBR1 functional involvement; no BUB1B-MVA allele or muscle-rescue experiment. The inference relevant to this project is that baseline aneuploidy alone cannot establish a therapeutic vulnerability. Restored tumor killing after deliberately increasing segregation errors is not restored chromosome fidelity. Neither KIF18A inhibition nor reversine earns an MVA-restoration lead from this study.

## Nab-paclitaxel: same-gene perturbation, adverse rescue meaning

Mondal et al. combined CRISPR discovery with stable **BUB1B, BUB3 and TTK shRNA knockdowns** in PANC-1 cancer cells. Knockdown increased nab-paclitaxel resistance by about tenfold and was associated with chromosome abnormalities. Direct **BAY-1217389** and **MPI-0479605** exposure antagonized nab-paclitaxel cytotoxicity. By contrast, **BAY-1816032**, targeting **BUB1 rather than BUB1B**, increased cytotoxicity with nab-paclitaxel. Overexpressing the three SAC genes did not significantly change drug sensitivity. [Journal full text, Figures 2–5](https://pmc.ncbi.nlm.nih.gov/articles/PMC11236977/).

**Classification:** same-gene genetic context plus direct pharmacology against other SAC components, in one tumor cell-line setting. “Rescue” here means cancer-cell survival during cytotoxic treatment, accompanied by compromised segregation control. It cannot establish rescue of healthy BUB1B-deficient muscle, general safety of SAC suppression, or an MVA-specific chemotherapy-response rule. The differing BUB1 result also prevents treating all checkpoint targets as interchangeable.

## Trophoblast kinases: developmental disruption, qualified target attribution

Kumar et al. exposed human CT27 trophoblast stem cells to **Antitumor agent 180 (ATA180), 250 nM**, described by the authors as a BUB1B inhibitor. Cells showed reduced proliferation/Ki67 and stem-state markers, flattened colonies, and induction of the syncytiotrophoblast marker CGB. Their functional validation used reduced histone-H3 Ser10 phosphorylation, an indirect mitotic/Aurora-B-linked readout. This is not direct proof of selective BUB1B binding or an allele-specific rescue. **GNE2861** in this study targets group-II PAK activity including PAK6; it is not evidence for pharmacological PAK2 rescue. [Journal full text, Table 1, section 3.5 and Figure 5](https://faseb.onlinelibrary.wiley.com/doi/10.1096/fj.202600627RR).

**Classification:** direct pharmacological developmental perturbation with qualified BUB1B attribution. It supplies an adverse-direction signal in a placental cell model, not measured pregnancy outcomes or muscle efficacy. Chorionic gonadotropin is a trophoblast product/marker, not a BUB1B-rescuing treatment in this experiment.

## Decision consequence

These records improve **version accounting and source-role interpretation**. They do not add three independent studies, three therapeutic leads, or evidence that niacin, apitegromab or everolimus works in BUB1B-MVA. They do not establish a reason to change the current niacin-restoration/apitegromab-compensation research comparison or promote conditional everolimus.

The methodological implication is to retain **checkpoint competence and chromosome-segregation fidelity** as distinct outcomes alongside viability and muscle function. Longer arrest, more surviving cells, or more tumor killing can each have different biological meanings. Drug identity, exact target, tissue, treatment sequence and functional direction must be checked before a lexical match becomes a candidate-efficacy claim.

## Remaining specific aliases and fourth version link

The full supplement contains 74 records, of which 48 PMIDs are absent from the
base corpus. Four preprint records have explicit XML `UpdateIn` links to journal
records already present, including [SPC105R PMID 38559067](https://pubmed.ncbi.nlm.nih.gov/38559067/)
→ [PMID 38865189](https://pubmed.ncbi.nlm.nih.gov/38865189/). Its saved abstract
describes a genetic domain experiment in fly oocytes, not an administered rescue
drug; only the three papers above received full-text appraisal in this pass.

All six specific-alias rows now have an interpreted role. Paclitaxel, BAY-1217389
and chorionic gonadotropin are covered above. Complete saved abstracts were also
read for the other three records: [iodine PMID 22162022](https://pubmed.ncbi.nlm.nih.gov/22162022/)
is a radioactive environmental exposure in an observational trisomy paper,
[protamine PMID 24369136](https://pubmed.ncbi.nlm.nih.gov/24369136/) is endogenous
sperm chromatin biology in a review, and [tamoxifen PMID 22884371](https://pubmed.ncbi.nlm.nih.gov/22884371/)
induces an engineered recombination system. These are not administered BUB1B
restoration treatments. This role classification does not independently validate
the observational paper's causal interpretation. Thirty ambiguous aliases remain
retrieval metadata; no therapeutic conclusion is assigned to them.
