# Next 128 saved-source revalidation

The additive review covers exactly the first 128 pending source/ID records in the original 3,598-record hash order after excluding the 186 accepted revalidations. `/root/revalidation_a64` and `/root/revalidation_b64` each read 64 full literal saved strings in 16 four-record batches. The 32 untruncated batches contain 195,082 abstract characters; every end marker was observed. All 128 saved content items are abstracts. These are source/ID records, not necessarily distinct studies or fully qualified papers.

Full evidence replay passed against both original retrieval files, the original selection and accepted-186 hashes, the earlier historical verifier, both unchanged reader input files, and every external batch log. The new review has 52 follow-up flags. The additive integration brings the outer ledger to 2,475 source/ID records and leaves 3,284 strings pending. The 314 accepted revalidations comprise the earlier 186 strings and these 128 abstracts; the older 2,161 annotations were not independently re-audited here. The 143 originally empty sources remain outside this queue and have separate recovery overlays. No historical ledger was rewritten.

## Evidence worth pursuing

| Source | What the saved abstract actually supports | Next primary-source question |
|---|---|---|
| MED:42305588 | Unnamed CDK1 inhibition attenuated fibrocystic liver disease in vivo in an ARPKD study. | Identify compound/exposure and whether the endpoint was cyst pathology, organ function, or another disease measure; retain the reported positive attenuation. |
| MED:41136336 | Increased sFRP3 expression preserved the neural stem-cell pool and improved spatial discrimination in AD-model mice. | Identify restoration modality and quantified stem-cell/behavioral outcomes; pharmacologic administration is not established. |
| MED:42054500 | ABBV-744/BET inhibition reduced plaque and inflammatory/foam-cell outcomes in an atherosclerosis model. | Determine model, dose, timing, effect sizes and adverse effects. |
| MED:41020815 | An unnamed COX-2 inhibitor partially reversed a SMIM25-associated HCC-cell phenotype. | Name the inhibitor and exact phenotype and distinguish cell testing from in vivo intervention. |
| MED:24561039; MED:17190794 | Active Rab11 improved spindle organization; exogenous Emi1 triggered meiotic entry. | Identify genetic/biological delivery and assess completed division and chromosome fidelity. |
| MED:31783707; MED:36675286 | SCI fibroblast proliferation changed with tiotropium/acetylcholine; Gja1 inhibition changed microglial number in opposite directions by coculture context. | Resolve compound-specific direction and perturbation identity; distinguish proliferation from survival. |
| MED:23719381; MED:34404896 | RAS–MAPK–MSK1 and MEK5/ERK5 pathways have pharmacologic evidence, with important functional benefits explicitly described for genetic perturbation. | Identify agents and assign each functional animal/combination result to the intervention actually used. |
| MED:29511940; MED:24748501 | rTM-treated DIC patients had 35.2% resolution; gastric-cancer methylation associated with docetaxel/oxaliplatin response and survival. | Qualify patient denominators, exposure/regimen and confounding; these observational results are not randomized or MVA rescue evidence. |

These benefits and quantitative outcomes are retained. None is upgraded to a direct MVA drug-rescue experiment by this batch, and current nominees niacin/apitegromab were not used to filter alternatives.

## Drug effects with adverse or cytotoxic division outcomes

CHR-6494 caused mitotic-catastrophe morphology, arrest and apoptosis along with antitumor effects (MED:21804608). PRT543/PRT811 and lenvatinib had anticancer activity (MED:39794830). BAY1816032 reduced hPheo1 tumor-cell viability/colonies/invasion (MED:42014906), while a separate source reports increased osteoclastogenesis under the same named Bub1 inhibitor (MED:38477771). Unnamed EZH1/2 inhibition impaired oocyte DNA repair/meiosis (MED:39511641). These directions are compatible with context-dependent effects and must not be summarized as generalized chromosome-fidelity restoration.

Additional explicit anti-proliferative/combination evidence includes BP-M345 analogs, tanshinone IIA, paclitaxel/212Pb-trastuzumab, PBOX-15/TRAIL, Ocoxin, BI2536/radiation, AK3/AK10/TNF and GSK-368/AZD-1390/AZD-6738. The annotations preserve named agents, assay/comparator roles, doses or ranges when supplied, and unreported identities. In particular, GSK-368 spelling and the nM EC50 ranges require source verification; Ocoxin and poplar propolis require formulation/mixture qualification. A low-cytotoxicity condition can still produce mitotic catastrophe and does not establish faithful productive division.

## Interpretation and source-quality questions

- Computational nominations are retained as such: cordycepin, CMap drugs, DGIdb interactions, expression/drug-sensitivity correlations and the 12 unnamed LVH candidates. The LVH abstract's literal “ruled out” wording leaves nomination versus exclusion ambiguous (MED:38994368). They are not discarded, but administration/efficacy is not inferred.
- PPR:PPR298223 is preprint version 1. Its reported 909 DEGs disagree with 453+457=910 and its hub/ROC lists differ. Remdesivir is a reused-dataset exposure, not a newly demonstrated candidate benefit. These discrepancies remain visible in the annotation.
- Three sources carry erratum links: MED:36781502 → MED:41484222; MED:34094924 → MED:38751807; MED:30553276 → MED:31783872. Correction texts were outside the fixed 128 and unread at that stage; the subsequent notice review below resolves what they changed. Source metadata has 17 records with comment/correction links: 6 “Comment in”, 3 “Erratum in” and 8 “Preprint in”; one record has a version list and is a preprint. Linked versions are not independent confirmation.
- The `Clinical Trial` indexing tag on MED:25502501 does not establish a human therapeutic trial or phase. Its observed lung-tissue species is unspecified in the abstract, while adoptive transfer to murine lung is explicit.
- Protein accumulation after Cdc20 silencing, gene-expression restoration after miR-34a/vGPCR knockdown, and engineered checkpoint arrest receive their reported mechanistic credit. Completed accurate division or absolute productive daughter-cell output is not supplied by those observations alone.
- `administered_agents` contains explicit role labels, including assay reagents, genetic constructs and predicted/background entities in some records. Counting nonempty arrays as administered therapeutic-drug evidence would be incorrect.

An independent 20-record literal-source sample audit identified three precision corrections. The [consolidated review](outer-revalidation-next128-review.json) records each old/new field with source hash, reason, readers and an additional literal reread log: lung-tissue species in MED:25502501; unspecified disease endpoint in MED:42305588; mitotic-catastrophe morphology versus completed-division fidelity in MED:21804608. The original part files remain byte-identical. The extra audit readings add zero new coverage credit.

## Linked corrections: subsequent notice review

All three notices were then read, separately from the fixed abstract batch. The
[structured follow-up](next128-errata-followup.json) records source identities,
reading depth and the effect on interpretation; it adds no original-abstract
coverage credit.

| Notice | What changed | Consequence for this search |
|---|---|---|
| [PMID 41484222](https://pmc.ncbi.nlm.nih.gov/articles/PMC12764733/) | Replaced a duplicated Nalm6 transwell image in Figure 5i of the MAD2L1 paper | Use the corrected figure. Its caption also identifies actual **sc-43 and garcinone D** drug additions in different engineered leukemia-cell backgrounds. Their identities, doses and selectivity require primary-methods follow-up; neither was named in the saved abstract. |
| [PMID 38751807](https://pmc.ncbi.nlm.nih.gov/articles/PMC11094535/) | Corrected an author affiliation in the pancreatic-cancer microtubule review | No scientific text or figure change in this notice; agent-level literature review remains open. |
| [PMID 31783872](https://pmc.ncbi.nlm.nih.gov/articles/PMC6883572/) | Replaced a duplicated Figure 4e apoptosis plot and repaired the c-Src paper's title | The corrected experiment uses genetic DMAP1/Bub3/Bcl-xL manipulations under nocodazole stress. It does not establish that an administered c-Src inhibitor produced the separate paclitaxel result. |

The authors state that their conclusions are unchanged. That statement is
retained as their claim, not independent validation of the underlying data. The
notices do not require reversal of these saved-abstract appraisals, but the two
figure corrections affect which evidence should support subsequent comparisons.
The new sc-43/garcinone D pointers remain research questions, not added efficacy
appraisals or independent replication.

## Verification boundary

`verify_revalidation_next128.py` calls the historical verifier, checks the exact first 128 after the accepted-186 exclusion, preserves every selection field, requires complete per-record reading declarations, matches both input annotations and the explicit revision ledger, validates external literal logs/end markers, and derives counts. With both corpora supplied, it replays the original full queue and matches every literal string and publication flag to the saved source. Nine meaningful corrupted fixtures were rejected, including absent/truncated logs, altered queue order, empty readers, untraced annotations, broken correction provenance, false metadata counts, changed accepted-186 bytes and wrong source hashes.

The full verifier requires external evidence and the historical workflow directory. It should be run separately from ordinary repository structural checks; it does not silently make those checks depend on private source logs. A pass establishes bookkeeping/provenance and declared-depth consistency, not scientific reading fidelity, full-paper qualification, patient phase, missense effect or clinical applicability. These findings provide research-source questions, not clinical recommendations.

Run the full evidence replay with the original retrieval files and saved external
reading evidence (the two part assignments/reviews and batch logs):

```sh
python3 workflow/track2-drug-census/verify_revalidation_next128.py \
  --evidence-root /path/to/revalidation-next128 \
  --corpus /path/to/literature-gene-unrestricted.json \
  --gene-corpus /path/to/literature-gene-title-abstract.json
```

The repository-only `verify_census.py` checks the new ledger's selection, source
identities and totals without claiming to replay external reading evidence.
The separate [three-paper qualification](vascular-enzyme-qualification.md) adds
six compound/tool appraisals, bringing the ledger to 122; it adds no abstract
coverage credit.
