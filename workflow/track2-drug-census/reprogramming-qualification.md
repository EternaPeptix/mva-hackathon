# A2S reprogramming: compound identity and outcome qualification

[Tran et al., 2019](https://doi.org/10.1016/j.celrep.2019.04.056) ([PMID 31067459](https://pubmed.ncbi.nlm.nih.gov/31067459/); [PMC6555151](https://pmc.ncbi.nlm.nih.gov/articles/PMC6555151/)) report genuine pluripotency outcomes from chemically assisted reprogramming of embryonic mouse fibroblasts. Preserve the expandable, withdrawal-resistant colonies and reported normal karyotype and three-germ-layer teratomas. These results belong to the tested genetic/cocktail context; they do not establish BUBR1 rescue or muscle treatment.

## Regimen and identities

Male and female E13.5 fibroblasts carried an inducible Oct4/Klf4/Sox2/c-Myc cassette at Col1a1 and rtTA at Rosa26. Doxycycline, 2 micrograms/mL from day 0, induced those factors. Cultures used ESC medium containing 15% FBS and LIF, irradiated MEF feeders and gelatin-coated coverslips. This was not chemical-only reprogramming.

| Component | Primary dose and start | Study catalog and identity limit |
|---|---|---|
| Paper's “ascorbic acid” (AA) | 50 micrograms/mL, day 0 | Sigma A8960. The current exact catalog identifies **L-ascorbic acid 2-phosphate sesquimagnesium salt hydrate**, also called magnesium ascorbyl phosphate. Historical lot/COA was not verified; do not silently substitute plain ascorbic acid or oral vitamin C. |
| SGC0946 | 5 micromolar, day 0 | APExBIO A4167; DOT1L inhibitor, chemically distinct from EPZ-5676/pinometostat and EPZ004777. |
| CHIR-99021 | 3 micromolar, 12 hours after doxycycline | Stemgent 04-0004-10; GSK3 inhibitor, one component of 2i. |
| PD-0325901 | 1 micromolar, 12 hours after doxycycline | Stemgent 04-0006-10; MEK inhibitor, the other 2i component; existing `new-mirdametinib` identity. |

Medium containing doxycycline and chemicals was replaced every two days. Colony assays were read at the specified day or four additional days after withdrawal of **both** doxycycline and the chemicals. Solvents, final vehicle fractions and vehicle matching are not specified in the accessible Methods; supplier solubility guidance cannot fill these gaps.

Current manufacturer listings confirm research catalog identities: [A8960](https://www.sigmaaldrich.com/US/en/product/sigma/a8960), [SGC0946](https://www.apexbt.com/sgc-0946.html), [CHIR99021](https://www.reprocell.com/product-catalog/small-molecules/stemolecule-chir99021) and [PD0325901](https://www.reprocell.com/product-catalog/small-molecules/stemolecule-pd0325901). These checks establish neither historical-lot authentication nor clinical availability.

## Positive outcomes and their denominators

| Experiment | Result | Qualification |
|---|---|---|
| Population colonies, Figures 1A-B | Approximately 1,900 NANOG-positive colonies and 42% efficiency by day 6 with A2S; approximately 3.2% in FBS, most independent of continued induction by day 12 | Initially 5,000 cells were plated, but the actual day-0 adherent count was the efficiency denominator. Colonies merged after day 6. Bars show SD from two replicate samples. |
| Withdrawal resistance | A2S colonies were independent of continued induction after day 4, compared with day 9 in FBS | Assessed after four additional days without doxycycline and chemicals; stronger than marker expression under continuing OSKM induction. |
| Single-cell reprogramming, Figure 1C | Approximately 40% AP-positive wells under A2S, in two independent experiments | One tdTomato-positive cell sorted per well; denominator was wells with tdTomato-positive cells. **Withdrawal was day 11 and scoring day 15**, not day 6. Exact per-experiment integers were not recovered without the image. |
| Established lines, Results referring to S1A-B | Day-6-derived colonies could be passaged in FBS, were reported karyotypically normal, and produced teratomas containing all three germ layers | Genuine reported differentiation/cytogenetic positives. Supplement images and clone, metaphase and animal denominators were not accessible. |
| Component omissions, Figure 6A | SGC0946+2i yielded roughly half the A2S colony count; AA+2i and AA+SGC0946 gave 6.6% and 10.4% efficiency on day 6 | Dual-combination colonies remained NANOG-positive after withdrawal. SD from two replicate samples. Individual-component effects were described as lower, with data not shown. |

The paper's three factors are AA, SGC0946 and **2i as a two-compound mixture**. Thus A2S contains four compounds, while the “dual” combinations omit one factor. CHIR99021 and PD0325901 were not individually separated within 2i. Preserve combination-level benefits without assigning the full effect to each compound. No colony-comparison p-value or formal drug-interaction synergy metric was found in the read main text/legends.

## BUB1B and functional boundaries

`Bub1b` is a single-cell RNA expression example alongside Mcm6 and Ccnb1 in cell-cycle gene groups. The authors explicitly consider selection of cells already cycling as an alternative to transcriptional induction. A2S retained cycling transcript patterns and Ki67 positivity more often than FBS. No BUBR1 protein amount, stability, checkpoint function or causal Bub1b perturbation was measured; the functional siRNA experiments concern Ehf, Eif4a1 and Phlda2.

Component omissions and inferred expression networks support suppression of somatic programs by 2i and cooperation between AA and SGC0946 in pluripotency-network activation. Network edges are statistical associations; study-specific drug target engagement and selectivity were not established. Ink4a remained activated and p53 transcripts persisted, so increased cycling is not proof of complete senescence removal or tissue rejuvenation.

Normal endpoint karyotype is a positive, but its sampling depth remains unresolved. It does not establish faithful segregation at every division, repair of pre-existing aneuploidy, absence of sequence-level mutations or functional health of all daughter cells. This study did not establish germline contribution, chimera competence or tetraploid complementation by its derived lines; general claims in the introduction refer to other studies.

Teratoma formation is the deliberate pluripotency test, not evidence of reduced tumor risk. The inducing cassette includes c-Myc. Neither tumor-free therapeutic use nor malignant-cancer causation was demonstrated. There is no BUB1B-variant model, skeletal-muscle strength/regeneration endpoint or clinical efficacy estimate.

## Source limits and canonical compounds

Complete continuous primary text, STAR Methods, key resources and all seven main figure legends were read. **No main figure image or supplement was successfully accessed.** In particular, normal karyotype and teratoma findings remain Results reports about S1A-B, not independently inspected panels. There is no raw-data reanalysis or additional abstract-reading credit. Assessed facts and source hashes are in [the companion JSON](reprogramming-qualification.json).

SGC0946 must not become an alias of `new-epz5676`. PD0325901 belongs to `new-mirdametinib`. The other compound identities are catalog-qualified ascorbate 2-phosphate, SGC0946 and CHIR99021. Each may be tracked as a protocol component, but the same complete-cocktail efficacy must not be counted as four independent compound successes.
