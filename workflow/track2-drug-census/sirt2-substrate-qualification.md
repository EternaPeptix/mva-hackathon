# NPD11033 and substrate-selective SIRT2 inhibition

**NPD11033 is a useful mechanistic probe with real pancreatic-cancer-cell growth/viability activity and a direct, separately measured reduction of BubR1 protein in HeLa cells.** Its inhibition of SIRT2 deacetylation while sparing the tested demyristoylation reaction shows why “SIRT2 inhibition” is not a single interchangeable intervention. It does not establish normal-cell safety, MVA rescue or a clinical effect of niacin.

## Source and exact reading depth

Kudo N, Ito A, Arata M, Nakata A, Yoshida M. 2018. *Identification of a novel small molecule that inhibits deacetylase but not defatty-acylase reaction catalysed by SIRT2.* [PMID 29685974](https://pubmed.ncbi.nlm.nih.gov/29685974/), [PMC5915714](https://pmc.ncbi.nlm.nih.gov/articles/PMC5915714/), [DOI 10.1098/rstb.2017.0070](https://doi.org/10.1098/rstb.2017.0070).

Complete main Introduction, Methods (a-h), Results (a-f), Discussion, Table 1 and all five main figure legends were read through the [publisher-provided full-text preview on ResearchGate](https://www.researchgate.net/publication/324705445_Identification_of_a_novel_small_molecule_that_inhibits_deacetylase_but_not_defatty-acylase_reaction_catalysed_by_SIRT2), which explicitly attributes the preview to the Royal Society. Overlapping primary PMC indexed sections corroborate the material findings. Three saved preview extracts cover the continuous main text through the end of the Discussion/data-access statement; the reference list was not used as a substitute for reading cited papers.

**No original figure pixels or complete PDF were successfully retrieved or visually inspected.** Figure-axis/lane values below are extracted text, not independently verified images. Normal PMC/publisher requests returned challenge HTML/403, Europe PMC XML failed and BioC returned no result; no challenge was bypassed. No separate supplement was identified in the accessible article or Europe PMC metadata (`hasSuppl=N`); no supplementary reading is claimed. PDB 5Y0Z/ligand 8K9 deposited metadata were checked, but crystal coordinates/raw diffraction data were not reanalysed.

## Identity, biochemical activity and selectivity

NPD11033 is a **synthetic cytisine-derived compound**, not cytisine itself. The study bought NPD11033 and inactive analogue RK-0310020 from Pharmeks through a local distributor. In RK-0310020, the bulky two 1,1-dimethylpropyl groups are replaced by methyl groups. Cytisine's pharmacology or clinical use cannot be transferred to NPD11033.

| Experiment | Finding | Controls and limits |
|---|---|---|
| SIRT2 deacetylation | IC50 **0.46 ± 0.056 micromolar**, using fluorescent acetyl-lysine peptide/electrophoretic mobility shift after a fluorogenic discovery screen | Figure 1: three independent assays; Methods derive IC50 means/SD from at least three independent curves. SIRT1/SIRT3 had only marginal inhibition; HDAC1/HDAC6 were not inhibited under tested conditions, with TSA as active control. This is a limited enzyme panel, not proteome-wide selectivity. |
| Inactive analogue | RK-0310020 lost deacetylase-inhibitory activity; extracted Figure 3 label gives IC50 >50 micromolar | Three independent assays. This matched-structure control supports mechanism but does not establish equal cellular uptake or rule out every off-target effect. |
| SIRT2 demyristoylation | NPD11033 did not inhibit the tested reaction even at **1 millimolar** | Fluorogenic myristoyl-lysine peptide assay; three independent assays. SirReal2 also failed to inhibit it, while suramin inhibited both deacetylation and demyristoylation. Do not generalize to every acyl substrate or cellular defatty-acylation pathway. |
| Structural binding | SIRT2 catalytic domain, residues 52-356, co-crystallized with NPD11033 at **2.0 Å** resolution; ligand induces a hydrophobic pocket and contacts His187 at a reported 2.8 Å hydrogen-bond distance | Purified protein 17 mg/mL with 1 mM ligand for crystallization. Deposited [PDB 5Y0Z](https://www.rcsb.org/structure/5Y0Z). This is direct purified-protein binding evidence, not cellular occupancy or a measured binding-residence time. |

The demyristoylation assay used recombinant SIRT2, 50 micromolar Ac-RHKK(Myr)-MCA, 0.1 mM NAD, 20 microlitres of 50 mM Tris-HCl pH 9.0, 4 mM MgCl2 and 0.2 mM DTT for 60 minutes at 37°C, followed by trypsin development and fluorescence. Exact deacetylase substrate/buffer details are largely delegated to earlier methods references rather than fully specified here; those prior papers were not newly appraised. Distinct substrates/assay conditions preclude treating the apparent potency ratio as a universal selectivity constant.

The authors propose alternative structural explanations for retained demyristoylation (ligand retention versus displacement by the acyl chain); neither is resolved. The structure suggests an effect on the cofactor-binding loop but does not itself establish a kinetic NAD-competitive mechanism.

## Cellular experiments: keep cell type, timing and co-treatment separate

| Arm | Exposure and comparator | Actual result and denominator |
|---|---|---|
| **PANC-1 eIF5A acetylation**, Figure 5A | NPD11033 concentration series; extracted lane labels 1, 3, 10, 20 micromolar. **1 micromolar TSA co-treatment for 8 h**. RK-0310020 20 micromolar and SirReal2 20 micromolar controls | Increased acetyl-eIF5A by immunoblot with NPD11033; inactive analogue did not reproduce it. Total eIF5A was measured. Authors explain that concurrent HDAC6 and SIRT2 inhibition is required for a marked signal. This is not evidence that NPD11033 alone produced the same acetylation change. Blot replication/densitometry denominator not stated. |
| **HeLa BubR1 protein**, Figure 5B | **24 h** NPD11033; extracted lane labels 3, 10, 30 micromolar. DMSO control label, RK-0310020 **30 micromolar**, and nicotinamide **20 millimolar** positive control | Dose-dependent BubR1 reduction by immunoblot, with alpha-tubulin loading control; inactive analogue did not reproduce reduction. This is direct protein-abundance evidence in HeLa, **not a PANC-1 measurement**. No quantitative densitometry/replicate n, half-life, acetyl-BubR1, ubiquitination, BUB1B mRNA or segregation assay reported. |
| **PANC-1 pharmacology**, Figure 5E | **72 h** NPD11033, SirReal2 or RK-0310020. Extracted concentration labels 0, 0.1, 0.3, 1, 3, 10, 30, 60 micromolar; no TSA co-treatment specified for this arm | ATPlite signal decreased; NPD11033 IC50 **7.5 ± 1.3 micromolar**. SirReal2 also reduced viability; inactive analogue did not reduce growth through 60 micromolar. Three independent assays, SD. Methods plate **5,000 cells/well**. |
| **PANC-1 SIRT2 knockdown**, Figures 5C-D | Two SIRT2 siRNAs (s22707/s22708) versus non-targeting AllStars control; RNAiMAX. Protein assessed day 5 and ATPlite day 7 | Reduced SIRT2 protein and viability; **four independent assays**, SD. Student's t-test. Methods label **p<0.01**, legend prints **p=0.01**; preserve this reporting discrepancy rather than invent an exact p-value. siRNA dose was not found. |

The article uses “in vivo inhibition” to describe **cultured-cell** experiments. It does not report drug-treated animals. Likewise, its claim of PANC-1 killing is based here on **ATP-based viable-cell signal**, which cannot alone separate cytostasis, altered per-cell ATP and cell death. Preserve the meaningful concentration-dependent anti-growth/viability result, but do not upgrade it to independently confirmed apoptosis, clonogenic eradication, tumour shrinkage or organism survival.

Exact panel doses are supported by extracted publisher figure text; original pixels remain uninspected. The main Results refer to Figure 5C when describing inactive-analogue pharmacology, while the legend assigns that dose-response experiment to **5E**. DMSO is labelled in the extracted figure, but final vehicle concentration and matching across all arms are not specified. No normal-cell comparator, exposure/PK, or therapeutic-index experiment was found. An inactive compound's lack of PANC-1 activity does not demonstrate NPD11033 safety.

Taken together, direct structural binding, enzyme assays, inactive analogue, active comparator and two-siRNA phenocopy support SIRT2-mediated activity. They do not close cellular target exclusivity: there is no resistant-target rescue, inhibitor-resistant mutant, cellular occupancy measurement or direct proof that BubR1 loss causes PANC-1 growth inhibition.

## Implication for the niacin/NAD/SIRT2/BUBR1 hypothesis

The HeLa result is consistent with the **direction** of a SIRT2-deacetylation/BUBR1-abundance relationship: pharmacological inhibition lowers BubR1. It strengthens the case for measuring BUBR1 protein when perturbing this pathway. It does **not** demonstrate the reverse chain from nicotinic acid to NAD, effective SIRT2 activity, usable variant BUBR1, faithful daughter cells and muscle function. No niacin treatment or BUB1B-variant rescue was performed. Here **“NA” means nicotinamide at 20 mM**, not nicotinic acid.

Cancer-cell anti-growth activity is a legitimate positive in its own model, not a reason to infer either benefit or harm in every tissue. Conversely, BubR1 reduction in HeLa does not make all SIRT2 inhibitors universally harmful, nor prove that SIRT2 activation is universally beneficial. Substrate, cell type, lesion, concentration, duration and intervention matter. The paper itself discusses both tumour-suppressive and tumour-promoting SIRT2 contexts from other studies; those citations were not newly qualified here. No numerical MVA efficacy or clinical safety estimate is defensible from this source.

## Current research supply and decision reconciliation

Current [Probechem PC-24697](https://www.probechem.com/products_NPD11033.aspx) lists NPD11033 for laboratory use with quotation-only packs. [MedKoo 564325](https://www.medkoo.com/products/30255) explicitly states currently out of stock with possible custom synthesis. These are **research supply listings**, not a verified immediately available batch, clinical product or study-lot authentication. No supplier was contacted and no purchase was made.

The [deposited ligand 8K9](https://www.rcsb.org/ligand/8K9) has formula C30H44N2O3, molecular weight 480.682, and fully specified stereochemical InChIKey `ITCUJNFERKWPKY-WMTXJRDZSA-N`. MedKoo's listing uses `ITCUJNFERKWPKY-JCVDRHSJSA-N`, with the hydroxyl-bearing centre unspecified in its InChI. Same name/connectivity does not authenticate stereochemical composition or biological equivalence of a current batch; this remains a COA/identity question, not proof that a vendor material is wrong. No verified CAS identifier was found in the checked primary/depositor/supplier records.

The ledger retains **`new-npd11033`** and a separate **`new-sirreal2`** entry. Suramin is a biochemical comparator; RK-0310020 is the inactive control. TSA remains the existing **`new-trichostatin-a`** entry. These control roles are preserved without treating every assay ingredient as a therapeutic nomination. Niacin remains provisional.

## SirReal2: a distinct active comparator

| Primary experiment | Qualified result and boundary |
|---|---|
| Figure 5E, PANC-1, 72 hours | SirReal2 reduced ATP-based viability in the concentration-response experiment. Three independent assays, SD; Methods plate 5,000 cells/well. Extracted axis values are 0, 0.1, 0.3, 1, 3, 10, 30 and 60 micromolar; original pixels uninspected. No numeric **SirReal2** cellular IC50 is reported. NPD11033 and inactive RK-0310020 are comparators. No TSA is specified in this arm. Final vehicle percentage is unstated. |
| Figure 5A, PANC-1, 8 hours | SirReal2 **20 micromolar plus TSA 1 micromolar** is the active comparator for acetyl-eIF5A. Total eIF5A is measured; blot replication/densitometry n is unstated. This supports the combined exposure, not an equal SirReal2-alone response. |
| Figure 4B-C, purified enzyme | SirReal2 inhibits SIRT2 deacetylation but not the tested demyristoylation reaction. This is substrate-qualified biochemical evidence, not every defatty-acylation substrate. |

The ATP result is a real cell-level positive but does not distinguish cytostasis, death and altered per-cell ATP. No normal-cell safety, drug-treated animal benefit, BUB1B-variant rescue, faithful chromosome segregation or muscle outcome is demonstrated here. The NPD11033 7.5 ± 1.3 micromolar cellular IC50 and NPD11033 co-crystal are not SirReal2 results.

**The HeLa BubR1-lowering arm is NPD11033, not SirReal2.** Figure 5B does not test SirReal2; do not transfer that observation to this decision. Supplier biological descriptions are not independently qualified primary experiments in this addendum.

The study names Sigma as the SirReal2 supplier without a catalog number or lot. The current [Sigma-Aldrich SML1514 listing](https://www.sigmaaldrich.com/US/en/product/sigma/sml1514) identifies **SirReal2 / SirReal 2**, also described as Sirtuin-rearranging ligand 2; CAS **709002-46-0**, formula **C22H20N4OS2**, molecular weight **420.55**, InChIKey **MENNDDDTIIZDDN-UHFFFAOYSA-N**. SML1514 is a catalog alias, not a second compound. The current listing is research-reagent evidence, not clinical availability or historical-lot authentication. Supplier DMSO solubility is not the study's final vehicle concentration.


A second reader checked Figure 5's extracted primary legend and the ATPlite method: HeLa and PANC-1 arms, timings, co-treatment and independent-assay counts remain distinct. Original image pixels were not available to either reader.
