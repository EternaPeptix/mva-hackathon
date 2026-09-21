# Public transfer result

Two worked sequence-to-mechanism-to-route examples: the existing BUB1B research hypothesis and an explicitly synthetic TRIP13 homozygous genotype from public literature. No second patient was analyzed, and no drug response was measured.

| Example | Input and computed consequence | Curated mechanism | Route result |
| --- | --- | --- | --- |
| BUB1B: existing hypothesis reconstructed from public reference sequence | Two heterozygous coding SNVs; trans hypothesis, phase unresolved; c.2210T>G: TTA→TGA, p.(Leu737Ter), 736 aa if translated; c.3006T>G: AAT→AAG, p.(Asn1002Lys), 1050 aa if translated | BUBR1 insufficiency is supported for selected MVA1 alleles. The actual N1002K abundance/function and phase remain unmeasured. | Niacin / NAD-SIRT2-BUBR1 stabilization: **conditional_research_route**; Arimoclomol / proteostasis restoration: **evidence_gap**; Everolimus / mTORC1 probe: **evidence_gap** |
| TRIP13: synthetic homozygous public genotype; no patient phenotype or reads | c.[1060C>T];[1060C>T] is stipulated for this example; c.1060C>T: CGA→TGA, p.(Arg354Ter), 353 aa if translated | Published Arg354Ter cells lacked detectable TRIP13 and had checkpoint/segregation defects rescued by WT TRIP13; mutant rescue failed. MAD2 remodeling supplies a distinct checkpoint mechanism. | Niacin / NAD-SIRT2-BUBR1 stabilization: **not_transferred**; Arimoclomol / proteostasis restoration: **not_transferred**; Everolimus / mTORC1 probe: **evidence_gap** |

## BUB1B: evidence gaps

NM_001211.6: 3153 nt CDS including termination; 1050 reference amino acids. RefSeq translation and complete Ensembl CDS agree.

**Niacin / NAD-SIRT2-BUBR1 stabilization — conditional_research_route.** Confirm trans and an abundance-limited, functional N1002K substrate; measure niacin exposure, NAD/SIRT2 engagement, BUBR1 mediation, segregation fidelity and productive muscle output. A full-length sequence is not a functional protein measurement.

- causal_protein_matches_stabilization_target: true — Conditional on the proposed pair, the affected protein is BUBR1; selected human MVA1 mutants show abundance defects (Suijkerbuijk 2010).
- full_length_canonical_product_possible: true — Computed canonical products; protein abundance and activity remain unmeasured.
- proximal_pathway_to_target_bridge: true — North 2014 supplies an NMN/SIRT2-to-BUBR1 bridge; Pirinen 2020 supplies niacin-to-human-muscle-NAD engagement. Joining these interventions is a hypothesis, not niacin/BUB1B rescue.

**Arimoclomol / proteostasis restoration — evidence_gap.** Demonstrate actual-mutant function at justified exposure before promoting proteostasis support beyond a reserve hypothesis.

- full_length_canonical_product_possible: true — Computed canonical products; protein abundance and activity remain unmeasured.
- actual_mutant_rescued_by_drug_at_relevant_exposure: none — No N1002K/arimoclomol rescue supplied. HSP90 dependence of selected BUBR1 mutants does not establish rescue by this medicine.

**Everolimus / mTORC1 probe — evidence_gap.** Qualify harmful excess and causal partial normalization first; then measure drug engagement, tissue output and renewal. Keep as a probe, not an automatically transferred rescue.

- harmful_excess_mtorc1_in_relevant_model: none — Current proposal has only qualitative carrier-mouse phosphoprotein bands and no confirming muscle RNA result; no affected/corrected human model establishes harmful excess.
- normalization_preserves_productive_tissue_output: none — No drug or genetic-normalization experiment has measured productive output in the proposed model.

Resolve trans; confirmed cis retires the proposed pair.

No qualified clinical drug lead is produced.

## TRIP13: evidence gaps

NM_004237.4: 1299 nt CDS including termination; 432 reference amino acids. RefSeq translation and complete Ensembl CDS agree.

**Niacin / NAD-SIRT2-BUBR1 stabilization — not_transferred.** Test TRIP13-specific correction and successful chromosome segregation as a positive control; a BUBR1 bypass would require direct evidence. Non-transfer is an evidence decision, not proof that a bypass is impossible.

- causal_protein_matches_stabilization_target: false — The sequence lesion and published genetic rescue concern TRIP13, not BUBR1. A shared checkpoint phenotype does not establish that raising BUBR1 substitutes for TRIP13.
- full_length_canonical_product_possible: false — Computed canonical products; protein abundance and activity remain unmeasured.
- proximal_pathway_to_target_bridge: none — The bounded source set supplies no NAD/SIRT2-to-TRIP13 rescue or validated BUBR1 bypass.

**Arimoclomol / proteostasis restoration — not_transferred.** Any truncation-rescue or readthrough mechanism needs its own transcript, protein-function and exposure evidence. It is not the full-length missense-stabilization route transferred from BUB1B.

- full_length_canonical_product_possible: false — Computed canonical products; protein abundance and activity remain unmeasured.
- actual_mutant_rescued_by_drug_at_relevant_exposure: none — Published WT rescue does not show pharmacological rescue of Arg354Ter; its expressed mutant failed checkpoint rescue.

**Everolimus / mTORC1 probe — evidence_gap.** Establish a relevant causal signaling abnormality and functional rescue before nominating this downstream drug. Shared MVA labeling is insufficient.

- harmful_excess_mtorc1_in_relevant_model: none — No TRIP13-genotype-specific harmful mTORC1 excess is established by the selected primary sources.
- normalization_preserves_productive_tissue_output: none — No relevant output-preserving rescue is supplied.

Homozygosity is stipulated for the synthetic example, not observed here.

No qualified clinical drug lead is produced.

## Evidence sources

- **yost2017:** [Yost et al., Biallelic TRIP13 mutations predispose to Wilms tumor and chromosome missegregation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5493194/) — Table 1, Figs. 1–3 and supplements: public Arg354Ter genotype, RNA/protein findings and genetic rescue; no repurposed-drug rescue.
- **kim2018:** [TRIP13 and APC15 drive mitotic exit by turnover of interphase- and unattached kinetochore-produced MCC](https://www.nature.com/articles/s41467-018-06774-1) — TRIP13 loss, MAD2 state, checkpoint activation and silencing; mechanistic context, not this synthetic example's experimental result.
- **clinvar431045:** [ClinVar VCV000431045.1](https://www.ncbi.nlm.nih.gov/clinvar/variation/431045/) — Pathogenic, one literature-only OMIM submission, no assertion criteria (zero stars); exact MANE variant. This is not expert-panel review or independent evidence from Yost.
- **suijkerbuijk2010:** [Suijkerbuijk et al., Molecular causes for BUBR1 dysfunction in the human cancer predisposition syndrome mosaic variegated aneuploidy](https://pubmed.ncbi.nlm.nih.gov/20516114/) — Selected human MVA1 mutant abundance defects and expression rescue; N1002K was not tested.
- **north2014:** [North et al., SIRT2 induces the checkpoint kinase BubR1 to increase lifespan](https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/) — NMN abundance experiments and separate genetic SIRT2 experiments; neither proves niacin/N1002K rescue.
- **pirinen2020:** [Pirinen et al., Niacin Cures Systemic NAD+ Deficiency and Improves Muscle Performance in Adult-Onset Mitochondrial Myopathy](https://pubmed.ncbi.nlm.nih.gov/32386566/) — Human muscle NAD engagement in another disease; strength also improved in treated healthy controls, so disease-specific functional rescue is unestablished.

Only sequence consequences and gate execution were recomputed. This is two worked examples, including one synthetic public genotype; it is not two-patient biological validation, a drug efficacy benchmark or an automated literature review.
