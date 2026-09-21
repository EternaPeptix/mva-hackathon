# jeweled — Track 2: testing two ways to preserve muscle

**Research submission draft, revised 21 September 2026; not submitted.**

We nominate **everolimus as the narrow research lead and niacin (nicotinic acid) as a serious comparator**, tested separately. The aim is to preserve useful muscle and its ability to renew in a child with recorded muscle atrophy and failure to thrive. Everolimus would reduce an excessive growth signal if that signal contributes to the disease. Niacin would try to increase the supply of usable BUBR1, a protein that helps cells distribute chromosomes correctly when they divide.

Neither drug has demonstrated benefit for the nominated genotype. Everolimus leads on established drug-to-target pharmacology and pediatric treatment experience; niacin has the stronger argument for addressing a possible shortage of functional BUBR1. This is a choice of hypotheses to investigate, not demonstrated superiority. The [route decisions](appendices/track2-route-decisions.md) explain the evidence, counterevidence and reconsideration conditions for every alternative reviewed.

## What the genetic result means

**BUB1B is the gene; BUBR1 is its protein product.** Our [Track 1 hypothesis](jeweled_track1_report_submission.md) involves two changes: Leu737Ter, which predicts an early stop and a 736-amino-acid protein, and Asn1002Lys, which changes one protein building block. The second change remains a variant of uncertain significance: its functional effect is unmeasured.

We also do not know whether the changes affect opposite gene copies or the same copy. In the latter arrangement, ordinary translation stops before reaching Asn1002Lys, while the other copy could produce normal protein. That does not establish a sufficient amount of working BUBR1. RNA decay, protein abundance and actual function require measurement. Results for other BUB1B mutations cannot be assigned to these two variants.

## Why everolimus leads narrowly

**mTORC1 is a cellular growth-and-nutrient signaling system.** Too much activity could impair muscle maintenance, but too little could also compromise growth and repair. In a [different BUB1B mouse model](https://www.jci.org/articles/view/126863/version/2), muscle showed increased chemical activation of two mTORC1 outputs. The key signaling measurement used only two mice per genotype. It establishes a small disease-relevant association, not that suppressing the pathway improves muscle.

[Everolimus](https://dailymed.nlm.nih.gov/dailymed/fda/fdaDrugXsl.cfm?setid=2150f73a-179b-4afc-b8ce-67c85cc72f04) inhibits this system and has pediatric formulations and monitoring experience in other diseases. The proposed test is partial normalization toward the corrected-cell range. A genetic reduction of pathway activity would first help distinguish a harmful signal from an accompanying marker.

Counterevidence matters. A [preliminary rapamycin observation](https://pmc.ncbi.nlm.nih.gov/articles/PMC3785294/) reported no benefit for age-related features in another BubR1 mouse model, without enough treatment detail for a decisive comparison. [Everolimus also lowered BUBR1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10149375/) in late-passage, nonmutant umbilical-cord cells. These observations make preserved muscle growth, renewal and BUBR1 function essential checks. An intermittent schedule is a testable option, not an established solution.

## Why niacin deserves the direct comparison

**NAD is a molecule cells use for energy chemistry and several regulatory enzymes.** Niacin can supply material for making it. In [North and colleagues' experiments](https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/), the NAD-dependent enzyme SIRT2 stabilized BUBR1. A different NAD precursor, NMN, increased BUBR1 abundance; genetically increasing SIRT2 improved selected cardiac outcomes and male survival in BubR1-deficient mice. Those results were not niacin treatment or demonstrated skeletal-muscle rescue.

Niacin nevertheless has a human muscle connection. [Pirinen's small open-label mitochondrial-myopathy study](https://pubmed.ncbi.nlm.nih.gov/32386566/) corrected patient muscle NAD, but strength also improved in niacin-treated healthy controls. Patient knee-extension strength rose about **1.1-fold after ten months**; other muscle groups improved more. With no placebo group, this does not establish a disease-specific strength benefit. We use it as evidence of human muscle NAD engagement, not demonstrated disease rescue. A [pediatric mitochondrial-myopathy case](https://pubmed.ncbi.nlm.nih.gov/39173541/) had biochemical correction with little clinical improvement and continued progression. Correcting NAD therefore cannot serve as the success endpoint.

The missing chain is specific: can niacin reach and change the relevant muscle chemistry, increase functional BUBR1, and preserve productive muscle? More defective protein might be unhelpful. A benefit through another mechanism could still be valuable, but would not prove BUBR1 restoration. Niacin, nicotinamide, NMN and NR are different substances; their results are not interchangeable. Prescription-niacin pediatric safety and effectiveness are unestablished.

## What the completed computation established

The [sequence and population analysis](../workflow/track2-model/README.md) reconstructed possible protein products and explored declared assumptions. In one synthetic selective-killing example, the fraction of cells with the expected chromosome count rose **16.8 percentage points while modeled differentiated-cell time fell 11.0%**. This is a useful warning about measurement: a healthier-looking percentage can coexist with less tissue. It is not a drug effect or a patient forecast.

We also [reanalyzed public muscle RNA counts](../workflow/track2-rna/README.md): nine mice, 12,989 retained genes and 150 pathway comparisons. None passed the global statistical threshold. The mTORC1 comparison had P = 0.219 and adjusted false-discovery rate = 0.983. RNA measurements did not corroborate pathway activation; they do not disprove the published protein-activation measurement. This reused study data and supplies no independent biological replication. Neither analysis ranks drug efficacy.

## A staged experiment a partner team could undertake

No laboratory partnership or patient-derived model access is established. The proposed partners are a genetics center, a laboratory that measures muscle formation and contraction, a clinical pharmacologist and the treating oncology team. The [detailed experimental plan](appendices/track2-evidence-pipeline.md) defines the proposed comparisons and decisions. The first partner milestone is a reproducible muscle defect in an interpretable affected-versus-corrected model; drug confirmation follows successful model and exposure qualification.

1. **Qualify the model.** Resolve phase where possible; measure each allele's RNA, protein and function. Compare affected cells with genetically corrected and procedural controls. Engineered models must retain the unresolved alternatives. Initial fibroblast work can screen allele mechanisms economically, but muscle benefit needs a muscle model.
2. **Test mechanisms and exposures.** Establish whether excessive mTORC1 activity contributes to the defect, and whether niacin improves working BUBR1 through NAD. Test the drugs separately at measured, pharmacologically justified concentrations. In parallel, qualify whether senescent cells—cells persistently withdrawn from division—harm muscle renewal; same-gene genetic clearance evidence warrants this reserve investigation.
3. **Measure useful muscle, then confirm.** Follow both an initial population of renewable muscle progenitors and an already formed muscle population. Measure absolute contractile output and repeated-stimulation performance, alongside productive cell replacement and chromosome accuracy. Freeze one regimen per drug after a pilot; confirm with independent preparations, masked measurements and the specified correction for multiple comparisons.

Treatment must improve affected-cell output against matched vehicle controls. Responses in genetically corrected cells help distinguish a general functional benefit from a genotype-selective effect; improvement in both can still be useful, but does not by itself prove BUBR1 restoration.

**Improvement in either muscle population can qualify, provided the other population, renewal capacity and chromosome accuracy are preserved.** Both need not improve. Fewer division errors because cells stopped dividing or died do not count as repair. Reduced cellular-aging markers alone do not count as muscle benefit. A precise negative applies to the tested genotype, exposure and assay; inadequate exposure is an uninterpretable nonresponse. Niacin could overtake everolimus if it meets the same functional standard more convincingly.

## Scope and contribution

This is a research nomination and future experimental plan, with no patient drug-response data or dosing recommendation. Recorded rhabdomyosarcoma makes [tumor and treatment compatibility](appendices/track2-evidence-pipeline.md#rhabdomyosarcoma-evidence-and-treatment-compatibility) consequential: current tumor status and regimen are unknown, and neither candidate is assumed to protect against cancer. Renal function, interactions, growth, immune and metabolic effects require individual assessment before translation.

If independently validated, this comparison could distinguish a rescuable BUBR1 deficiency from harmful downstream growth signaling and identify which muscle-cell state benefits.

The innovation is the comparison: connect the actual allele products to two distinct intervention routes, then require greater muscle output without sacrificing replacement capacity or chromosome accuracy. The scalable component is a reusable sequence-to-mechanism-to-function decision process, supported by reproducible public-data analyses. Drug choices and exposures must be requalified for each genotype and tissue. The code and decision rules are reusable; new tissue models remain the principal experimental bottleneck. This addresses the challenge's [rigor, impact, innovation and scalability criteria](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/tabs/about.py) without claiming a new drug discovery.

Claude Code and Codex assisted literature review, analysis and writing on personal subscriptions, with training disabled according to the user; exact historical model/plan names are unverified. [Methods](jeweled_methods_form_answers.md#track-2-drug-repositioning) and [data provenance](PROVENANCE.md#dataset-provenance) distinguish public resources from controlled challenge data. Raw case reads and genome-scale payloads are excluded from this private repository. Team-authored reports carry CC BY 4.0; code also has an MIT option, with third-party and controlled-data exceptions in the [license](../LICENSE).

The [new pitch](jeweled_track2_pitch_reconciled.md) uses stock Skylar narration from Cartesia Sonic 3.6 on its Free plan, verified on 21 September 2026. Cartesia training/retention settings remain unverified. OpenAI image generation through Codex supplied a conceptual cover illustration, not biological evidence. Human review and public video upload are pending.
