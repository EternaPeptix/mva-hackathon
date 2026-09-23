# EGFR inhibitors and naproxen: bladder-tumor prevention, not BUBR1 restoration

**23 September 2026.** [Lubet and colleagues](https://doi.org/10.3233/BLC-200423)
report actual drug benefits in a chemically induced rat bladder-cancer model.
The paper was published in **2021**; its later PubMed indexing is not a new
experiment. These results support a tumor-prevention research comparison, but
do not displace niacin or apitegromab for the proposed MVA muscle question.

## What was measured

Female F344 rats received OH-BBN and then daily or weekly drug treatment. Drugs
started after carcinogen exposure; follow-up ended 28 weeks after the last
carcinogen treatment. Endpoints were bladder-plus-lesion mass and development
of large, palpable tumors, operationalized with a bladder-weight threshold
above 200 mg. These are not survival, checkpoint fidelity, normal-muscle output
or regression of confirmed baseline tumors. The paper's inference that microscopic
carcinomas already existed at later treatment starts comes from earlier work.

| Intervention | Actual result retained | Decision for the MVA comparison |
|---|---|---|
| Lapatinib | In the late-start daily arm, mean bladder-plus-lesion mass was 197 versus 356 mg and large-tumor incidence 26% versus 69% in controls. Weekly dosing also reduced mass; not every incidence comparison was significant. | Tumor-prevention comparator. No BUBR1 restoration or constitutional-cell benefit demonstrated. |
| Naproxen | Daily monotherapy: mean mass 185 mg and large-tumor incidence 23%, versus 356 mg and 69% in controls. | Retain this common medicine's actual prevention result. Its availability does not establish an MVA mechanism or suitability for the child. |
| Lapatinib plus naproxen | Table 2 reports 0% large tumors among 29 rats and mean mass 137 mg. Figure 3B remains at 100% free of the large-tumor endpoint. | A real combination result, but not proof that every histological lesion disappeared, that the combination is synergistic, or that benefit persists after treatment. |
| Erlotinib | Daily late-start treatment: mean mass 243 versus 394 mg; large-tumor incidence 27% versus 74%. | Tumor-prevention comparator with a tumor Bub1b RNA decrease, not an abundance-restoration route. |
| Gefitinib | This paper reuses earlier prevention data and includes short-treatment expression data. | Keep the measured RNA finding separate from newly performed prevention experiments; earlier efficacy data are not an independent replication here. |
| NO-naproxen | Separately tested for short-treatment RNA changes; earlier prevention results are summarized in Supplementary Table 1. | Distinct intervention from naproxen. Neither the naproxen combination outcome nor present-study efficacy is assigned to it. |
| Zileuton | Short-treatment RNA changes are reported. Lack of prevention is described as unpublished data, without a new efficacy table here. | Mechanistic comparator only. A pathway annotation does not establish functional benefit. |

The rat schedules are source descriptors, not proposed human regimens. The
combination used the same daily lapatinib and naproxen regimens as the respective
monotherapies. No formal interaction test establishes synergy. The control and
lapatinib groups repeated across Tables 1 and 2 are not independent replications.

## The gene-specific result is in the downward direction

The RNA experiment was a separate five-day treatment experiment in animals with
palpable tumors. Methods state five arrays per group; discussion says five to
six. Gefitinib data came from earlier work. There was **no lapatinib–naproxen
combination RNA arm**, and the NSAID RNA arm used **NO-naproxen**, not naproxen.

Targeted inspection of the original supplementary workbooks gives:

| Treatment | Native workbook and cells, Sheet1 | Control → treated expression means | Listed P value | Listed signed fold change |
|---|---|---:|---:|---:|
| Erlotinib | Supplementary Table 3, `A20:F20` | 6.3032 → 5.5789 | 0.0053614 | −1.6521 |
| Gefitinib | Supplementary Table 4, `A19:F19` | 6.3031915 → 5.407627 | 0.0003605458 | −1.8603377 |

These are **Bub1b RNA decreases in tumors**. The reported signed fold-change
convention must not be read as negative expression. The stored means, P values
and fold changes are retained without relabeling the P values as adjusted.
Methods, results and table footnotes describe differing significance/filtering
rules; the exact correction applied to these per-gene spreadsheet P values is
not resolved. No BUBR1 protein, kinase/checkpoint function, direct target
attribution or MVA genotype was tested.

Bub1b is absent from the selected altered-gene lists for lapatinib, NO-naproxen
and zileuton. **Absence from a filtered list is not evidence of no effect.**
Supplementary Table 6 is a union of genes across three separate EGFR-inhibitor
treatments, not a combination-drug arm. Its Bub1b entry cannot be assigned to
all three drugs. Sarcomere/myofibril pathway labels in tumor RNA do not measure
muscle performance.

## Source and transfer limits

- Main Table 2 gives control/lapatinib group sizes of 34/29, while Supplementary
  Table 1 prints 24/35 and 8/30 as the corresponding large-tumor fractions.
  These denominators remain unresolved; no integer counts are reconstructed from
  the rounded percentages and the two tables are not silently reconciled.
- Methods and Supplementary Table 1 describe 150 mg OH-BBN per gavage twice
  weekly for eight weeks. Figure 3's caption instead prints 150 mg/kg weekly
  and inconsistent age/timing wording. This prevents presenting a silently
  repaired dosing protocol.
- Figure 3 plots freedom from the large-tumor endpoint, despite its printed
  y-axis wording suggesting the opposite direction. The caption and results
  support the tumor-free interpretation; this label issue is retained.
- Randomization and blinding are not specified in the inspected report.
  Unshown body-weight observations do not establish comprehensive safety.
  Weekly-versus-daily comparisons are not a demonstrated clinical safety or
  noninferiority result.
- Basal bladder tumors induced by a carcinogen are not constitutional MVA,
  rhabdomyosarcoma, or nonmalignant patient muscle. Tumor suppression and lower
  cell-cycle RNA can be desirable in this model without repairing BUBR1
  insufficiency. This is a reason to keep the indication and endpoint explicit,
  not evidence that these drugs necessarily harm MVA cells.

## Reading scope

Root read the complete scientific main XML, original supplementary Tables 1–2,
native Figure 3, and the headers and Bub1/Bub1b-matched cells in all six original
supplementary workbooks. Workbook search is not a claim to have scientifically
reviewed every unrelated gene row. Raw files and extraction outputs remain
outside the repository; compact source hashes and exact cell references are
recorded in the companion JSON. Linked earlier papers are not credited as read.
This qualification adds no abstract-reading or independent-experiment count.
