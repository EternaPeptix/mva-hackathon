# Completing the outer abstracts: what exhaustion returns

**22 September 2026.** Fifty-seven disjoint assignments of up to 64 records each
covered every remaining saved outer abstract, 3,598 complete saved abstract
strings at unfiltered indices 1729–5326. Outer-corpus coverage of saved
abstracts is now **complete: 5,759 of 5,759 records that carry a saved abstract
string have been read across eight selections**, and 0 remain unread. That
sentence is narrower than it sounds and four limits belong next to it. **143**
remaining outer source/ID records carry no saved abstract text at all; nothing
was read for them and no annotation exists. Source/ID records include reviews,
preprints and publication versions, so 5,759 is not 5,759 independent
experiments. Abstract-level appraisal is not full-text review, and one record in
this tranche was taken to full text. And the count is the outer corpus only: the
1,460-record gene-query set is excluded from it and was read separately. This
pass flags **323 source questions**; a flag is not a qualified medicine.

## What an exhaustive search returned

3,598 abstracts produced **62 leads that survived source verification**, of which
**19 were confirmed as filed, 24 were overclaims, 17 were underclaims and 2 were
misattributions. Both misattributions are named: the
[trichostatin A record, PMID18326024](https://pubmed.ncbi.nlm.nih.gov/18326024/),
where the spindle collapse belongs to HDAC3 and N-CoR knockdown and only a
similar-defect statement attaches to the drug (house rule 2); and
[SP600125, PMID22438244](https://pubmed.ncbi.nlm.nih.gov/22438244/), where the
abstract identifies MPS1 and not JNK as the mediating target, so the record
corroborates no JNK hazard and its cross-reference to a second record is
withdrawn**; 44 carried at least one house-rule flag. By route class they
are 29 constraints, 18 fidelity records, 6 restoration, 5 compensation and 4
clearance. Five were filed tier A. Verification against the saved strings
downgraded four of them:
[C-type natriuretic peptide, PMID38136917](https://pubmed.ncbi.nlm.nih.gov/38136917/),
where every fidelity-relevant readout in goat cumulus-oocyte complexes is a
stained-morphology proportion or a transcript and the one functional endpoint is
oocyte maturation rate (house rule 1);
[Potentilla fulgens extract, PMID31740710](https://pubmed.ncbi.nlm.nih.gov/31740710/),
chemoprevention of areca-nut-induced stomach cancer in wild-type mice, reported
as a single 20% figure shared between two preparations, with p53 and Securin
measured as gene expression (house rule 6);
[metformin, PMID30022161](https://pubmed.ncbi.nlm.nih.gov/30022161/), where AMPK
phosphorylation stabilises TET2 but the abstract attaches no species, dose or
route to the metformin arm and it is the hyperglycaemia arm, not the drug arm,
that is qualified as shown in vivo (house rule 8); and
[resveratrol, PMID24194906](https://pubmed.ncbi.nlm.nih.gov/24194906/), where
protection of mouse oocytes from methylglyoxal rests on stains and fluorescent
probes and the abstract makes no sirtuin or NAD claim at all (house rule 1). One
tier-A record survived verification: UMK57.

Three triage lenses — beat niacin on protein supply, beat apitegromab on muscle
output, and disqualify for this patient in particular — were run over the same
62 leads. All three returned `displaces_current_nomination = false`. The
restoration lens records that **not one of the 62 leads measures BUBR1 abundance
at all**. The muscle lens records that a direct scan of the raw records for all
62 returns **zero occurrences** of muscle, myotube, myoblast, myofibre, atrophy,
body weight, lean mass, motor function, grip strength, contractility or
myostatin in any title or abstract. The safety lens records that every one of
the 62 is a constraint or an exclusion, and that the single tier-A record is
itself constrained rather than disqualifying.

That is the finding, and it should not be dressed as anything else. The value of
a negative changes with the size of the search: the nomination previously
survived a partial read of the outer corpus and now survives an exhaustive one.
What that licenses is a bounded claim — no administered medicine in the saved
outer abstracts restored faithful division in a BUBR1-limited cell, and none
supplied an absolute force gain in a growing animal. What it does not license is
a claim about the literature. Abstract depth can hide a full-text positive,
which is what the 323 flags exist to record, and the 143 records without saved
abstract text were never searched at all.

## What this tranche does to the restore/compensate distinction

On the restoration side the census closes with niacin's mechanism neither
corroborated nor contradicted by any administered agent. A regex count over the
3,598 saved strings finds BubR1 or BUB1B named in **two** records, neither a
drug-abundance experiment, and niacin or nicotinic acid named in **one**, a
colon-carcinogenesis review listing dietary compounds. The closest structural
match to niacin's own node is
[CX-4945, PMID34050247](https://pubmed.ncbi.nlm.nih.gov/34050247/): CK2 activity
promotes the interaction between CDCA3 and Cdh1, and pharmacological CK2
inhibition disrupts that degradation and elevates CDCA3. It is the only lead in
which a named administered small molecule raises the steady-state level of an
APC/C-Cdh1 substrate at the substrate-recognition step, and no mitotic block is
invoked anywhere, so it escapes house rule 4 in a way the blockade family does
not. It is still not a candidate: the protein raised is CDCA3, the CK2 effect is
described as substrate-specific, the abstract carries no quantitation of any
kind, and the stated purpose is to raise a genome-instability-correlated protein
to make tumour cells more killable by platinum, which is an adverse direction in
this child rather than merely an unhelpful framing.

The other end of the same axis is now closed rather than open.
[MLN4924, PMID31111756](https://pubmed.ncbi.nlm.nih.gov/31111756/) disables the
cullin-RING ligases, duly stabilises Emi1, and produces metaphase arrest,
impaired spindles and spindle-checkpoint over-activation in mouse oocytes: the
arrest signature in its purest form, with no functional benefit anywhere. With
it, the machinery-level strategy of blocking degradation to raise a mitotic
protein stands at five probes and no functional benefit — `prior-13`,
`new-mg132`, `new-protame`, `new-apcin` and this record. The inference is not
that each was a bad probe but that the level is wrong: a viable supply lever has
to be substrate-selective, which is niacin's shape and CX-4945's shape, not the
ligase's or the proteasome's. Support for that premise arrives by paralogy from
[BAY-320 and BAY-524, PMID26885717](https://pubmed.ncbi.nlm.nih.gov/26885717/),
which compare Bub1 catalytic inhibition against Bub1 protein depletion and find
that inhibition impairs chromosome arm resolution while exerting only minor
effects on mitotic progression or SAC function. In that kinase family checkpoint
competence tracks protein presence rather than catalytic output, which is the
axis the restoration route assumes. It is Bub1 and not BUB1B, in HeLa and RPE1
with normal BUBR1, so it is an argument about the premise and must not be banked
as a BUBR1 result; and it carries its own administered harm, sensitisation to
low-dose paclitaxel with impaired segregation and proliferation. A second, independent census
record on the same family runs against that comfort.
[BAY1816032, PMID38477771](https://pubmed.ncbi.nlm.nih.gov/38477771/) applied
the Bub1 inhibitor to primary bone-marrow-derived macrophages from wild-type
mice and increased osteoclastogenesis, reproducing pharmacologically, in
non-transformed primary cells, the direction of the paper’s myeloid Bub1
conditional-knockout arm, where femoral bone mineral density fell through
increased osteoclastic resorption. No chromosome-segregation or aneuploidy
readout appears and no animal was dosed, so this is a non-mitotic tissue harm
rather than a fidelity result — but bone is a live concern in a growing child,
and it qualifies the catalysis-versus-scaffolding comfort drawn from BAY-320.

The compensation side is not constrained by this tranche so much as untouched by
it. Across the 62 leads nothing measures muscle output; across all 3,598 saved
strings, 55 mention muscle at all and the word myostatin appears in none of
them. The single myostatin-pathway mention is GDF8, named alongside GDF11 and
GDF-15 in a proteomics review of aging skeletal-muscle satellite-cell
biomarkers, with no myostatin-pathway agent administered anywhere. The single record in which an
administered agent raised a tissue's regenerative output is
[Wy-14643, PMID42019507](https://pubmed.ncbi.nlm.nih.gov/42019507/), which
accelerated hepatocyte proliferation in mouse liver 48 h after 70% partial
hepatectomy; it is a rodent tool peroxisome proliferator, the tissue is adult
liver recovering from resection, the benefit is proliferation, and the stated mechanism is HIF1α-FOXM1 with post-hepatectomy recovery as
the improved outcome, neither of which reaches muscle, so what survives is a question about clinically usable PPARα agonists and a post-mitotic
force endpoint, not a lead. Two records land on standing compensation entries
and both run against them.
[Ursolic acid, PMID34858165](https://pubmed.ncbi.nlm.nih.gov/34858165/)
down-regulated PLK1 and CCNB1 protein by mediating p53 signalling in vitro and
in an in vivo arm with no stated species, dose, route or duration; the targets
were nominated computationally (house rule 3) and the readout is expression
(house rule 1), so this does not demote `new-ursolic-acid`, but it does pose the
dose-overlap question — whether mitotic-protein suppression occurs at or below
the myotrophic exposure — to a muscle-compensation candidate in a child whose
lesion is mitotic.
[Vismodegib, PMID35723318](https://pubmed.ncbi.nlm.nih.gov/35723318/) reduced
the liver regeneration rate after partial hepatectomy, decreasing Cyclin D and
accumulating Cyclin B: pharmacological Hedgehog blockade suppresses an
organ-renewal endpoint, while the implied compensation direction is Hedgehog
agonism, which is contraindicated in a child with recorded rhabdomyosarcoma.
Both directions are closed, which is a result rather than a gap.

The everolimus probe is unchanged for the plainest possible reason: no
mTOR-family agent appears among the 62 verified leads at all. That is a coverage
fact and not evidence in either direction.

## Source questions from this tranche

| Source question | Evidence credited at abstract depth | What remains to establish |
|---|---|---|
| [UMK57, PMID27829147](https://pubmed.ncbi.nlm.nih.gov/27829147/) | A named small molecule potentiating MCAK suppresses chromosome mis-segregation in CIN cancer cells at sub-lethal doses; mis-segregation rebounds within a few days through Aurora B changes that hyper-stabilise k-MT attachments, reversibly. The complete saved string, 1,121 characters, states the agent, the target, the mechanism, the rebound and the reversibility. | Taken to full text below. At abstract depth the open questions are direction and amount in a BUBR1-limited non-transformed cell, durability past 72 hr, and clinical eligibility. |
| [CX-4945, PMID34050247](https://pubmed.ncbi.nlm.nih.gov/34050247/) | CK2 activity promotes the CDCA3–Cdh1 interaction; CK2 inhibition disrupts CDCA3 degradation and elevates its level, with no mitotic block invoked. | Whether CK2 licenses Cdh1 recognition beyond CDCA3, and whether CX-4945 raises BUBR1 protein in a BUB1B-hypomorphic cell without raising mis-segregation. A broad Cdh1 lever would stabilise cyclin B and securin, which is the rule-4 signature. No fold-change, n, cell line or in vivo arm appears. |
| [MLN4924, PMID31111756](https://pubmed.ncbi.nlm.nih.gov/31111756/) | Neddylation inhibition stabilises Emi1 and arrests mouse oocytes at metaphase with impaired spindles and SAC over-activation. | Whether any sub-arresting dose raises an APC/C substrate without impairing the spindle. If not, machinery-level blockade is closed rather than open. The saved string drops the qualifier after "during", so metaphase I versus II is not determinable. |
| [BAY-320 / BAY-524, PMID26885717](https://pubmed.ncbi.nlm.nih.gov/26885717/) | Bub1 catalytic inhibition, compared directly against Bub1 depletion, impairs chromosome arm resolution with only minor effects on mitotic progression or SAC function, and sensitises cells to low-dose paclitaxel. | Whether the catalysis-versus-scaffolding dissociation holds for BUB1B: is SAC function preserved in a cell carrying a kinase-domain-mutant BUBR1 at normal abundance? This is paralogy, not a BUBR1 result. |
| [AICAR, PMID21315436](https://pubmed.ncbi.nlm.nih.gov/21315436/) | AICAR induces p53-mediated apoptosis in primary mouse embryonic fibroblasts trisomic for individual chromosomes; AICAR and an unnamed protein-folding inhibitor, especially combined, are efficacious against aneuploid human cancer lines. | The AMPK-activating exposure at which trisomic-cell killing begins, stated against metformin's paediatric plasma range, and whether the same apoptosis occurs in constitutionally aneuploid patient cells. Chloroquine is named only as a member of the identified set; the third compound stays unnamed under house rule 7. |
| [Apcin, PMID32152539](https://pubmed.ncbi.nlm.nih.gov/32152539/); [proTAME, PMID23145059](https://pubmed.ncbi.nlm.nih.gov/23145059/) | Apcin binds one Cdc20 site required both for substrate ubiquitination and for MCC-dependent APC/C inhibition, so it causes net inhibition and prolongs mitosis when SAC activity is low, and net activation and shortens mitosis when SAC activity is high, cooperating with p31comet. proTAME prolongs metaphase, and the prolonged metaphase causes loss of sister chromatid cohesion, which re-activates the SAC. | Which side of the inversion a BUBR1-hypomorphic cell with residual MCC occupies, measured as segregation error in daughter cells rather than as mitotic duration; and whether proTAME's cohesion loss requires an arrest this patient's cells cannot mount. High and low SAC activity were experimenter-set states, not a genetic MCC deficiency. |
| [Oxindole-1 with ZM447439 and ABT-263, PMID26188358](https://pubmed.ncbi.nlm.nih.gov/26188358/) | Two polyploidy-inducers; Ox-1 lowers CDK1/cyclin B1 activity, inhibiting Bcl-xL phosphorylation and conferring resistance to apoptosis. Adding ABT-263 to either gives synergistic loss of viability with greater sustained tumour-growth inhibition in AML lines and primary blasts, and knockdown of Bcl-xL but not Bcl-2 reproduces it. | Whether constitutionally aneuploid, non-transformed cells are Bcl-xL-dependent for survival at navitoclax, DT2216 and PZ15227 exposures. Drug-induced polyploidy in transformed leukaemia is not constitutional mosaic aneuploidy. |
| [Zoledronic acid, PMID20015195](https://pubmed.ncbi.nlm.nih.gov/20015195/) | Loss of prenylated Cenp-F from kinetochores and a significant increase in lagging chromosomes in three breast cancer lines, reversed by a mevalonate-pathway intermediate; no effect on Cenp-F in MCF-10A. | Whether paediatric dosing exposures raise mis-segregation in non-malignant and BUBR1-low cells. The same question is owed to `new-lonafarnib` and `new-FTI-277`. This is a realistic paediatric co-medication, not a nomination. |
| [Nutlin-3a, PMID18922897](https://pubmed.ncbi.nlm.nih.gov/18922897/) | A transient 24 h exposure followed by withdrawal produced endoreduplication and stable tetraploid clones more resistant to ionising radiation and cisplatin than their diploid parents, dependent on p53–p21. | Frequency in non-transformed BUBR1-low cells. Any p53-raising route needs a ploidy readout after drug withdrawal, not only a proliferation or apoptosis readout during exposure. Tetraploid clone frequency is never quantified. |
| [Genotoxic panel, PMID16868089](https://pubmed.ncbi.nlm.nih.gov/16868089/) | Brief exposures at nontoxic levels stimulated loss of heterozygosity at all 53 scored loci in normal diploid mouse ES cells, 1–8 × 10⁻³ per cell, 10–123 times untreated, comparable to Bloom-null cells. | Whether LOH induction is higher in BUBR1-hypomorphic cells at the same exposure. The figures are pooled across six exposures and 53 loci, so no single agent carries a frequency. The absence of persistent instability appears only in the title. |
| [Vitamin A diet, PMID40119798](https://pubmed.ncbi.nlm.nih.gov/40119798/); [CNP, PMID38136917](https://pubmed.ncbi.nlm.nih.gov/38136917/); [bFGF, PMID24014875](https://pubmed.ncbi.nlm.nih.gov/24014875/) | Three accessible or endogenous-sounding exposures with proliferation in the wrong direction: a retinol-enriched diet raised urothelial Ki-67 and the cell-cycle and PI3K-Akt pathways with no change in differentiation markers; CNP pre-treatment raised goat oocyte maturation with higher spindle-assembly and SAC gene expression; bFGF reconstitution normalised patient MSC proliferation, collagen and haematopoietic support ex vivo. | Retinol dose against paediatric tolerable upper intake limits; whether any CNP/NPR-B agonist changes somatic mis-segregation at growth-promoting exposures; an FGF delivery route without systemic FGFR exposure. None of the three carries a muscle, BUBR1 or fidelity endpoint, and each is easier for a family or prescriber to reach than a tool compound. |
| [BAL0891, PMID39184047](https://pubmed.ncbi.nlm.nih.gov/39184047/); [MPS1-IN-3, PMID23940287](https://pubmed.ncbi.nlm.nih.gov/23940287/); [reversine, PMID20624901](https://pubmed.ncbi.nlm.nih.gov/20624901/); [hesperadin, PMID12707311](https://pubmed.ncbi.nlm.nih.gov/12707311/); [GSK923295, PMID20167803](https://pubmed.ncbi.nlm.nih.gov/20167803/); [BAY1816032, PMID42014906](https://pubmed.ncbi.nlm.nih.gov/42014906/) | Six named administered agents whose stated cellular mechanism is checkpoint override, alignment failure or blocked error correction: a dual TTK/PLK1 inhibitor built to force aberrant mitotic exit faster than a TTK-specific agent; checkpoint override with increased aneuploidy in combination; ejection of Mad1 and ROD-ZWILCH-ZW10 from unattached kinetochores; anaphase entry with monooriented chromosomes; alignment failure with arrest; and loss of H2A-T120 phosphorylation. | Nothing that would convert them. Reversine's own abstract argues that error correction and the checkpoint are co-regulated and share catalytic machinery, so an MPS1 inhibitor is not expected to separate checkpoint suppression from blocked correction. These are recorded as direction, not as leads. |
| [Dapagliflozin, PMID39126109](https://pubmed.ncbi.nlm.nih.gov/39126109/) | In hypoxia-cultured H9c2 ventricular cells, dapagliflozin raised SIRT1 and SIRT3 expression and conferred resistance to hypoxic injury, stated three times and surviving the string damage. | SIRT2 — the only isoform the BUBR1 argument rests on — is never measured, and neither is BUBR1. Obligate glycosuria is deliberate caloric wasting in a child with failure to thrive, and no paediatric exposure data appears. Recorded so a later pass cannot resurrect it on the sirtuin resemblance. |

## What the completed census closes off

Five boundaries follow, and they are boundaries rather than shrugs: each is a
class that has now been searched to the end of this corpus and returned nothing.

**Machinery-level degradation blockade is closed as a supply strategy.** Five
probes — bortezomib, MG132, proTAME, apcin and MLN4924 — raise a mitotic protein
by disabling the proteasome, the ligase or the APC/C, and all five produce a
held cell rather than a corrected one. The consolidation is the change: these
should be read as one closed boundary, not five open probes.

**Direct mitotic-kinase and kinetochore pharmacology does not open a third
route.** Every named agent in this tranche that reaches a chromosome-fidelity
endpoint moves fidelity the wrong way, by design — TTK, Aurora, PLK1, Haspin,
CENP-E, Hec1, BUB1 and MPS1 agents alike are built to break the checkpoint this
patient has already lost. The single exception is UMK57, and it loses its own
benefit in three days.

**Aneuploidy-selective clearance runs against this patient by construction.**
The one record in the census that tests a clearance agent in primary,
non-transformed aneuploid cells kills them, and the killing is p53-mediated.
Selectivity for aneuploid cells is selectivity against the compartment this
child's tissue mass may partly depend on. Two overstatements are refused here:
the record supplies no p53-deficient comparator and no outgrowth experiment, so
whether that exposure would spare and thereby enrich p53-incompetent aneuploid
cells is a hazard to resolve and not a demonstrated selection pressure; and
engineered single-chromosome trisomic MEFs are not MVA tissue.

**BUBR1 protein supply was never measured by anyone.** No administered agent
anywhere in the outer corpus raises BUBR1 protein in a deficient cell, in either
direction. Niacin's mechanism is therefore unreplicated and unchallenged, and
the missing experiment is now demonstrably missing rather than merely unread.

**Muscle output is absent, not contradicted.** No verified lead measures force,
mass or motor function; the gap that apitegromab's rationale turns on — absolute
force in a growing animal — is untouched, and this corpus could not have closed
it because it did not contain it.

What stays open is equally specific: 323 flagged source questions, 143 outer
records with no saved abstract text, full-text depth for everything read at
abstract depth, the inner gene-query corpus, and the entire world outside these
retrievals.

## Primary qualification adds one entry and revises four

One record was taken to full text. The ledger goes to **109 route/compound
appraisals**, previously 108; the new agents that surfaced in this tranche — CNP,
CX-4945, MLN4924, BAY-320, Wy-14643, vismodegib, bFGF, zoledronic acid,
Nutlin-3a, the genotoxic panel — enter as constraints in the record-level
assessments and in the table above, not as routes.

| Ledger entry | What the record establishes | Effect on the nomination |
|---|---|---|
| `new-umk57` (new) | [Orr 2016, PMID27829147](https://pubmed.ncbi.nlm.nih.gov/27829147/) / [PMC5130158](https://pmc.ncbi.nlm.nih.gov/articles/PMC5130158/), read in full. 100 nM UMK57 potentiates MCAK, reduces metaphase k-MT attachment stability by more than 35%, significantly lowers anaphase lagging chromosomes in U2OS, HeLa and SW-620 and chromosome non-disjunction in HCT116 — and the gain is gone within 72 hr with target engagement intact. | Retained as the ledger's clearest case that chemically potentiating a protein does not inherit the durable benefit of overexpressing it, and as a duration and cross-resistance constraint on every attachment-balancing route. Not a repurposing candidate: zero matches in the 14,882-entity clinical and approval catalog. Does not displace niacin, apitegromab or the everolimus probe. |
| `new-aurora-balancing` (revised) | The same source runs ZM447439 as a timing experiment rather than a single-timepoint balancing arm. At 250 nM, partial inhibition verified on early-prometaphase histone H3 Ser10 phosphorylation, it neither raised lagging-chromosome rates nor blocked UMK57's effect in U2OS within 1 hr; it prevented the 72 hr rebound in U2OS, delayed resistance in SW-620, did not do so in HeLa, and was overrun by 96 hr. Clones grown continuously under UMK57 became specifically less sensitive to ZM447439, with no change to doxorubicin, taxol, vinblastine or BI-2536. | Qualification, not counter-evidence; disposition unchanged. It bounds the entry's own timed-treatment caveat on three axes — combination-dependent, line-dependent, durability under four days — and adds the entry's only endpoint scored in cycling anaphase rather than under MG132. [PMID26188358](https://pubmed.ncbi.nlm.nih.gov/26188358/) separately names ZM447439 as a polyploidy-inducer whose polyploid product resists apoptosis through Bcl-xL. Both sources added. |
| `new-apcin` (revised) | [PMID32152539](https://pubmed.ncbi.nlm.nih.gov/32152539/): apcin's sign is set by regulatory context, not by dose — net APC/C inhibition prolonging mitosis at low SAC activity, net activation shortening it at high SAC activity, from a common Cdc20 site, cooperating with p31comet to relieve MCC-dependent inhibition. | Disposition unchanged; the entry can no longer be carried as a monotonic APC/C inhibitor. Net activation shortens mitosis, which is the direction of the premature exit this lesion already produces, so no apcin exposure has a knowable direction here. The inversion is **not** transferred to `new-protame`: the record concerns a Cdc20 site targeted by apcin and never mentions proTAME. |
| `prior-27` (revised) | The entry declines "AICAR reference experiments" without ever citing them. [PMID21315436](https://pubmed.ncbi.nlm.nih.gov/21315436/) is that work: the killing is p53-mediated apoptosis in primary MEFs trisomic for chromosome 1, 13, 16 or 19, with AICAR and 17-AAG especially efficacious in combination against aneuploid human cancer lines. | Disposition unchanged — the entry already reads decline lead, and nothing here raises a decline to an exclusion. What is new to the ledger is that the demonstrated aneuploidy selectivity is conditioned on an intact p53 response, that killing the vulnerable compartment is a tissue-supply harm this entry did not carry, and that engineered trisomic MEFs are not MVA tissue. Source added. |
| `prior-23` (revised) | The same AML record names ABT-263 directly: added to Oxindole-1 or to ZM447439 it produced synergistic loss of viability with greater sustained tumour-growth inhibition in AML lines and primary blasts, and Bcl-xL knockdown but not Bcl-2 reproduced the synergy. | Records the hazard against the senotherapeutic reserve rather than changing it: the exposure is selective for cells that survive an abnormal mitosis, which is the compartment this child may live on. The model is drug-induced polyploidy in transformed leukaemia, not constitutional mosaic aneuploidy, so the dependency is untested here rather than transferred. |

Four entries were examined and deliberately **not** revised. `prior-14` already
cites the same AICAR/17-AAG source, already declines the lead, and already
carries the tissue-supply harm verbatim; the p53 statement is AICAR-specific and
does not attach to the HSP90 entry. `prior-26` already declines on
chloroquine-specific label myopathy, which is more decisive for this patient
than anything in the record, and the record describes no chloroquine experiment
at all. `new-protame` is unchanged and its boundary is reaffirmed. `prior-01` is
unchanged: UMK57 does not compete with niacin, because it acts on attachment
dynamics and not on BUBR1 amount.

### UMK57: what it establishes, and why it is not a candidate

[UMK57](https://pmc.ncbi.nlm.nih.gov/articles/PMC5130158/) was used at 100 nM for
every assay. The paper's stated reason is a potency floor — the lowest
concentration giving maximal suppression of lagging chromosomes — and it
separately notes that this is about five-fold below the calculated U2OS IC50, a
value never given; 100–2000 nM was used for proliferation assays only, and
controls were 0.1% DMSO. Specificity is unusually well built. No kinesin ATPase
tested is inhibited, MCAK-dependent microtubule depolymerisation is enhanced in
two independent in vitro assays, ultracentrifugation sedimentation and
microscopy, and the one-group analog UMK95 is inert on depolymerisation,
proliferation and segregation. At the working dose the compound does not alter
total MCAK levels, MCAK localisation or spindle organisation at any mitotic
stage, so potentiation of activity without change of amount or position is
measured rather than inferred. Mis-segregation falls in control and
Kif2b-depleted cells and does not change in MCAK-depleted cells.

The endpoints are segregation outcomes rather than arrest biochemistry. Anaphase
lagging-chromosome rates fall significantly in U2OS, HeLa and SW-620 and are
unaffected in non-transformed diploid RPE-1 and BJ; because the paper notes that
many lagging chromosomes still reach the correct daughter, chromosome
non-disjunction was measured in HCT116 released from a 15 h 100 µM monastrol
arrest, with mitotic shake-off and 8 h of replating to complete cytokinesis
before scoring. Only the second is a completed division, and neither counts
absolute faithful daughters per starting cell. The only effect size stated
numerically anywhere in the text is a greater than 35% reduction in metaphase
k-MT attachment stability; prometaphase stability is unchanged acutely, and
every photoactivation measurement was made under 5 µM MG132 to prevent mitotic
exit, so house rule 4 bounds those numbers although not the anaphase and
non-disjunction endpoints. That is what makes the record procedurally new: the
attachment-balancing entries the ledger already holds, `new-aurora-balancing`
and `new-reversine-balancing`, buy their improvement under MG132 exit blockade,
and this is the first of that class whose benefit is scored in cells that
finished.

Then it disappears, and the disappearance is the transferable finding. Lagging
chromosome rates return to control levels within 72 hr in U2OS, and a
time-dependent decrease in responsiveness was also observed in multiple CIN
cancer lines with no timing stated for those lines. The paper closes the trivial
explanations one at a time. Serially transferred culture medium retains UMK57
activity to 120 hr and renewing the medium every 24 hr does not prevent the
rebound, so it is not compound decay. Verapamil at 6 µM, positive-controlled in
the same cells by enhancing the doxorubicin response, did not significantly
alter the loss of responsiveness, from which the authors conclude that
resistance occurs independently of multidrug efflux pump activity — one
inhibitor, so efflux is made unlikely rather than excluded. Metaphase
destabilisation is equivalent at 1 hr and 72–120 hr and centromeric
phospho-MCAK is unchanged between 1 hr and greater than 72 hr, so target
engagement persists; note that it is verified to 120 hr and not across the
colony assay, since no metaphase value is reported at 6–8 weeks. Under extended
treatment the stable versus unstable microtubule fractions and the
unstable-population turnover are also unchanged, which is what confines the
adaptation to prometaphase k-MT stability.

The long-horizon endpoint is negative, and it is the one that matters
clinically. Continuous UMK57 through 6–8 weeks of single-cell colony growth did
not reduce the fraction of U2OS cells deviating from the modal chromosome copy
number, scored by FISH with centromeric α-satellite probes for chromosomes 2 and
3 in two independent colonies per condition; RPE-1 colonies grown 4–6 weeks stay
stably diploid and insensitive. A per-mitosis fidelity gain did not become a
better karyotype. What returns instead is the same defect relocated to a
different mitotic stage: prometaphase k-MT attachment stability increases
significantly at 72–120 hr and at 6–8 weeks and is reversible on UMK57 removal,
active Aurora B at centromeres falls substantially after 72 hr, and substrate
phosphorylation does not follow it — phospho-KNL1 rises specifically at
unaligned kinetochores while phospho-Hec1 is unchanged. The authors read this as
an alteration of the Aurora B signalling network, the spatially regulated
activities of both the kinase and the counteracting phosphatases, rather than as
loss of kinase output. For `new-aurora-balancing` that has a direct consequence:
centromeric active-kinase staining and substrate phosphorylation dissociate in
the same cells, so that entry's standing requirement for checkpoint surveillance
cannot be discharged by centromeric kinase staining.

The overexpression-versus-potentiation distinction is the entry's headline and
must carry its qualifier. The authors state that the adaptive response was not
observed previously when CIN was suppressed by MCAK or Kif2b overexpression, and
name the distinction between protein overexpression and chemical enhancement of
protein activity. That is a comparison against prior published work; no
overexpression arm was run through this paper's FISH assay, so it is a
cross-study inference rather than a controlled contrast. It remains the ledger's
clearest named case for the distinction, which is house rule 2 read from the
drug side: any restoration route whose supporting evidence is a genetic
gain-of-function must state that the chemical implementation is a separate claim
and must measure durability over weeks rather than at a single timepoint. Inside
one experiment the same asymmetry runs the other way — transient GFP-MCAK makes
cells hypersensitive at the same 100 nM dose, with a substantial fraction
showing chromosome-alignment and spindle-organisation defects while GFP and
GFP-Kif2b controls show very few — so the response scales with the amount of the
potentiated protein, not with dose alone. The authors' own closing reading of
the resistance is that cells achieve it through multiple pathways, which they
take as strong selective pressure to maintain high mis-segregation in CIN cancer
cells; blocking the Aurora B route therefore does not bound the problem.

For this patient the agent fails on three separate points, each stated at the
strength the source supports and no higher. **BUBR1.** The only BubR1
measurement in the paper is kinetochore BubR1 level at the working dose, cited
together with Aurora B and Astrin localisation as the basis for the conclusion
that UMK57 does not affect chromosome bi-orientation or the spindle assembly
checkpoint; an unchanged kinetochore pool is the implied reading and no value is
stated in the text. Total BUBR1 abundance was never measured and no
BUBR1-deficient cell was tested, so the record is silent on BUBR1 supply rather
than affirmative against it. This is an ordinary null, not a stronger one.
**Axis.** The defect UMK57 corrects, hyper-stable k-MT attachments in CIN cancer
cells, sits on a different axis from low BUBR1 abundance with a weak SAC, not at
the opposite end of the same one; whether BUBR1-hypomorphic cells have hyper- or
hypo-stable attachments is not addressed anywhere in the paper, which is exactly
why the direction in the deficient cell is untested. **Direction.** Two
harm-direction signals are on record and untested in a checkpoint-weak
background: the GFP-MCAK hypersensitivity above, and a subtle but significant
reduction in the inter-kinetochore distance of aligned kinetochore pairs, with
no value stated. That tension drop occurred without significantly altering
mitotic progression, so in a checkpoint-competent cell it was not large enough
to delay mitosis; whether the same drop matters in a SAC-weak cell is untested
in either direction, and the paper states no safeguard.

Three scope limits complete the record. The work is cell culture only — three
CIN cancer lines, a fourth cancer line used for the non-disjunction assay, and
two non-transformed diploid lines — with no animal, no muscle, no
BUBR1-deficient and no non-transformed CIN model. UMK57 is not inert in
non-transformed cells even though its segregation effect there is null: RPE-1
and U2OS show differential sensitivity to UMK57 but
not to vinblastine, which separates it from a direct tubulin agent, and neither
the readout, the more-sensitive line nor the magnitude is given. And the paper
twice calls the fidelity-experiment doses sub-lethal while reporting that UMK57
inhibits proliferation in a dose-dependent manner over the separate 100–2000 nM
proliferation range, with no effect size at 100 nM and no cell line stated; for
a child with recorded rhabdomyosarcoma and failure to thrive the difference
between a sub-lethal working dose and a dose-dependent antiproliferative range
is not cosmetic. All figure panels cited here are named by the running text but
were not viewed, so no panel value is credited and every magnitude other than
the greater-than-35% metaphase result is figure-resident and unstated.

Eligibility is what actually decides the disposition, and it should be visible
as a separate question from data quality. UMK57 has zero matches in the saved
14,882-entity clinical and approval catalog, and before this tranche it appeared nowhere in this repository; its only other
trace is its own record in this tranche’s review file. It is a
research-tool compound with no human exposure in the reviewed record, which is
the same boundary that holds proTAME and Apcin as probes rather than candidates.
It is also the strongest case for that boundary being about eligibility rather
than about experimental quality: proTAME was declined partly because its readout
was arrest biochemistry with preserved securin and cyclin B1, whereas UMK57 has
MCAK-dependence, an inactive analog control and endpoints scored in anaphase and
after cytokinesis — and is still declined, for absence of any clinical
existence, absence of any demonstrated segregation effect in a non-cancer cell,
and mismatch to the BUBR1 lesion.

Three things would change the call, and the first is decisive. A MCAK
potentiator applied to a non-transformed BUBR1-hypomorphic cell — ideally MVA
patient fibroblasts carrying a truncating allele in trans with a kinase-domain
missense allele — scoring anaphase lagging chromosomes, non-disjunction and
absolute faithful daughters per starting cell, with a proficient-genotype arm
and a washout arm. A benefit there, with a direction and an amount, reopens the
route; a harm there converts this entry from a duration constraint into an
active contraindication for the whole attachment-destabilising class. Second,
any evidence that the fidelity gain survives past 72 hr in any cell: the rebound
is adaptive and reversible on washout with target engagement preserved, so an
intermittent schedule that restores responsiveness and holds fidelity across
weeks, measured by the same FISH modal-deviation assay with more than two
colonies per condition, is the specific missing experiment. Chronic exposure in
a growing child is the only schedule that matters clinically, and it is exactly
the one shown to fail. Third, eligibility: a catalog-eligible MCAK potentiator
with human exposure data, a paediatric-relevant window and a muscle-output
endpoint. Absent the third, even a positive first and second enter as a
mechanism route and a follow-up pointer, not as a repurposing candidate.

## Reading and source-quality boundaries

The [fixed selection](outer-unfiltered-final-selection.json) takes **all**
remaining available outer abstracts in the same SHA-256 order under seed
`mva-outer-unfiltered-v1`, excluding the full 1,460-record gene-query set and
the 2,161 previous outer readings, so the selection rule is "all of them in hash
order" rather than a head count. No drug name, outcome, publication date,
publication type or biological filter selected the batch; the queue resumes
strictly after the last record of the preceding tranche; each of the 57
assignments was disjoint; provenance fields were preserved and raw abstracts
remain outside Git. Of the 3,741 records that remained, 3,598 carry a saved
abstract string and 143 do not.

`verify_final.py` replays that selection, the exclusions and the content hashes
against the original retrievals. It **PASSED with full source replay**,
`source_selection_replayed: true`, so every selected record's saved abstract
hash was reproduced from the original corpus rather than taken from the review
file, and the 3,598 reviewed records are exactly the records that remained. It
does not prove biological efficacy, experimental reproducibility or independent
study counts.

Three source-quality items constrain what may be quoted. Systematic
word-dropping is present in a subset of saved abstract strings, with terminal
nouns and numerals elided mid-clause; readers were required to declare it per
record rather than reconstruct missing text, and no omitted value was inferred.
Most of that loss is in the reading channel rather than in the saved string.
2,073 of the 3,598 filed `source_abstract_quality` fields assert dropping,
elision or truncation, yet every saved string opened during verification was
grammatical, complete and hash-matched to `abstract_sha256`. The
lead-verification pass demonstrated the artifact in at least six records —
UMK57 (3368), GSK923295 (2485), TC-A2317 (2471), Nutlin-3a (2352), YLT-11
(4490) and AR-12 (5184) — and on YLT-11 two consecutive runs of the same reader
tool dropped different tokens. Any elision note in any chunk should therefore be
re-checked against the input JSON before it is trusted, and two examples are
already corrected above. Concretely: the UMK57 record's saved
string is complete at 1,121 characters, and the AML polyploidy record — filed
with a house-rule-7 caveat that the Bcl-2-family partner was elided — in fact
names **ABT-263** in full on re-reading, which is why `prior-23` rather than an
unnamed agent carries that hazard. Preprint-source records remain in scope, and
reading an abstract does not establish peer review or final publication
identity. Records whose saved string is a Significance section, a conference
abstract or is truncated upstream are flagged individually in
`source_abstract_quality` and are not treated as complete reports. The
corpus-wide counts quoted above — no myostatin mention and one GDF8 mention,
two BubR1/BUB1B mentions, one niacin mention, one mosaic-variegated record, 55
records mentioning muscle — were computed over the saved strings with
case-insensitive regex. The single mosaic-variegated record is
[PMID36333305](https://pubmed.ncbi.nlm.nih.gov/36333305/), which names mosaic
variegated hyperploidy in SLF2/SMC5-mutant patient cells (Atelis syndrome); it
is an observational genetics report with no administered agent. So the corpus
contains one constitutional chromosomal-instability syndrome and no drug tested
in one.

The complete [record-level assessments](outer-unfiltered-final-review.json)
retain genetic rescues, negative observations, exposure hazards and
computational nominations alongside the leads, in input order, one annotation
per abstract. The lead-verification and triage files are held with the raw
tranche material outside Git. No unnamed agent is assigned to an approved
medicine, and the table above holds the original abstract-stage questions;
primary qualification is recorded separately.

Nothing in this tranche displaces the niacin/apitegromab comparison, and with
the outer abstracts complete that sentence now has a wider scope than it did:
not "nothing in these 512" but nothing in the saved outer corpus. No record
supplies an administered medicine that restored faithful division in a
BUBR1-limited cell; none supplies an absolute force gain in a growing animal;
and the only administered agent in the ledger whose fidelity benefit is scored
after mitosis rather than under an exit blockade is a research-tool compound
with no clinical existence, whose benefit expires in three days and whose own
resistance mechanism is a rewiring of the Aurora B node another ledger route
depends on. Unanswered scientific and coverage questions remain — 323 flagged
source questions and 143 records with no saved abstract text among them.
Completion of the census does not trigger publication or Track 2 submission.
