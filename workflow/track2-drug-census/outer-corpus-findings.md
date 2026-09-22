# Bounded outer-corpus review: 60 additional complete abstracts

**Scope:** 60/60 explicitly selected saved abstracts read on 21 September 2026; no full-paper claim. The broad corpus contains 7,362 records. Subtracting the exact source/ID union of 1,460 gene-title/abstract records leaves 5,902, of which 5,759 have saved abstracts. This pass reviews only 60 of those 5,902. It makes no efficacy or irrelevance claim about the remaining records.

The 60 records have no source/ID overlap with the authored gene-abstract notes, recovered-missing-abstract notes, or explicitly appraised PubMed-supplement identifiers. None of their PMIDs appeared in the inspected census Markdown narratives or decisions at commit `103980903eef11a267c728148f15cef215790aad`. This is local prior-appraisal accounting, not a claim that a route is scientifically novel. The unrestricted corpus SHA256 is `e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9`.

Selection used muscle/progeroid/fidelity/aging terms together with benefit/intervention terms, plus a separate beneficial-title scan and selected mechanistic/harm comparators. No clinical-catalog alias was required. The discovery ranking contains 419 candidates; those candidates were **not** all read. The controlling manual ledger, exact patterns, score formula, selected IDs, hashes, and review-depth information are in the accompanying JSON files. Six batches of ten were read and persisted. Notes are original judgments rather than reproduced abstracts.

**Primary-review decision:** TM5484/PAI-1 and hesperetin/CISD2 are substantive
experimental muscle-compensation alternatives. TM5484 has a direct contraction
endpoint; hesperetin has a small motor-performance result. Their primary
qualifications below explain why they merit comparison without displacing the
current niacin/apitegromab nomination on evidence of MVA benefit.

## Most actionable additions

| Route and exact MED ID | What the complete saved abstract establishes | Competitive implication and next primary check |
|---|---|---|
| **Cortical-tension modulation: ConA**, [35399525](https://pubmed.ncbi.nlm.nih.gov/35399525/) | ConA at 10 µg/mL during culture of vitrified mouse oocytes improves kinetochore–microtubule attachment, SAC activity, polar-body extrusion and embryo development, while reducing aneuploidy. | A particularly useful new fidelity comparator: actual administered-agent, chromosome and developmental outcomes. Confirm denominators, exposure window, completed segregation and harm boundaries. Vitrification injury and oocyte mechanics differ from BUBR1-deficient muscle; this is not an existing-medicine nomination. |
| **Chaetocin / nuclear-envelope budding**, [39352925](https://pubmed.ncbi.nlm.nih.gov/39352925/) | Chaetocin limits progerin-associated envelope budding/chromatin loss, ameliorates HGPS cell/mouse defects and extends HGPS mouse lifespan; sustained ERK signaling is implicated. | A concrete pharmacological progeroid-rescue route worth primary qualification. Verify dose, lifespan magnitude, functional tissues and genomic/toxicity outcomes. Mechanism targets progerin pathology, not demonstrated BUBR1 insufficiency. |
| **Hesperetin / CISD2**, [35871686](https://pubmed.ncbi.nlm.nih.gov/35871686/), [38263133](https://pubmed.ncbi.nlm.nih.gov/38263133/) | Late-life oral treatment improves aged-mouse healthspan, body composition and metabolic/organ phenotypes with Cisd2 dependence. A separate study gives older-donor keratinocyte and five-month late-life mouse-skin benefit. | Stronger organismal/causal evidence than an expression-only lead. Determine exact survival and muscle-strength endpoints, dose/exposure and independence of cohorts. No BUBR1 correction is established by these abstracts. Hesperetin must not be silently substituted with hesperidin. |
| **HDAC6 inhibition: tubacin/tubastatin A**, [39834092](https://pubmed.ncbi.nlm.nih.gov/39834092/) | Both agents improve spindle morphology and polar-body extrusion in TUBB8-D417N oocytes, after a defined pathogenic-variant perturbation. | A functional compensation precedent with an actual mutation, but of another gene. Test whether primary data include euploid embryos/offspring rather than morphology alone. [33195286](https://pubmed.ncbi.nlm.nih.gov/33195286/) additionally reports SIRT2/HDAC6-inhibition rescue after RASSF1A depletion, without naming the inhibitors. |
| **Indole acetic acid / AhR–Sirt2**, [40197001](https://pubmed.ncbi.nlm.nih.gov/40197001/) | Supplementation extends fly lifespan and stress resistance; AhR/Sirt2 mutant controls attenuate the effect. | A new upstream Sirt2-linked administered route. Primary mammalian transfer, exposure and BUBR1/fidelity evidence are needed; fly lifespan is not mammalian checkpoint restoration. |
| **Liraglutide**, [38331723](https://pubmed.ncbi.nlm.nih.gov/38331723/) | Improves viability and muscle protein-catabolic markers in high-glucose C2C12 myoblasts, with AMPK-associated/transcriptomic changes. | Direct muscle-cell benefit warrants a reserve entry. Determine dose, myotube/contractile outcomes, and whether benefit persists outside hyperglycemia; the abstract gives no BUBR1 mediation. |
| **Calcitriol**, [31880094](https://pubmed.ncbi.nlm.nih.gov/31880094/) | Exogenous 1,25-dihydroxyvitamin D3 improves deficient/naturally aged mouse bone phenotypes and passaged human stromal-cell osteogenesis through VDR–Ezh2–p16 signaling. | Defined administered compound with organismal and human-cell outcomes, but bone-focused. Do not substitute ordinary vitamin D supplementation or infer muscle/fidelity efficacy. |

None should be rejected solely for lacking the actual MVA alleles, since the current nominations share that gap. On the evidence read here, none yet establishes superior expected muscle benefit to the niacin-restoration/apitegromab-compensation comparison. ConA offers the clearest additional fidelity endpoint; chaetocin and hesperetin offer the strongest new organismal-aging signals. These remain different experimental questions, not interchangeable scores.

## Particularly informative full-text blind spots

Several abstracts describe a potentially useful intervention but leave its implementation unnamed. They deserve directed primary-source identification rather than exclusion by catalog matching:

- **PAI-1**, [24778222](https://pubmed.ncbi.nlm.nih.gov/24778222/): the saved abstract reports deficiency-mediated organ/survival rescue in Klotho mice. It does not identify a drug. This review's genetic classification intentionally remains at abstract depth; any subsequently verified TM5441 drug arm must be recorded as separate primary evidence.
- **p53 plus LSD1 inhibition**, [39809271](https://pubmed.ncbi.nlm.nih.gov/39809271/): reactivates aged-human-cell centromeres/CENP-A loading and reports mitotic rejuvenation; chemical versus genetic implementation, faithful division and tumor-suppressor consequences require checking.
- **CDK1 inhibition**, [35910534](https://pubmed.ncbi.nlm.nih.gov/35910534/): unnamed inhibitor and siRNA accompany proliferation/senescence benefit in old tonsil-derived stromal cells. Fidelity is an essential missing outcome.
- **Ferroptosis inhibition**, [39910323](https://pubmed.ncbi.nlm.nih.gov/39910323/): unnamed inhibitors improve aged-oocyte quality with Fth1 genetic corroboration. Named RSL3/FAC are the damaging induction agents, not the proposed treatment.
- **Five-compound hepatocyte reprogramming**, [38144457](https://pubmed.ncbi.nlm.nih.gov/38144457/): systemic 5C improves liver regeneration/function, but the abstract does not name its ingredients. Avoid converting a mixture result to an individual-drug claim.

Lower-priority administered observations include nectandrin B in old fibroblasts (31199782), fucoxanthin during oxidative injury of placenta-derived stromal cells (40999517), and DZNep in inflamed chondrocytes (28744016). Their abstracts supply cytoprotection or senescence/inflammation outcomes, without chromosome or muscle-function rescue. Existing navitoclax/clearance and AICAR routes gain additional contextual evidence rather than becoming newly discovered routes.

## Primary follow-up: PAI-1 produces a substantive muscle challenger

**TM5484 merits an experimental muscle-function comparison.** The primary
manuscript reports stronger nerve-evoked gastrocnemius twitch/tetanic force and
preserved fiber area after dietary treatment from six to twelve months
(nine mice per group). Chow contained **0.01% drug, or 100 mg/kg feed**;
this is not a 100 mg/kg bodyweight/day regimen. A later-start arm
(18→19 months, ten mice/group, 10 mg/kg/day oral gavage) improved fiber area,
but twitch was unchanged and tetanic force showed only a trend. Body/muscle wet
weight did not significantly improve in the prevention comparison. Thus a
prevention result should not be advertised as reversal of established weakness.

Complete primary text/captions were read, but original figure images, numerical
force effects and force-normalization units were unavailable. Specific-force
improvement is unverified. The later-arm caption conflicts with its methods on
vehicle/duration; the regimen above follows Methods. Proposed IGF1/plasmin
mechanisms are not demonstrated necessities for the muscle effect.
[TM5484 primary study](https://pubmed.ncbi.nlm.nih.gov/33213843/).

Two supporting studies test different interventions. **TM5441** delayed
mortality in Klotho-deficient mice (n=11; 100 mg/kg/day in chow); the quantified
genetic survival gains and open-field/organ findings belong to genetic PAI-1
loss. The drug median and starting age were not recovered. A later PAI-1-knockout
study preserved muscle/grip outcomes in aged females, without demonstrated male
benefit; it administered no therapeutic drug. These are not replications of
TM5484 or evidence of its sex-specific response.
[TM5441 primary study and Figure 1](https://www.mmt.med.tohoku.ac.jp/en/pdf/24.pdf),
[genetic muscle study](https://www.aging-us.com/article/206318/text).

PAI-1 compensation now has more than a biomarker rationale. Its direct mouse
force endpoint gives it a meaningful advantage over several other experimental
reserves, while human efficacy/exposure, approved-product eligibility and
transfer to developmental MVA muscle remain unestablished. None of these
studies treats a BUB1B model. **Retain niacin/apitegromab as the nominated pair
and elevate TM5484 for subsequent functional comparison**, using delayed-start
as well as preventive conditions, absolute and normalized force, cell renewal
and chromosome fidelity. The same actual-allele uncertainty applies to all
three candidates; it is not a selective exclusion rule.

## Primary follow-up: hesperetin has motor benefit, with a small denominator

Original main and supplementary figures confirm a functional result beyond
body composition. Male C57BL/6 mice received dietary hesperetin, nominal
100 mg/kg/day. Rotarod performance improved after six months from age 20 months
(**four mice per arm**); the three-month comparison was nonsignificant. This
measures composite motor performance, not isolated muscle force. Femoris
histology improved; gastrocnemius histology was nonsignificant despite stronger
wording in the main text. Lean percentage increased, without a marked
significant old-group difference in absolute lean volume.

A separate muscle/heart Cisd2 knockout model failed to improve despite
comparable exposure, supporting a context-specific requirement. In a separate survival schedule starting
at 21 months, median lifespan was 28.2 versus 25.95 months for vehicle
(n=19 versus 8); survival and rotarod schedules/sample sizes differ. No BUB1B
model was tested.
**Elevate an experimental muscle-function comparator**, while retaining the
male-only, small-cohort and composite-endpoint limits. Human functional evidence
and a relevant MVA-muscle bridge remain unestablished.
[Primary paper and supplement](https://link.springer.com/article/10.1186/s12929-022-00838-7).

## Primary follow-up: chaetocin and an important model boundary

The chaetocin primary Results, Methods and Figure 6 caption were read from
accessible full-text XML. In HGPS fibroblasts, 50 nM for 24 hours moved progerin
away from the nuclear envelope without lowering its abundance. In
Lmna-G609G homozygous mice, 0.1 mg/kg intraperitoneally every three days from
four weeks of age increased median survival **140→160 days versus vehicle**
(+14.3%). The untreated median was 95 days; using that as the drug comparator
would inflate the treatment-specific effect. Survival-group sizes and original
figure images were not independently verified here. Vascular/tissue and survival
benefit does not establish skeletal-muscle contractility or BUBR1 restoration.
The paper discusses UPCP as a separate progerin–BUBR1 intervention; chaetocin
must not inherit that mechanism. **Retain a progerin-specific experimental
precedent, not a replacement for either nominated MVA route.**
[Primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11474064/).

Another full-text check resolves an easy model-transfer error. The MCP-1 paper
includes BubR1-hypomorphic mice for a biomarker comparison, but its
dasatinib/quercetin intervention is in **Ercc1-deficient mice**, and its
rapamycin intervention is in old wild-type mice. The reported biomarker response
therefore supplies no drug-mediated rescue of the BubR1 mice. Its D+Q serum
analysis reuses a prior treatment cohort; it is not an independent replication
of that treatment's lifespan result.
[Primary Results and Figure 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC5847863/).

## Counterevidence that affects interpretation

- Muscle relevance can cut against indiscriminate checkpoint suppression. [25960061](https://pubmed.ncbi.nlm.nih.gov/25960061/) links Mps1/Mad2 disruption to aneuploid, poorly differentiating satellite cells and protective p21 arrest. [38806674](https://pubmed.ncbi.nlm.nih.gov/38806674/) shows that excessive checkpoint activation also damages gut and marrow. The PLK1–PP2A balance paper (37712330) supports testing a functional window, not maximization of arrest or attachment strength alone.
- A lower aneuploid fraction can reflect killing: p53/p21 suppression during arsenite exposure (20000476) reduces centromere-positive micronuclei in a setting where the authors invoke mitotic-arrest apoptosis. DL-PDMP kills aneuploid primary fibroblasts (28775166). CCB02/PJ34 improve cystic-kidney function through selective apoptosis (38385746), not repair of the abnormal cells.
- GSK3 inhibitors cause chromosome alignment failure/nondisjunction with RNAi corroboration (17697341). Conversely, AICAR/A769662 restore spindle orientation in PRL-overexpressing cells (36336348); that does not yet establish corrected aneuploidy or BUBR1 deficiency.
- SIRT2 direction depends on compartment, isoform and lesion: hippocampal SIRT2.3 lacks clear benefit (33803627), low Sirt2 participates in insect diapause (29769432), and SIRT2 inhibition rescues RASSF1A-deficient oocyte morphology (33195286). SIRT2-mediated sensory-neuron protection against cisplatin is disrupted by spironolactone (32134743); spironolactone is the adverse counterintervention.

## Deliverables and coverage limits

`outer-corpus-review.json` contains all 60 original, source-hashed assessments: 31 primarily administered-agent records, 17 genetic, five basic-mechanism, four implementation-unspecified/other, two reviews and one biomarker-association record; 27 receive a primary-follow-up flag. These are record classifications, not 31 beneficial drugs or 27 novel routes. Agent arrays explicitly distinguish damaging challenges, engineered triggers and tracers.

`outer-corpus-selection.json` controls the bounded selection; `outer-corpus-candidate-ranking.json` preserves discovery-only metadata; `outer-corpus-overlap-audit.json` preserves reviewed-set overlap checks and source hashes. Every selected source/ID/title/hash and complete-saved-abstract flag was validated. These annotations preserve the abstract-only reading depth. The separately labeled primary follow-ups below do not retroactively upgrade every abstract to a full-paper review.
