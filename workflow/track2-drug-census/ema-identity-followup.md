# Three EMA name gaps: identity and evidence follow-up

**The development aliases do not close these three census gaps.** Official sources establish the named medicine identities, but **37 documented names, identifiers and separately labeled typography variants** produce **zero exact cross-references** to main preferred names/aliases, the DrugCentral unverified supplement, the FDA2026 supplement or Broad preferred names. They also produce **zero title/abstract matches across all 7,362 saved gene-corpus records** (7,191 with abstracts). This is a bounded missing-name finding, not certified worldwide unique-drug identity or proof of null efficacy.

| Medicine | Official identity recovery | Evidence relevant to the proposal |
|---|---|---|
| **Andembry / garadacimab** | CSL312, CSL-312; FDA UNII 32W6AJL0DY and named garadacimab-gxii subconcept. EU initial authorisation 2025-02-10. | Actual-drug randomized HAE attack prevention. No retrieved BUB1B-restoration or muscle-function experiment. |
| **Kostaive / zapomeran** | Formulated vaccine ARCT-154; active substance mRNA-2105. EU initial authorisation 2025-02-12. | COVID-19 clinical protection and antibody endpoints. Actual ARCT-154 mouse muscle study measures spike/replicase clearance, not muscle function. |
| **Waskyra / etuvetidigene autotemcel** | Telethon003, OTL-103, GSK2696275; FDA adds etu-cel, Telethon-003, TLT003. EU initial authorisation 2026-01-09. | WAS-gene HSPC product with hematopoietic, infection and bleeding outcomes in WAS; it supplies WAS, not BUB1B. |

Sources: [EMA Andembry identity](https://www.ema.europa.eu/en/medicines/human/orphan-designations/eu-3-21-2532), [FDA garadacimab substance record](https://precision.fda.gov/ginas/app/ui/substances/a858d516-cf0b-4cae-a854-6ec6f70e98cb), [EMA Kostaive assessment pp12,14](https://www.ema.europa.eu/en/documents/assessment-report/kostaive-epar-public-assessment-report_en.pdf), [EMA Waskyra assessment p12](https://www.ema.europa.eu/en/documents/assessment-report/waskyra-epar-public-assessment-report_en.pdf), [FDA Waskyra review p4](https://www.fda.gov/media/190281/download?attachment=).

The name recovery is specific: **ARCT-021 is a different predecessor; LUNAR-COV19/STARR describe platforms**. They were not admitted as product aliases. Trial identifiers were not treated as additional drugs. Source naming links do not independently prove formulation or active-moiety equivalence.

## Actual evidence and comparison

Garadacimab's VANGUARD publication reports 65 randomized HAE patients and 64 treated/analyzed (39 drug, 25 placebo); six-month HAE attack rates were 0.27 versus 2.01 per month. It offers clinical evidence for a different disease and endpoint, not BUB1B rescue. The primary abstract and complete registry were inspected. [PMID36868261](https://pubmed.ncbi.nlm.nih.gov/36868261/), [NCT04656418](https://clinicaltrials.gov/study/NCT04656418).

Kostaive's pivotal trial and booster trial concern vaccine protection, neutralizing antibodies and tolerability. Its actual-product mouse study measured protein clearance after intramuscular administration; the presence of a muscle tissue endpoint does not mean strength or compensation was tested. Replacing spike/replicase instructions with BUB1B would create another intervention. [Pivotal original trial](https://www.nature.com/articles/s41467-024-47905-1), [booster primary abstract](https://pubmed.ncbi.nlm.nih.gov/38141632/), [mouse clearance study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11999312/).

Waskyra's clinical program concerns WAS: the original 2019 single-arm interim report followed eight boys, and FDA integrates 27 patients across two single-arm trials and expanded access. The 2013 and 2019 reports overlap, so those cohorts must not be added. WASP engraftment, reduced infections and bleeding do not establish muscle compensation or BUB1B restoration. [Ferrua 2019 original full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6494976/), [FDA integrated review](https://www.fda.gov/media/190281/download?attachment=).

The broader search found adjacent evidence worth preserving: **donidalorsen** was associated with a CK reduction in an incidental HAE analysis, with no demonstrated between-arm muscle-function result; **3F7 anti-FXII/FXIIa antibody** improved obstructed mouse-kidney histology and tubular repair. Neither can be attributed to garadacimab or interpreted as MVA rescue. The donidalorsen paper explicitly lacks exercise/fitness data and longitudinal CK–fatigue linkage. [Donidalorsen CK study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12955664/), [3F7 kidney study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12775217/).

**Research disposition:** none of the three currently earns promotion over the separate niacin restoration hypothesis and apitegromab muscle-compensation comparator. The existing comparison has specific, though incomplete, anchors in NAD/SIRT2/BUBR1 biology and disease-associated motor-function evidence; the new medicines' retrieved evidence concerns HAE, COVID-19 or WAS. This is an evidence-priority judgment, not a claim that niacin/apitegromab are effective in MVA or that the three medicines have no effect. Actual BUB1B restoration, MVA benefit and muscle compensation remain UNKNOWN for all three.

## Scope, reproducibility and limitations

`ema_identity_screen.py` replays `ema-recovered-aliases.json` against explicit catalog/corpus paths, writing `ema-identity-match-rows.tsv`/JSON, a header-only `ema-identity-mentions.tsv`, and `ema-identity-screen-summary.json`. Every catalog input and the complete corpus is hashed. Corpus SHA256 is `e5b8b017bb68af38818f1cc8c1bee47adb16afe0367981104a6b17886d2954d9`. The screen inspects saved titles/abstracts, not every full text.

Direct disease/gene searches returned zero in valid responses, with original query translations retained. PubMed drops or expands some unindexed development names; Waskyra's 11 broad muscle results were not literal product matches. Kostaive and garadacimab broader muscle/target hit sets were used for discovery and were not exhaustively adjudicated. Garadacimab's first Europe PMC disease response lacked `hitCount`; that failed response is explicitly separated from the successful repeat.

The current work does not amend the historical **180 unmatched EMA product rows**, change the **15,357 census denominator**, or certify all European/national/global coverage. It resolves documented names and establishes a transparent follow-up screen. Raw records remain outside Git in the three medicine folders; `ema-identity-evidence.json` and each source manifest provide exact official/primary links and evidence depth.


Replay the saved inputs without a network request:

```sh
python3 workflow/track2-drug-census/ema_identity_screen.py \
  --catalog-dir workflow/track2-drug-census \
  --corpus /path/to/literature-gene-unrestricted.json \
  --aliases workflow/track2-drug-census/ema-recovered-aliases.json \
  --output-dir /tmp/mva-ema-identity-replay
```

The per-alias table and mention table are deterministic. The summary records a
new generation timestamp on replay; compare all other fields and input hashes.
The historical EMA product table retains its original unassessed-efficacy status;
this separate follow-up supplies the three medicine-specific assessments.
