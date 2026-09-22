# EMA direct inventory: bounded source-coverage audit

The current census has a useful regulatory coverage gap: its DrugCentral EMA list has 456 rows from an unversioned export, while the official current EMA feed supplies direct product status and initial-authorisation dates. This layer should supplement the source-bounded census. It does not establish coverage of every drug, every national European approval, or unique active substances.

Source: European Medicines Agency (EMA), accessed September 2026. The official [JSON guide](https://www.ema.europa.eu/en/about-us/about-website/download-website-data-json-data-format) identifies the [medicines feed](https://www.ema.europa.eu/en/documents/report/medicines-output-medicines_json-report_en.json) as centralised-procedure medicine-page data. It is published twice daily; the alternative XLSX carries the same data, so downloading both would add no source coverage. This audit evaluated this one inventory and selected individual linked EPAR pages; no gene literature search was performed.

## Measured result

The preserved feed has **2,746 product records**, all with distinct nonblank EMA product IDs: **2,351 human** and **395 veterinary**. The derived table contains exactly **1,568 currently `Authorised` human product rows**. It excludes 783 human rows with other statuses and all veterinary rows. There are 32 `Opinion` and two `Opinion under re-examination` human rows in the source. Withdrawn, refused, suspended, revoked, expired and lapsed statuses stay outside the authorised table.

Among the 1,568 authorised human products, 1,388 have at least one exact whole-field name cross-reference to the existing census/supplements; **180 do not**. This is a triage result, not 180 missing drugs: combinations, descriptive vaccine/cell fields, new brands and missing synonyms contribute. The existing **15,357 entity denominator remains unchanged**.

The comparison uses NFKC/casefold/whitespace normalization and separately preserves matches for product name, active substance and INN/common name against (a) main preferred names and aliases, (b) status-unverified DrugCentral names/aliases, (c) FDA2026 supplement ingredient/brand/search aliases, and (d) Broad preferred names. No combination splitting, salt removal, fuzzy matching or identity merging occurs in the complete product table. All equivalence remains `UNKNOWN_NOT_EVALUATED`.

## Small manually checked comparison

| Product / EMA product ID | Official active substance or INN | Initial EU authorisation | Coverage finding |
|---|---|---|---|
| [Andembry](https://www.ema.europa.eu/en/medicines/human/EPAR/andembry), EMEA/H/C/006116 | garadacimab | 2025-02-10 | No exact product/substance/INN name match in any of the four compared layers. Plausible missing-name/substance candidate; identity remains UNKNOWN. |
| [Kostaive](https://www.ema.europa.eu/en/medicines/human/EPAR/kostaive), EMEA/H/C/006207 | zapomeran; feed active-substance field says COVID-19 sa-mRNA Vaccine | 2025-02-12 | Product page explicitly identifies zapomeran. No exact match in any layer. Vaccine identity must remain distinct from a generic COVID-19 vaccine class. |
| [Waskyra](https://www.ema.europa.eu/en/medicines/human/EPAR/waskyra), EMEA/H/C/006525 | etuvetidigene autotemcel; feed also retains the long modified-cell description | 2026-01-09 | No exact match in any layer. Gene/cell therapy; no active-moiety equivalence assumed. |
| [Kygevvi](https://www.ema.europa.eu/en/medicines/human/EPAR/kygevvi), EMEA/H/C/005119 | doxecitine; doxribtimine | 2026-03-26 | Whole product field unmatched, but separately inspected components already occur as CHEMBL66115 and CHEMBL52609. This is an authorisation/status update opportunity, not two new census drugs. |
| [Zemcelpro](https://www.ema.europa.eu/en/medicines/human/EPAR/zemcelpro), EMEA/H/C/005772 | non-expanded allogeneic umbilical cord-derived CD34- cells; dorocubicel | 2025-08-25 | Dorocubicel already occurs as CHEMBL5314825. Full product/component equivalence is unresolved. |
| [Frehemgo](https://www.ema.europa.eu/en/medicines/human/EPAR/frehemgo), EMEA/H/C/006344 | denecimig | None in source | Positive opinion dated 2026-09-17, not an authorisation; denecimig already matches CHEMBL5314670. Excluded from authorised layer. |
| [Alofisel](https://www.ema.europa.eu/en/medicines/human/EPAR/alofisel), EMEA/H/C/004258 | darvadstrocel | 2018-03-23 | Source status is Withdrawn, withdrawal date 2024-12-13. Already matches CHEMBL4298180. Excluded from authorised layer. |

The first three official product pages independently agree with the feed's substance/INN and initial-authorisation dates. Their downloaded HTML is retained with hashes. These are source coverage findings, not evidence of BUB1B benefit or a treatment recommendation.

## Date and source traps

Use `marketing_authorisation_date` for initial authorisation. `european_commission_decision_date` can concern a later change: Andembry has initial 10 February 2025 but EC field 9 July 2026; Kostaive has initial 12 February 2025 but EC field 27 November 2025. Conversely, Kygevvi's EC field is blank despite an initial-authorisation date confirmed on its EPAR page. Six authorised human rows have blank initial-authorisation dates; they stay blank.

The feed was retrieved at **2026-09-22T05:54:31.996026+00:00** (21 September in Phoenix), with HTTP `Last-Modified: Tue, 22 Sep 2026 04:02:48 GMT`. Its raw `meta.timestamp` says **2026-09-22T06:02:48Z**, which is later than retrieval and differs from Last-Modified by two hours. This clock discrepancy is preserved and flagged, not corrected or interpreted as evidence of a future update. No explicit source schema/version is provided; observed root/meta keys and record fields are recorded.

## Source reuse and reproduction

Source: European Medicines Agency (EMA), accessed September 2026. EMA's
[legal notice](https://www.ema.europa.eu/en/about-us/about-website/legal-notice)
permits reuse of EMA-owned public information with attribution in each copy;
third-party material has separate terms. This is not an asserted Creative
Commons licence or EMA endorsement. Attribution accompanies each derived row.

The [product table](ema-authorised-human-products.tsv),
[coverage summary](ema-coverage-summary.json), [checked examples](ema-candidate-review.json)
and [retrieval metadata](ema-retrieval-metadata.json) preserve source fields and
original observations. Raw feed/page responses stay outside Git. The public
[comparison script](ema_coverage.py) makes no network calls and rebuilds the table
from a supplied original feed and the existing catalog inputs:

```sh
python3 workflow/track2-drug-census/ema_coverage.py \
  --ema-json /path/to/saved-ema-feed.json \
  --catalog-dir workflow/track2-drug-census \
  --retrieval-metadata workflow/track2-drug-census/ema-retrieval-metadata.json \
  --output-dir /tmp/mva-ema-replay
```

A live download can change; the archived retrieval metadata only identifies the
saved feed with SHA256 `a2acab9410d80517e58a131d45103e8184a36dd23752cea0453e674bd2ab0cfa`.
Fresh inputs need their own retrieval record. The script checks the feed hash,
source IDs and fields, records catalog-input hashes, and preserves product-level matches.
An independent reconstruction confirmed the 1,568/1,388/180 counts and the
combination/date traps. These checks establish coverage, not substance identity
or efficacy. No new literature search or BUB1B efficacy decision was assigned
to the EMA-only rows.

Next, resolve the three named candidate identities against authoritative
substance identifiers and triage the other unmatched product fields. The
comparison does not prove that any entity is absent under every possible alias,
and 180 whole-field mismatches are not 180 new medicines. Publication and
submission remain pending while scientific and source-coverage work continues.
