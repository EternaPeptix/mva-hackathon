# Food-substance discovery layer

**Retrieved 23 September 2026; source page updated 28 July 2026.** This layer
preserves all **3,971 rows** in the retrieved [FDA Substances Added to Food
table](https://www.hfpappexternal.fda.gov/scripts/fdcc/index.cfm?order=ASC&search=&set=FoodSubstances&showAll=true&sort=Sortterm_ID&type=basic).
It broadens discovery beyond marketed medicines and specialist research probes.
It is not a list of safe treatments or every food compound in existence.

## Coverage and source status

| Measure | Verified snapshot |
|---|---:|
| Preserved source rows / source record IDs | 3,971 / 3,971 |
| Exact clinical-catalog name candidates | 298 rows |
| Exact curated-appraisal name candidates | 7 rows |
| Literal mention rows in the 7,362-record saved corpus | 234 |
| Food source rows with at least one literal hit | 76 |
| Matched literature source/ID records | 213 |
| Common/short/class-name ambiguous mention rows | 161 |
| Literal-policy rows, now covered by separate abstract overlay | 73 |
| New manual-reading / curated-appraisal credit from this inventory | 0 / 0 |

The table retains **17 PROHIBITED**, **3 PROHIBITED WITH EXCEPTIONS**,
**8 DELISTED** and **8 NLFG** suffixes. NLFG means no longer FEMA GRAS in the
source's terminology. The other 3,935 names lack those explicit suffixes;
that is not a finding that they are currently permitted, approved or safe.
The [FDA scope statement](https://www.fda.gov/food/food-additives-petitions/substances-added-food-formerly-eafus)
explains that this is a partial inventory and that non-FDA evaluations are not
FDA approvals. Many GRAS notices and most indirect food-contact substances
are in other inventories.

Source names, technical-effect text, CAS-or-other-ID fields, regulation links
and detail URLs are preserved. Mixtures, extracts, enzymes and generic classes
remain source records, not chemically resolved single molecules. Two ID cells
are blank; other identifiers must not all be relabeled as CAS registry numbers.
No current stock or individual food-use authorization was checked.

## What a match means

Matching uses full literal source names after case/Unicode/whitespace
normalization. Only an explicit regulatory suffix is removed for retrieval;
the original label and status stay intact. There is no synonym, fuzzy,
structure, CAS or ingredient-component equivalence. The 3,971 rows therefore
overlap existing catalogs and do not increase the 15,357-entity denominator.

These examples show why the [retrieval rows](mentions.tsv) are questions, not
automatically identified interventions:

- **CASEIN** can occur inside “casein kinase,” a protein name.
- **L-CYSTEINE** can occur inside N-acetyl-L-cysteine or S-trityl-L-cysteine,
  which are different compounds.
- **HEXYL ISOTHIOCYANATE** occurs inside the name of a modified wasabi compound.
- **Collagen**, fatty acids and amino acids can be endogenous molecules,
  matrix/media components or experimental readouts.
- A name cross-reference to a mixture appraisal does not assign the mixture's
  effect to the isolated component.

The common-name flags help triage but do not detect every ambiguity. A hit can
be a harmful exposure, control, computational prediction or cited background.
A miss can reflect an unlisted synonym or a name only present in full text.
All entries retain unknown BUB1B efficacy at this inventory layer. Prior reading
in another layer is neither erased nor counted again.

## Replay and validation

From the repo root, using the externally saved source bytes:

```sh
python3 workflow/track2-drug-census/food_substance_inventory.py \
  --html /path/food-substances-all.html \
  --corpus /path/literature-gene-unrestricted.json \
  --out /path/replayed-food-inventory
```

[summary.json](summary.json) records input, code and output hashes. The parser
checks the source-declared row count and table header, rejecting a partial page.
It preserves duplicate labels rather than merging them and stages outputs
before rechecking inputs. An input-change failure leaves the previous package
untouched. Publishing the three files is not a filesystem transaction: the
summary is written last, and its hashes detect an interrupted replacement.

Five regression tests cover status suffixes, chemical-name separators,
partial/changed tables, duplicate/blank identities, ambiguity flags and failed
rebuild preservation. A separate agent independently reconciled all source
cells and literal match keys; its scope and fix recheck are recorded in
[review.json](review.json). These are data-handling checks, not independent
clinical review or efficacy validation.

The [lunasin qualification](../lunasin-qualification.md) is a separate,
source-driven food-peptide result recovered through literature review. It is
not evidence that this FDA table contains every relevant food-derived candidate.

## Abstract reading and two primary qualifications

The [separate reading overlay](abstract-findings.md) covers all 67 source/ID records behind the 73 literal-policy rows. It retains derivatives, solvents, endogenous molecules, retraction and version distinctions. The 161 common-name rows remain outside this reading pass. Eight records overlap prior accepted revalidations and two overlap batch08; 26 others occur in the remaining queue. No automatic outer-corpus credit is added.

[Astaxanthin](../astaxanthin-qualification.md) and [citral](../citral-qualification.md) received separate primary qualifications. These add two appraisals and two name cross-references, not demonstrated MVA treatments. The original independent parser review remains bound to its 161-entry snapshot; the parent recheck verifies only the current two-link update, with source fields and all 234 mention rows unchanged.
