# Cell selection and tissue output

Open [index.html](index.html) directly in a browser after cloning or downloading
this repository. Keep `index.html`, `style.css`, `app.js` and `data.js` together.
The tool works offline without installation, a web server, uploads or analytics.
It has not been deployed to a public host.

The default retained example shows a **16.8-percentage-point rise in the
chromosome-normal fraction alongside an 11.0% fall in raw differentiated-cell
time**. Absolute accumulated output is **51,110.95 → 45,513.34 cell × model-time
units**. This is an assumed population-selection result, not measured tissue
function, a drug effect or a prediction for a patient.

## What the controls do

The explorer looks up a slice of [retained results](../population-results.json).
It does not run a new simulation or interpolate between settings. Four controls
select renewal, abnormal-to-normal progenitor division ratio (0.7 or 1.8), error
probability (0.005, 0.02 or 0.04), and horizon (20 or 80 arbitrary model-time units).
An additional selector chooses one of the eight existing hypothetical arms.
Every arm is compared with its matched untreated reference.

The slice has 24 scenarios × eight arms, with the following original assumptions
fixed: initial senescence 0.25, senescence impairment 1, senescent contribution
0.25, senescent niche occupancy 0, abnormal-cell contribution 0.5, off-target
fraction 0.02 and crowding capacity 2,000. The model starts with 4,000 cells,
20% abnormal. These are assumptions; none is estimated from the case genome.

The two main endpoints are chromosome-normal fraction among all surviving cells
and raw differentiated-cell time. Absolute surviving counts and optional weighted
functional cell-time are also shown. All three population compartments retain the
constitutional genotype. These are aggregate compartments, **not individual
lineages or pedigrees**. A higher chromosome-normal share cannot establish genetic
correction or functional rescue. The tool makes no claim of unique methodological
novelty. See the [model documentation](../README.md) for the scope and limitations.

## Reproduce the packaged data

From the repository root, with Python 3.9+ and no third-party packages:

```sh
python3 workflow/track2-model/explorer/build_data.py
```

This regenerates only `explorer/data.js` from the existing retained JSON. It checks
the model source hash, verifies all 192 selected rows' fractions, counts and ratios,
and compares the default against the separately retained illustrative result. It
fails if the published +16.8-point / −11.0% example changes. It preserves source
numeric precision and embeds source/result SHA-256 values. Repeated generation
from the same source is byte-identical; it does not run `population.py` or modify
its source results.

**Save scenario JSON** exports both full result rows, the selected/fixed assumptions,
the intervention definition and source hashes. It contains no uploaded data.

## Optional browser verification

With Node.js, Playwright and its Chromium runtime already installed:

```sh
node workflow/track2-model/explorer/check_browser.cjs
```

This checks all 192 selectable comparisons against the stored rows, the default
headline and source example, exported values, keyboard operation, absence of
external requests, and mobile horizontal overflow. Desktop and mobile screenshots
are written to a temporary directory printed in the result for visual review.
These checks validate the explorer, not the biological model. Playwright is used
only for verification and is not needed by the offline tool.
Set `EXPLORER_CHROMIUM` to use an existing Chromium executable if it is not in
Playwright's default browser location.
