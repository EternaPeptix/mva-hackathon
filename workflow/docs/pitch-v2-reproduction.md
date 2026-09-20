# Reproduce the second narrated pitch template

This is an MVA-specific schematic template and a local media-integrity workflow.
It produces a three-minute, 1920×1080, 24 fps research draft. A render is not a
scientific result, an independently reviewed narration, or an uploaded submission.
The source snapshot omits the case-specific storyboard, clinical reports, movie and
case production receipts. Supply an authorized, reviewed storyboard separately.

## Environment and inputs

The completed production used an Apple Silicon Mac Studio, macOS 26.5.1 (25F80),
Swift 6.3.3, Python 3.12.14, Pillow 12.3.0 and NumPy 2.3.5. It used the installed
Samantha voice, `/usr/bin/say`, `/usr/bin/afconvert`, and the system Arial and Arial
Bold fonts under `/System/Library/Fonts/Supplemental/`. Swift uses Apple's native
AVFoundation, AppKit, CoreGraphics, CoreVideo and ImageIO frameworks. It does not
require a GPU service, Spark, FFmpeg, a network call or an additional voice download.

Use a separate Python environment with the stated libraries available. Verify the
voice, fonts and compiler before production. Changing the OS, voice resources,
compiler or libraries can change narration timing and media bytes; fresh checks are
required. Recorded hashes identify one production, not universal reproducible bytes.

The storyboard is UTF-8 JSON with:

- `title`, `persistent_footer`, and `duration_seconds: 180`;
- `synthetic_narration.voice` and `synthetic_narration.initial_words_per_minute`;
- six `scenes`, each with `id`, `title`, `kicker`, `display_source`, and exactly four
  strings in `sentences`.

The executed scene order is `purpose`, `genotype`, `mechanism`, `candidates`,
`measurement`, `oncology`. The preparer also retains the older `decision` layout.
The diagrams and several labels are fixed in the Python source, so changes to a
storyboard require checking its agreement with those diagrams. This is not an
arbitrary slide-layout engine. Optional `sources`, per-scene `source_keys`, status
and release notes can preserve a claim audit, but the renderer does not validate
their scientific content. The visible draft and synthetic-narration labels must
remain accurate. Gene/protein names are expanded only in synthesized speech;
captions preserve the supplied text.

## Produce and check a fresh output

Set absolute paths to the checkout, interpreter, authorized storyboard and a new
output directory outside the checkout. These commands are sequential: stop on a
nonzero exit status and retain failed outputs before choosing another new path.

```sh
MVA_REPO=/absolute/path/to/source-checkout
MVA_PY=/absolute/path/to/python3
MVA_STORY=/absolute/path/to/authorized-storyboard.json
MVA_PITCH=/absolute/path/to/new-pitch-production

mkdir "$MVA_PITCH"
"$MVA_PY" "$MVA_REPO/workflow/scripts/prepare_track2_pitch_v2.py" \
  --storyboard "$MVA_STORY" --output "$MVA_PITCH/prepared"

/usr/bin/swiftc -O "$MVA_REPO/workflow/scripts/render_track2_pitch_v2.swift" \
  -o "$MVA_PITCH/render_track2_pitch_v2"
"$MVA_PITCH/render_track2_pitch_v2" \
  "$MVA_PITCH/prepared/video-plan.json" "$MVA_PITCH/rendered"

/usr/bin/swiftc -O "$MVA_REPO/workflow/scripts/inspect_track2_pitch.swift" \
  -o "$MVA_PITCH/inspect_track2_pitch"
"$MVA_PITCH/inspect_track2_pitch" "$MVA_PITCH/prepared/video-plan.json" \
  "$MVA_PITCH/rendered/jeweled_track2_pitch_draft_20260920.mp4" \
  "$MVA_PITCH/inspection"

"$MVA_PY" "$MVA_REPO/workflow/scripts/verify_track2_pitch.py" \
  --storyboard "$MVA_STORY" --prepared "$MVA_PITCH/prepared" \
  --rendered "$MVA_PITCH/rendered" --inspection "$MVA_PITCH/inspection" \
  --prepare-source "$MVA_REPO/workflow/scripts/prepare_track2_pitch_v2.py" \
  --output "$MVA_PITCH/verified"
```

Each production stage requires a previously absent output directory. The render's
fixed filename is a template identifier, not a claim about a later execution date.
Plans use absolute artifact paths; preserve that directory tree or regenerate it
in the new location. Do not rewrite a completed plan without recording the change.

The preparer records voice commands, per-sentence audio and measured timing,
WebVTT captions, image bounds, source/font hashes and the contiguous frame plan.
The native renderer encodes H.264 video, combines audio, checks the actual media
sample count and duration, and decodes the audio. Empty compressed control buffers
are recorded separately from video samples.

The independent inspection step decodes one midpoint for every planned static
segment. The final verifier checks prepared-file identities, caption text/order and
timing, 180 seconds of decoded audio, all 24 sentence windows, and every segment
midpoint against its source image. Audio requires correlation at least 0.98 and
per-window RMS ratios of 0.8–1.2. Image comparison resizes both sides to 480×270 RGB
and requires mean absolute error at most 4, RMSE at most 10, and the 99th percentile
absolute error at most 32 on the 0–255 scale. These are media-integrity thresholds,
not readability or scientific-quality scores.

Retain `preparation.json`, `video-plan.json`, `render.json`, `inspection.json` and
`verification.json`, source and executable hashes, runtime/voice/font identities,
and the original movie, captions and audio. Inspect the decoded contact sheet and
full-size frames for layout. Listen to the complete movie for pronunciation and
intelligibility and review scientific claims separately. Midpoint checks do not
exhaustively inspect every frame, and waveform preservation does not establish
that the spoken words are correct or understandable. No step uploads the movie.
