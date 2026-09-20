#!/usr/bin/env python3
"""Check a local pitch's captions, decoded audio and all static segment midpoints.

This is media-integrity QA, not validation of scientific claims or pronunciation.
The native renderer supplies decoded audio; inspect_track2_pitch.swift supplies
one decoded midpoint per planned segment. Existing output directories are refused.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import wave

import numpy as np
from PIL import Image


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def require(value, message):
    if not value:
        raise ValueError(message)


def read_wave(path):
    with wave.open(str(path), "rb") as reader:
        require((reader.getnchannels(), reader.getsampwidth(), reader.getframerate()) == (1, 2, 44100), "Unexpected audio format")
        return np.frombuffer(reader.readframes(reader.getnframes()), dtype="<i2").astype(np.float64)


def vtt_seconds(value):
    h, m, s = value.split(":")
    return int(h)*3600 + int(m)*60 + float(s)


def run(args):
    require(not args.output.exists(), "Existing output preserved; choose a new QA directory")
    story = read_json(args.storyboard)
    prep = read_json(args.prepared / "preparation.json")
    plan = read_json(args.prepared / "video-plan.json")
    rendered = read_json(args.rendered / "render.json")
    inspection = read_json(args.inspection / "inspection.json")
    require(prep["status"] == rendered["status"] == inspection["status"] == "complete", "Incomplete production stage")
    require(prep["storyboard_sha256"] == sha(args.storyboard), "Storyboard changed since preparation")
    require(prep["code_sha256"] == sha(args.prepare_source), "Preparation source identity differs")
    for filename, info in prep["files"].items():
        path = args.prepared / filename
        require(path.stat().st_size == info["bytes"] and sha(path) == info["sha256"], "Prepared artifact changed: " + filename)
    require(story["duration_seconds"] == plan["duration_seconds"] == 180, "Wrong requested duration")
    require(plan["draft"] and plan["fps"] == 24 and plan["frames"] == 4320, "Wrong draft plan")
    require(rendered["width"] == 1920 and rendered["height"] == 1080 and rendered["video_frames"] == 4320, "Wrong encoded video dimensions or sample count")
    require(abs(rendered["duration_seconds"]-180) < .001 and abs(rendered["audio_track_duration_seconds"]-180) < .001, "Wrong media duration")
    require(abs(rendered["nominal_frame_rate"]-24) < .001, "Wrong frame rate")
    movie = Path(rendered["video"])
    require(movie.is_file() and Path(inspection["movie"]).resolve() == movie.resolve(), "Inspected movie identity differs")

    captions = plan["captions"]
    expected = [(i, j, text) for i, scene in enumerate(story["scenes"]) for j, text in enumerate(scene["sentences"])]
    require([(c["scene"], c["sentence"], c["caption"]) for c in captions] == expected, "Caption text/order differs from storyboard")
    require(len(captions) == prep["captions"] == 24, "Missing caption")
    require(prep["narrated_words"] == sum(len(t.split()) for _, _, t in expected), "Word count differs")
    blocks = (args.prepared / "captions.vtt").read_text().strip().split("\n\n")
    require(blocks[0] == "WEBVTT" and len(blocks)-1 == len(captions), "Wrong WebVTT structure")
    for number, (block, caption) in enumerate(zip(blocks[1:], captions), 1):
        lines = block.splitlines()
        start, end = lines[1].split(" --> ")
        require(lines[0] == str(number) and "\n".join(lines[2:]) == caption["caption"], "WebVTT text/order differs")
        require(abs(vtt_seconds(start)-caption["start_seconds"]) <= .00051 and abs(vtt_seconds(end)-caption["end_seconds"]) <= .00051, "WebVTT timing differs")

    original = read_wave(Path(plan["audio"]))
    decoded = np.fromfile(rendered["decoded_audio_PCM"], dtype="<i2").astype(np.float64)
    require(len(original) == len(decoded) == rendered["decoded_audio_samples"] == 180*44100, "Decoded audio sample count differs")
    require(np.max(np.abs(decoded)) < 32767, "Integer clipping detected")
    global_correlation = float(np.corrcoef(original, decoded)[0, 1])
    require(np.isfinite(global_correlation) and global_correlation >= .98, "Overall audio waveform changed")
    audio_rows = []
    previous_end = 0
    for caption in captions:
        start, end = caption["start_sample"], caption["end_sample"]
        require(0 <= previous_end <= start < end <= len(decoded), "Audio windows overlap or escape media")
        previous_end = end
        clip = Path(caption["wav"])
        require(sha(clip) == caption["sha256"], "Speech source changed")
        require(len(read_wave(clip)) == end-start, "Speech source length differs")
        a, b = original[start:end], decoded[start:end]
        correlation = float(np.corrcoef(a, b)[0, 1])
        original_rms, decoded_rms = float(np.sqrt(np.mean(a*a))), float(np.sqrt(np.mean(b*b)))
        ratio = decoded_rms / original_rms if original_rms else float("nan")
        require(np.isfinite(correlation) and correlation >= .98 and .8 <= ratio <= 1.2, "Missing or altered caption audio window")
        audio_rows.append({"scene":caption["scene"], "sentence":caption["sentence"], "start_seconds":start/44100, "end_seconds":end/44100, "waveform_correlation":correlation, "decoded_to_source_RMS_ratio":ratio})

    rows = inspection["segments"]
    require(len(rows) == len(plan["segments"]), "A planned segment was not sampled")
    image_rows = []
    for index, (observed, intended) in enumerate(zip(rows, plan["segments"])):
        require(observed["segment_index"] == index, "Segment order differs")
        for key in ("start_frame", "end_frame"):
            require(observed[key] == intended[key], "Segment interval differs")
        midpoint = (intended["start_frame"] + intended["end_frame"] - 1)//2
        require(observed["sampled_frame"] == midpoint and abs(observed["actual_seconds"]-midpoint/24) < 1/48, "Sampled midpoint differs")
        require(Path(observed["source_image"]).resolve() == Path(intended["image"]).resolve(), "Wrong intended image")
        with Image.open(intended["image"]) as src, Image.open(observed["decoded_image"]) as dst:
            require(src.size == dst.size == (1920,1080), "Image dimensions differ")
            a = np.asarray(src.convert("RGB").resize((480,270), Image.Resampling.LANCZOS)).astype(np.float64)
            b = np.asarray(dst.convert("RGB").resize((480,270), Image.Resampling.LANCZOS)).astype(np.float64)
        delta = np.abs(a-b)
        mae, rmse, p99 = float(np.mean(delta)), float(np.sqrt(np.mean(delta*delta))), float(np.percentile(delta,99))
        require(mae <= 4 and rmse <= 10 and p99 <= 32, "Encoded image differs materially from planned segment")
        image_rows.append({"segment_index":index, "sampled_frame":midpoint, "source_sha256":sha(Path(intended["image"])), "decoded_sha256":sha(Path(observed["decoded_image"])), "RGB_mean_absolute_error":mae, "RGB_RMSE":rmse, "RGB_99th_percentile_absolute_error":p99})

    args.output.mkdir()
    # Actual decoded caption frames, for a visual pass in addition to numeric QA.
    sheet = Image.new("RGB", (1920,2160), "white")
    for i, caption in enumerate(captions):
        match = [row for row in rows if Path(row["source_image"]).resolve() == Path(caption["caption_image"]).resolve()]
        require(len(match) == 1, "Caption frame not uniquely inspected")
        with Image.open(match[0]["decoded_image"]) as frame:
            sheet.paste(frame.convert("RGB").resize((480,270)), ((i%4)*480, (i//4)*360))
        # The 90-pixel margin retains room between the six rows of scenes.
    sheet.save(args.output / "decoded-caption-contact-sheet.png")
    receipt = {"status":"complete", "verified_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
               "movie":str(movie), "movie_bytes":movie.stat().st_size, "movie_sha256":sha(movie),
               "storyboard_sha256":sha(args.storyboard), "preparation_source_sha256":sha(args.prepare_source), "verifier_source_sha256":sha(Path(__file__)),
               "duration_seconds":180, "video_frame_samples":4320, "audio_samples":len(decoded),
               "caption_count":len(captions), "word_count":prep["narrated_words"], "decoded_audio_waveform_correlation":global_correlation,
               "decoded_audio_peak_dBFS":float(20*np.log10(np.max(np.abs(decoded))/32767)),
               "audio_windows":audio_rows, "segment_image_checks":image_rows,
               "limits":["Decoded midpoint per static segment, not exhaustive frame-by-frame inspection", "Waveform preservation does not independently validate pronunciation or speech intelligibility", "Scientific claims require separate review", "Local research draft, not uploaded and not clinical evidence"]}
    (args.output / "verification.json").write_text(json.dumps(receipt,indent=2,allow_nan=False)+"\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("storyboard", "prepared", "rendered", "inspection", "prepare-source", "output"):
        parser.add_argument("--"+key, required=True, type=Path)
    result = run(parser.parse_args())
    print(json.dumps({key:result[key] for key in ("status","duration_seconds","caption_count","word_count","decoded_audio_waveform_correlation","movie_sha256")}))
