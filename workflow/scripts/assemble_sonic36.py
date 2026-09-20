#!/usr/bin/env python3
"""Assemble pinned, complete scene speech without synthesis, trimming or resampling."""
import argparse
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
import unicodedata
import wave

sys.dont_write_bytecode = True
SR, FPS, DURATION = 44100, 24, 180
EPS = 0.00001


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pin(path):
    path = Path(path)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha(path)}


def bind(info):
    path = Path(info["path"])
    require(path.is_file() and path.stat().st_size == info["bytes"] and sha(path) == info["sha256"], "Pinned input changed: " + path.name)
    return path


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def spoken(text):
    return text.replace("BUB1B", "Bub one B").replace("BUBR1", "Bub R one").replace("HSP70", "H S P seventy")


def tokens(text):
    # Keep punctuation and spelling; normalization is Unicode and whitespace only.
    return unicodedata.normalize("NFKC", text).split()


def align_sentences(canonical, said, timestamps, duration):
    require(len(canonical) == len(said) == 4, "Exactly four canonical sentences required")
    require(said == [spoken(text) for text in canonical], "Spoken text differs from declared name expansion")
    words, starts, ends = (timestamps[k] for k in ("words", "start", "end"))
    require(len(words) == len(starts) == len(ends) > 0, "Word/timestamp length mismatch")
    flat, provenance = [], []
    for index, (word, start, end) in enumerate(zip(words, starts, ends)):
        require(isinstance(word, str) and bool(tokens(word)), "Empty/nontext timestamp word")
        require(type(start) in (int, float) and type(end) in (int, float) and math.isfinite(start) and math.isfinite(end), "Invalid timestamp scalar")
        require(0 <= start <= end + EPS and end <= duration + EPS, "Timestamp outside raw scene audio")
        if index:
            require(start + EPS >= starts[index-1] and end + EPS >= ends[index-1], "Nonmonotone provider timestamps")
        for token in tokens(word):
            flat.append(token)
            provenance.append(index)
    expected = [token for sentence in said for token in tokens(sentence)]
    require(flat == expected, "Provider word stream does not exactly match requested spoken text")
    result, cursor = [], 0
    for number, (caption, sentence) in enumerate(zip(canonical, said)):
        count = len(tokens(sentence))
        first, last = provenance[cursor], provenance[cursor+count-1]
        require(cursor == 0 or provenance[cursor-1] != first, "A provider segment crosses a sentence boundary")
        require(cursor+count == len(flat) or provenance[cursor+count] != last, "A provider segment crosses a sentence boundary")
        start, end = starts[first], ends[last]
        require(start < end, "Sentence has no positive timestamp extent")
        if result:
            require(start + EPS >= result[-1]["source_end_seconds"], "Sentence timestamp ranges overlap")
        result.append({"sentence": number, "caption": caption, "spoken_text": sentence,
                       "word_start_index": first, "word_end_index_exclusive": last+1,
                       "token_start_index": cursor, "token_end_index_exclusive": cursor+count,
                       "source_start_seconds": start, "source_end_seconds": end})
        cursor += count
    require(cursor == len(flat), "Unassigned spoken tokens")
    return result


def frame_boundaries(samples, maximum_padding=45):
    require(len(samples) == 6 and all(type(x) is int and x > 0 for x in samples), "Six nonempty full-scene audio streams required")
    added = DURATION - sum(samples)/SR
    require(6*0.7 <= added <= maximum_padding, "Speech duration leaves insufficient or excessive reading padding")
    boundaries = [0]
    running = 0.0
    for count in samples:
        running += count/SR + added/6
        boundaries.append(round(running*FPS))
    boundaries[-1] = DURATION*FPS
    for index, count in enumerate(samples):
        begin = round((boundaries[index]/FPS + .7)*SR)
        require(begin + count <= round(boundaries[index+1]/FPS*SR), "Complete scene would exceed assigned reading window")
    return boundaries, added


def load_inputs(path):
    import numpy as np
    config = json.loads(path.read_text())
    require(config["status"] == "ready" and config["schema"] == "sonic36-pitch-assembly-inputs-v1", "Unready or wrong input schema")
    require(config["duration_seconds"] == DURATION and config["sample_rate"] == SR and config["fps"] == FPS, "Unexpected media geometry")
    all_pins = [config["storyboard"]] + config["request_manifests"]
    for row in config["vendor"]:
        require(row["original"]["sha256"] == row["copy"]["sha256"], "Vendor copy changed")
        all_pins.extend([row["original"], row["copy"]])
    for source in config["sources"]:
        all_pins.extend(source["files"].values())
    for info in all_pins:
        bind(info)
    story = json.loads(bind(config["storyboard"]).read_text())
    require(story["duration_seconds"] == DURATION and len(story["scenes"]) == 6, "Wrong canonical story")
    requests = []
    for info in config["request_manifests"]:
        manifest = json.loads(bind(info).read_text())
        require(manifest["storyboard_sha256"] == config["storyboard"]["sha256"], "Request storyboard binding differs")
        requests.extend(manifest["requests"])
    require(len(requests) == len({r["name"] for r in requests}) == 6, "Duplicate/missing scene requests")
    reqs = {r["name"]:r for r in requests}
    rows = []
    require([r["id"] for r in config["sources"]] == [s["id"] for s in story["scenes"]], "Scene order changed")
    for scene, source in zip(story["scenes"], config["sources"]):
        files = source["files"]
        request = json.loads(bind(files["request.json"]).read_text())
        receipt = json.loads(bind(files["receipt.json"]).read_text())
        stamp = json.loads(bind(files["word-timestamps.json"]).read_text())
        spec = reqs[source["id"]]
        require(spec == source["request"] and spec["canonical_sentences"] == scene["sentences"], "Scene canonical/source manifest mismatch")
        require(request == spec["payload"] and request["transcript"] == " ".join(spec["spoken_sentences"]), "Requested transcript differs")
        require(request["model_id"] == "sonic-3.6" and request["voice"] == config["synthetic_voice"]["voice_id"], "Voice/model changed")
        require(request["output_format"] == {"container":"raw", "encoding":"pcm_s16le", "sample_rate":SR}, "Wrong requested audio format")
        require(request["generation_config"]["speed"] == config["synthetic_voice"]["speed"], "Provider speed differs")
        require(receipt["status"] == "complete" and receipt["requested_model"] == "sonic-3.6", "Incomplete synthesis receipt")
        require(receipt["request_sha256"] == files["request.json"]["sha256"] and receipt["voice_id"] == request["voice"], "Synthesis request binding differs")
        for name, digest in receipt["files"].items():
            require(name in files and files[name]["sha256"] == digest, "Synthesis file binding differs")
        audio = bind(files["narration.wav"])
        with wave.open(str(audio), "rb") as stream:
            require((stream.getframerate(), stream.getnchannels(), stream.getsampwidth()) == (SR,1,2), "Wrong raw scene format")
            data = np.frombuffer(stream.readframes(stream.getnframes()), dtype="<i2").copy()
        require(len(data) == receipt["samples"] > 0 and abs(len(data)/SR-receipt["duration_seconds"]) <= 1/SR, "Raw audio count differs")
        peak = int(np.max(np.abs(data.astype(np.int32))))
        require(0 < peak < 32767, "Silent or integer-clipped raw scene")
        align = align_sentences(scene["sentences"], spec["spoken_sentences"], stamp, len(data)/SR)
        require(receipt["timestamp_words"] == len(stamp["words"]), "Receipt word count differs")
        rows.append({"id":source["id"], "raw":data, "audio":str(audio), "samples":len(data), "peak":peak,
                     "alignment":align, "timestamps":stamp, "source_files":files})
    return config, story, rows, all_pins


def run(input_path, output):
    import numpy as np
    from PIL import Image, ImageDraw
    require(not output.exists(), "Existing output preserved; choose a new directory")
    config, story, rows, input_pins = load_inputs(input_path)
    boundaries, padding = frame_boundaries([r["samples"] for r in rows], config["maximum_added_silence_seconds"])
    visual_info = next(r["copy"] for r in config["vendor"] if r["copy"]["path"].endswith("prepare_track2_pitch_v3.py"))
    spec = importlib.util.spec_from_file_location("pinned_old_visuals", bind(visual_info))
    visuals = importlib.util.module_from_spec(spec); spec.loader.exec_module(visuals)
    output.mkdir(parents=True)
    tick = time.monotonic()
    footer = ["AI: Claude Code + Codex • Personal subscriptions • Training disabled (user-confirmed)",
              "RESEARCH DRAFT • No drug-response experiments • Synthetic voice: Cartesia Sonic 3.6 • Cartesia plan/training settings unverified"]
    audio = np.zeros(DURATION*SR, dtype="<i2")
    scenes, captions, segments, all_bounds, base_images = [], [], [], [], []
    for i, (scene, row) in enumerate(zip(story["scenes"], rows)):
        slide = visuals.make_slide(scene, i, "")
        for text, y in zip(footer, [1018,1048]):
            slide.text((96,y), text, visuals.font(18), visuals.MUTED, max_width=1728)
        p = output/f"scene-{i+1:02d}.png"; slide.image.save(p)
        base_images.append(slide.image); all_bounds.extend(slide.bounds)
        start_frame, end_frame = boundaries[i:i+2]
        begin = round((start_frame/FPS+.7)*SR); end = begin+row["samples"]
        require(not np.any(audio[begin:end]), "Scene speech overlaps")
        audio[begin:end] = row["raw"]
        require(np.array_equal(audio[begin:end], row["raw"]), "Raw scene payload changed on placement")
        sr = {"id":row["id"], "image":str(p), "title":scene["title"], "start_frame":start_frame, "end_frame":end_frame,
              "start_seconds":start_frame/FPS, "end_seconds":end_frame/FPS,
              "audio_start_sample":begin, "audio_end_sample":end, "raw_audio":row["audio"],
              "raw_samples":row["samples"], "raw_audio_sha256":sha(row["audio"]), "source_files":row["source_files"],
              "raw_timestamp_payload":row["timestamps"], "reading_padding_before_seconds":begin/SR-start_frame/FPS,
              "reading_padding_after_seconds":end_frame/FPS-end/SR}
        scenes.append(sr)
        last_frame = start_frame
        for aligned in row["alignment"]:
            absolute_start = begin/SR+aligned["source_start_seconds"]
            absolute_end = begin/SR+aligned["source_end_seconds"]
            cstart, cend = math.floor(absolute_start*FPS), math.ceil(absolute_end*FPS)
            require(last_frame <= cstart < cend <= end_frame, "Rounded caption windows overlap or escape scene")
            if last_frame < cstart:
                segments.append({"start_frame":last_frame,"end_frame":cstart,"image":str(p)})
            image = slide.image.copy(); draw = ImageDraw.Draw(image); font = visuals.font(38)
            lines, line = [], ""
            for word in aligned["caption"].split():
                candidate = (line+" "+word).strip()
                if draw.textlength(candidate,font=font)>1636:
                    lines.append(line); line=word
                else:
                    line=candidate
            lines.append(line); require(len(lines)<=3, "Caption too dense")
            text="\n".join(lines); box=draw.multiline_textbbox((140,874),text,font=font,spacing=9)
            require(box[2]<=1784 and box[3]<=995, "Caption does not fit")
            draw.rounded_rectangle((96,856,1824,1005),radius=18,fill=visuals.INK)
            draw.multiline_text((140,874),text,font=font,fill=visuals.WHITE,spacing=9)
            cp=output/f"caption-{len(captions)+1:02d}.png"; image.save(cp)
            segments.append({"start_frame":cstart,"end_frame":cend,"image":str(cp)})
            captions.append({**aligned,"scene":i,"start_seconds":absolute_start,"end_seconds":absolute_end,
                             "start_sample":round(absolute_start*SR),"end_sample":round(absolute_end*SR),
                             "start_frame":cstart,"end_frame":cend,"caption_image":str(cp),
                             "display_lead_seconds":absolute_start-cstart/FPS,"display_tail_seconds":cend/FPS-absolute_end})
            last_frame=cend
        if last_frame<end_frame:
            segments.append({"start_frame":last_frame,"end_frame":end_frame,"image":str(p)})
    require(len(captions)==24 and segments[0]["start_frame"]==0 and segments[-1]["end_frame"]==4320, "Incomplete media plan")
    require(all(a["end_frame"]==b["start_frame"] for a,b in zip(segments,segments[1:])), "Frame gap/overlap")
    wav=output/"narration.wav"
    with wave.open(str(wav),"wb") as writer:
        writer.setnchannels(1);writer.setsampwidth(2);writer.setframerate(SR);writer.writeframes(audio.tobytes())
    vtt="WEBVTT\n\n"+"".join(f"{n+1}\n{visuals.caption_time(c['start_seconds'])} --> {visuals.caption_time(c['end_seconds'])}\n{c['caption']}\n\n" for n,c in enumerate(captions))
    (output/"captions.vtt").write_text(vtt)
    sheet=Image.new("RGB",(1280,1080),"white")
    for i,image in enumerate(base_images):
        sheet.paste(image.resize((640,360)),((i%2)*640,(i//2)*360))
    sheet.save(output/"contact-sheet.png")
    plan={"title":story["title"],"draft":True,"width":1920,"height":1080,"fps":FPS,"duration_seconds":DURATION,
          "frames":4320,"audio":str(wav),"segments":segments,"scenes":scenes,"captions":captions,"persistent_footer":footer}
    save(output/"video-plan.json",plan)
    for info in input_pins:
        bind(info)
    receipt={"status":"complete","completed_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"elapsed_seconds":time.monotonic()-tick,
             "code_sha256":sha(__file__),"storyboard_sha256":config["storyboard"]["sha256"],"input_manifest":pin(input_path),"input_files":input_pins,
             "captions":24,"narrated_words":sum(len(c["caption"].split()) for c in captions),"spoken_tokens":sum(len(tokens(c["spoken_text"])) for c in captions),
             "audio_samples":len(audio),"audio_seconds":len(audio)/SR,"raw_scene_audio_seconds":sum(r["samples"] for r in rows)/SR,
             "added_silence_seconds":padding,"global_audio_gain":1.0,"raw_scene_payloads_preserved":True,"raw_integer_clipping_samples":0,
             "peak_dBFS":float(20*np.log10(np.max(np.abs(audio.astype(float)))/32767)),"synthetic_voice":config["synthetic_voice"],
             "timestamp_numeric_tolerance_seconds":EPS,"VTT_rounding_tolerance_seconds":.0005,"display_rounding_bound_seconds":1/FPS,
             "max_display_lead_seconds":max(c["display_lead_seconds"] for c in captions),"max_display_tail_seconds":max(c["display_tail_seconds"] for c in captions),
             "diagram_text_bounds":all_bounds,"fonts":{str(p):sha(p) for p in [visuals.FONT,visuals.BOLD]},
             "files":{p.name:{"bytes":p.stat().st_size,"sha256":sha(p)} for p in output.iterdir() if p.is_file()},
             "limits":["No speech time scaling, trimming, resynthesis or within-scene rearrangement", "Provider timestamps are not independently measured phoneme timing", "Waveform checks do not establish pronunciation", "Full listening review remains separate", "No upload/publication"]}
    save(output/"preparation.json",receipt)
    return receipt


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs",required=True,type=Path);parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args(); os.nice(10)
    result=run(args.inputs.resolve(),args.output.resolve())
    print(json.dumps({k:result[k] for k in ["status","captions","narrated_words","raw_scene_audio_seconds","added_silence_seconds","peak_dBFS"]}))
