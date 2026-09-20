"""Generate only the declared narration requests; credential read without echo, never saved."""
import argparse
import base64
import datetime
import getpass
import hashlib
import json
from pathlib import Path
import time
import urllib.error
import urllib.request
import wave
import warnings


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(path, data):
    with path.open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False)
        f.write('\n')


def synthesize(key, payload, output):
    output.mkdir(exist_ok=False)
    save(output / 'request.json', payload)
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    tick = time.monotonic()
    req = urllib.request.Request(
        'https://api.cartesia.ai/tts/sse',
        data=json.dumps(payload).encode(),
        headers={'Authorization': 'Bearer ' + key, 'Cartesia-Version': '2026-08-14',
                 'Content-Type': 'application/json', 'Accept': 'text/event-stream'},
    )
    chunks, events, words, starts, ends = [], [], [], [], []
    done = False
    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            assert response.status == 200
            assert 'text/event-stream' in response.headers.get('Content-Type', '')
            safe_headers = {k: v for k, v in response.headers.items()
                            if k.lower() in {'content-type', 'date', 'x-request-id', 'request-id'}}
            for raw in response:
                if time.monotonic() - tick > 120:
                    raise RuntimeError('Per-request time limit reached; do not automatically retry')
                line = raw.decode().strip()
                if not line or not line.startswith('data:'):
                    continue
                event = json.loads(line[5:].strip())
                if done:
                    raise ValueError('Events after completion')
                if event.get('type') == 'error':
                    raise RuntimeError('Provider error: ' + json.dumps(event))
                assert event.get('status_code', 500) in (200, 206)
                if event['type'] == 'chunk':
                    data = base64.b64decode(event['data'], validate=True)
                    chunks.append(data)
                    events.append({k: v for k, v in event.items() if k != 'data'} |
                                  {'audio_bytes': len(data), 'audio_sha256': hashlib.sha256(data).hexdigest()})
                elif event['type'] == 'timestamps':
                    ts = event['word_timestamps']
                    assert len(ts['words']) == len(ts['start']) == len(ts['end'])
                    words.extend(ts['words']); starts.extend(ts['start']); ends.extend(ts['end'])
                    events.append(event)
                elif event['type'] == 'done':
                    assert event.get('done') is True
                    done = True
                    events.append(event)
                else:
                    raise ValueError('Unexpected SSE event type')
        assert done and chunks and words
        pcm = b''.join(chunks)
        assert len(pcm) % 2 == 0 and any(pcm)
        duration = len(pcm) / 2 / 44100
        assert 0 < duration < 90
        for i, (a, b) in enumerate(zip(starts, ends)):
            assert 0 <= a <= b <= duration + 0.1
            if i:
                assert a >= starts[i - 1] and b >= ends[i - 1]
        audio = output / 'narration.wav'
        with wave.open(str(audio), 'wb') as f:
            f.setnchannels(1); f.setsampwidth(2); f.setframerate(44100); f.writeframes(pcm)
        save(output / 'events.json', events)
        save(output / 'word-timestamps.json', {'words': words, 'start': starts, 'end': ends})
        receipt = {'status': 'complete', 'started_at_utc': started,
                   'completed_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                   'elapsed_seconds': time.monotonic() - tick,
                   'endpoint': 'https://api.cartesia.ai/tts/sse', 'api_version': '2026-08-14',
                   'requested_model': payload['model_id'], 'resolved_snapshot': None,
                   'voice_id': payload['voice'], 'request_sha256': sha(output / 'request.json'),
                   'response_headers': safe_headers, 'duration_seconds': duration,
                   'samples': len(pcm) // 2, 'sample_rate': 44100, 'channels': 1, 'sample_width': 2,
                   'timestamp_words': len(words), 'files': {p.name: sha(p) for p in output.iterdir() if p.is_file()},
                   'credential_saved': False, 'human_listening_review': 'pending',
                   'plan_and_training_settings': 'Not verified for Cartesia; do not inherit Claude/Codex settings.'}
        save(output / 'receipt.json', receipt)
        return {'directory': str(output), 'status': 'complete', 'seconds': duration, 'words': len(words)}
    except Exception as error:
        message = str(error).replace(key, '[REDACTED]')
        if isinstance(error, urllib.error.HTTPError):
            message += ': ' + error.read().decode('utf-8', 'replace').replace(key, '[REDACTED]')[:1000]
        save(output / 'failure.json', {'status': 'failed', 'error_type': type(error).__name__,
                                     'message': message, 'automatic_retry': False,
                                     'started_at_utc': started, 'received_audio_chunks': len(chunks)})
        raise RuntimeError(message) from None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    manifest = json.loads(a.manifest.read_text())
    assert 1 <= len(manifest['requests']) <= 6
    assert not a.output.exists()
    for row in manifest['requests']:
        assert row['payload']['model_id'] == 'sonic-3.6'
        assert row['payload']['output_format'] == {'container': 'raw', 'encoding': 'pcm_s16le', 'sample_rate': 44100}
        assert row['payload']['add_timestamps'] is True
        assert row['name'].isalnum() or all(c.isalnum() or c in '-_' for c in row['name'])
    with warnings.catch_warnings():
        warnings.simplefilter('error', getpass.GetPassWarning)
        key = getpass.getpass('Cartesia credential (not saved): ')
    assert key.startswith('sk_car_') and '\n' not in key
    a.output.mkdir(parents=True)
    for row in manifest['requests']:
        print(json.dumps(synthesize(key, row['payload'], a.output / row['name'])), flush=True)
    key = None


if __name__ == '__main__':
    main()
