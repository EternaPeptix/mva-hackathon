#!/usr/bin/env python3
"""Finite owned-group guard for fictional qualification or explicitly requested public audit."""
import argparse
import collections
import datetime
import hashlib
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path('/Users/jeweled/Documents/mva-hackathon')
NAME = 'population-frequency-all-source-audit-20260920'
METHOD = ROOT / 'methods' / NAME
RESULT = ROOT / 'results' / NAME
LIMITS = dict(wall_seconds=600, group_RSS_bytes=4*1024**3,
              output_bytes=64*1024**2, minimum_free_bytes=40*1024**3)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def save(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def fingerprint(path, expected=None):
    path = Path(path).resolve()
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            digest.update(block)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'input changed while hashing')
    result = dict(path=str(path), bytes=after.st_size, mtime_ns=after.st_mtime_ns, sha256=digest.hexdigest())
    if expected:
        require(result['sha256'] == expected['sha256'] and result['bytes'] == expected['bytes'], 'input pin mismatch: '+str(path))
    return result


def process_group(pgid):
    result = []
    for line in subprocess.check_output(['ps', '-axo', 'pid=,pgid=,rss='], text=True).splitlines():
        values = line.split()
        if len(values) == 3 and int(values[1]) == pgid:
            result.append(dict(pid=int(values[0]), rss_bytes=int(values[2])*1024))
    return result


def terminate_owned(process):
    require(process.pid != os.getpgrp(), 'refuse parent group')
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    until = time.monotonic()+5
    while time.monotonic() < until:
        process.poll()
        if not process_group(process.pid):
            break
        time.sleep(.1)
    if process_group(process.pid):
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    process.wait(timeout=5)


def total_size(directory):
    return sum(p.stat().st_size for p in directory.rglob('*') if p.is_file())


def caller_snapshot():
    p = subprocess.run(['ps', '-p', '79898', '-o', 'command='], text=True, capture_output=True)
    return dict(pid=79898, present=p.returncode == 0, command=p.stdout.strip() or None,
                interpretation='Read-only snapshot. The audit never signals this process; only its own new session can be terminated.')


def read_config(path, expected_sha):
    pin = fingerprint(path)
    require(pin['sha256'] == expected_sha, 'configuration SHA mismatch')
    config = json.loads(Path(path).read_text())
    require(config['scope'] == 'read-only 25-exact-key public frequency field audit', 'scope changed')
    require(config['root'] == str(ROOT) and config['limits'] == LIMITS, 'root/limits changed')
    require(platform.machine() == 'arm64' and platform.node().split('.')[0] == '512S2', '512S2 only')
    require(config['expected_key_count'] == 25 and config['expected_fictional_controls'] == 18, 'scope/control count changed')
    require(Path(config['database']['path']) == ROOT/'exomiser/2512_hg38/2512_hg38_variants.mv.db', 'database path changed')
    require(set(config['sources']) == {'AllFrequencySources.java', 'run-audit.py'}, 'source set changed')
    require(set(config['inputs']) == {'queries.tsv', 'selected-keys.json', 'runtime-pins.json'}, 'input set changed')
    return config, pin


def worker(args):
    os.nice(10)
    out = args.output.resolve()
    config, config_pin = read_config(args.config, args.expected_config_sha256)
    pins = [config_pin]
    try:
        for name, expected in (config['sources'] | config['inputs']).items():
            pins.append(fingerprint(METHOD/name, expected))
        runtime = json.loads((METHOD/'runtime-pins.json').read_text())
        require(Path(sys.executable).resolve() == Path(runtime['runtime']['python']['path']), 'Python runtime changed')
        for expected in list(runtime['runtime'].values()) + runtime['dependencies']:
            pins.append(fingerprint(expected['path'], expected))
        selected = json.loads((METHOD/'selected-keys.json').read_text())
        keys = [r['key'] for r in selected['records']]
        query_keys = [line.replace('\t', ':') for line in (METHOD/'queries.tsv').read_text().splitlines()]
        require(len(keys) == len(set(keys)) == 25 and keys == query_keys, 'exact key set/order mismatch')
        prefix = [runtime['runtime']['java']['path'], '-Xmx2g', '-XX:ActiveProcessorCount=2',
                  '-cp', runtime['runtime']['lirical_jar']['path'], str(METHOD/'AllFrequencySources.java')]
        if args.worker_mode == 'qualify':
            command = prefix + ['--self-test']
        else:
            qualification = fingerprint(args.qualification, config['qualification_receipt'])
            receipt = json.loads(args.qualification.read_text())
            require(receipt['status'] == 'complete' and receipt['mode'] == 'qualify', 'qualification incomplete')
            require(receipt['source_sha256'] == {n: p['sha256'] for n, p in config['sources'].items()}, 'qualification source mismatch')
            require(receipt['qualified_config_core'] == {k: v for k, v in config.items() if k != 'qualification_receipt'}, 'qualified configuration changed')
            require(receipt['controls']['tests_passed'] == 18 and receipt['controls']['tests_failed'] == 0
                    and receipt['controls']['tests_skipped'] == 0 and receipt['controls']['database_opened'] is False, 'qualification did not pass')
            pins.append(qualification)
            pins.append(fingerprint(config['database']['path'], config['database']))
            command = prefix + ['--database', config['database']['path'], '--queries', str(METHOD/'queries.tsv'), '--output', str(out/'all-frequency-sources.json')]
        save(out/'command.json', dict(command=command, mode=args.worker_mode, started_at=stamp(), inputs=pins))
        with (out/'java.stdout.log').open('xb') as stdout, (out/'java.stderr.log').open('xb') as stderr:
            completed = subprocess.run(command, stdout=stdout, stderr=stderr)
        require(completed.returncode == 0, 'Java exit '+str(completed.returncode))
        result = dict(status='complete', completed_at=stamp(), mode=args.worker_mode,
                      source_sha256={n: p['sha256'] for n,p in config['sources'].items()}, inputs=pins,
                      qualified_config_core={k: v for k, v in config.items() if k != 'qualification_receipt'})
        if args.worker_mode == 'qualify':
            controls = json.loads((out/'java.stdout.log').read_text())
            require(controls['status'] == 'complete' and controls['tests_passed'] == 18 and controls['tests_failed'] == 0
                    and controls['tests_skipped'] == 0 and controls['database_opened'] is False, 'fictional controls failed')
            result['controls'] = controls
        else:
            data = json.loads((out/'all-frequency-sources.json').read_text())
            require(data['status'] == 'complete' and [r['key'] for r in data['records']] == keys, 'audit key/output mismatch')
            for actual, old in zip(data['records'], selected['records']):
                require(actual['current_probe_projection'] == old['prior_probe_response'], 'current-prefix projection differs from saved response: '+old['key'])
                require(actual['stored_frequency_row_count'] == len(actual['all_stored_frequency_rows_in_order']), 'row count mismatch')
            result['explanation_counts'] = dict(collections.Counter(r['empty_export_explanation'] for r in data['records']))
            result['all_25_saved_prefix_projections_exact'] = True
            result['scope'] = 'Resource-field explanation only; original inventories and frequency classifications unchanged'
        # Bind runtime origins to the freshly hashed dependency set; no unpinned reader class.
        data_for_origins = controls if args.worker_mode == 'qualify' else data
        from urllib.parse import unquote, urlparse
        pinned_paths = {p['path'] for p in pins}
        for url in data_for_origins['class_origins'].values():
            parsed = urlparse(url)
            require(parsed.scheme == 'file' and str(Path(unquote(parsed.path)).resolve()) in pinned_paths, 'unbound reader class origin')
        for old in pins:
            current = fingerprint(old['path'], old)
            require(current['mtime_ns'] == old['mtime_ns'], 'input modified during job')
        result['outputs'] = [fingerprint(p) for p in sorted(out.iterdir()) if p.is_file() and p.name not in ['guard.stdout.log', 'guard.stderr.log']]
        save(out/'receipt.json', result)
        save(out/'status.json', dict(status='complete', completed_at=stamp(), receipt_sha256=fingerprint(out/'receipt.json')['sha256']))
    except BaseException as error:
        save(out/'status.json', dict(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error), retained_without_retry=True))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--expected-config-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--qualify', action='store_true')
    mode.add_argument('--execute', action='store_true')
    mode.add_argument('--worker-mode', choices=['qualify', 'query'])
    parser.add_argument('--qualification', type=Path)
    args = parser.parse_args()
    if args.worker_mode:
        worker(args)
        return
    config, _ = read_config(args.config, args.expected_config_sha256)
    out = args.output.resolve()
    require(out.is_relative_to(METHOD/'qualifications' if args.qualify else RESULT), 'output outside new owned scope')
    require(not out.exists(), 'refuse overwrite or automatic retry')
    if args.execute:
        require(args.qualification is not None and config.get('qualification_receipt') is not None, 'reviewed qualification binding required')
    require(shutil.disk_usage(ROOT).free >= LIMITS['minimum_free_bytes'], 'insufficient free disk')
    out.mkdir(parents=True)
    command = [sys.executable, '-B', str(Path(__file__).resolve()), '--config', str(args.config),
               '--expected-config-sha256', args.expected_config_sha256, '--output', str(out),
               '--worker-mode', 'qualify' if args.qualify else 'query']
    if args.qualification:
        command += ['--qualification', str(args.qualification)]
    start = time.monotonic()
    state = dict(status='running', started_at=stamp(), command=command, limits=LIMITS, caller_before=caller_snapshot())
    process = None
    try:
        with (out/'guard.stdout.log').open('xb') as stdout, (out/'guard.stderr.log').open('xb') as stderr:
            process = subprocess.Popen(command, stdout=stdout, stderr=stderr, start_new_session=True)
            state['owned_group'] = process.pid
            peak = 0
            while True:
                group = process_group(process.pid)
                peak = max(peak, sum(p['rss_bytes'] for p in group))
                require(time.monotonic()-start < LIMITS['wall_seconds'], 'whole-job time limit')
                require(sum(p['rss_bytes'] for p in group) < LIMITS['group_RSS_bytes'], 'whole-job group RSS limit')
                require(total_size(out) < LIMITS['output_bytes'], 'whole-job output limit')
                require(shutil.disk_usage(ROOT).free >= LIMITS['minimum_free_bytes'], 'minimum free disk')
                code = process.poll()
                if code is not None:
                    require(not process_group(process.pid), 'owned descendant left running')
                    break
                time.sleep(.25)
        require(code == 0, 'guarded worker failed')
        terminal = json.loads((out/'status.json').read_text())
        receipt = fingerprint(out/'receipt.json')
        require(terminal['status'] == 'complete' and terminal['receipt_sha256'] == receipt['sha256'], 'worker terminal linkage failed')
        state.update(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic()-start,
                     peak_sampled_group_RSS_bytes=peak, receipt_sha256=receipt['sha256'])
    except BaseException as error:
        if process is not None:
            terminate_owned(process)
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error), retained_without_retry=True)
        raise
    finally:
        state['caller_after'] = caller_snapshot()
        save(out/'guard-status.json', state)


if __name__ == '__main__':
    main()
