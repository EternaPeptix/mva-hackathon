#!/usr/bin/env python3
"""Export three committed Markdown drafts and their unchanged prediction CSV.

Usage: python3 workflow/scripts/export_submission_reports.py NEW_OUTPUT_DIR
       [--revision HEAD]

This local, stdlib-only export does not upload or publish anything. Relative file
links become GitHub blob links pinned to the selected commit; private repository
access is still required. Only the current reports' inline, untitled links with
single-token destinations are supported (labels may span lines). Unsupported
reference-style, image and HTML links fail explicitly, rather than being omitted.
"""

import argparse
import hashlib
import os
from pathlib import Path, PurePosixPath
import posixpath
import re
import subprocess
from urllib.parse import quote, unquote, urlsplit


REPORTS = (
    'reports/jeweled_track1_report_submission.md',
    'reports/jeweled_track2_report_submission.md',
    'reports/jeweled_methods_form_answers.md',
)
CSV = 'reports/jeweled_evidence-reconciled-v2.csv'
INLINE_LINK = re.compile(r'(?<!!)\[[^\[\]]+\]\((?P<target>[^\s()<>]+)\)')
SCHEME = re.compile(r'^[A-Za-z][A-Za-z0-9+.-]*:')


def git(root, *args):
    result = subprocess.run(
        ['git', '-C', str(root), '--literal-pathspecs', *args],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if result.returncode:
        raise ValueError('Git could not read the requested repository metadata or committed file: '
                         + result.stderr.decode('utf-8', errors='replace').strip())
    return result.stdout


def github_origin(root):
    origins = git(root, 'remote', 'get-url', '--all', 'origin').decode().splitlines()
    if len(origins) != 1:
        raise ValueError('Exactly one known GitHub origin URL is required')
    origin = origins[0]
    if origin.startswith('git@github.com:'):
        path = origin[len('git@github.com:'):]
    else:
        url = urlsplit(origin)
        if (url.hostname != 'github.com' or url.port is not None or url.query or url.fragment
                or url.password or not (
                    (url.scheme == 'https' and url.username is None)
                    or (url.scheme == 'ssh' and url.username == 'git'))):
            raise ValueError('Origin must be an unambiguous HTTPS or SSH github.com repository URL')
        path = url.path.lstrip('/')
    path = path.removesuffix('.git')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9-]*/[A-Za-z0-9_.-]+', path):
        raise ValueError('Origin must identify exactly one GitHub owner/repository')
    if path.split('/')[1] in ('.', '..'):
        raise ValueError('Invalid GitHub repository name')
    return 'https://github.com/' + path


def require_file(root, revision, path, checked):
    if path in checked:
        return
    entry = git(root, 'ls-tree', '-z', revision, '--', path)
    rows = entry.rstrip(b'\0').split(b'\0') if entry else []
    if len(rows) != 1:
        raise ValueError('Referenced file is absent or ambiguous in the selected commit: ' + path)
    metadata, name = rows[0].split(b'\t', 1)
    mode, kind, _ = metadata.split()
    if (name.decode('utf-8') != path or kind != b'blob'
            or mode not in (b'100644', b'100755')):
        raise ValueError('Referenced path is not a regular committed file: ' + path)
    checked.add(path)


def check_link_syntax(text, source):
    # These fixed documents contain no link-bearing code examples. Refuse them
    # if later introduced, so literal examples cannot silently become live links.
    # Bracketed genotype notation alone is not a Markdown link.
    for code in re.finditer(r'(`+)(.*?)\1|(?ms:^\s*~~~.*?^\s*~~~\s*$)', text, re.DOTALL):
        if '](' in code.group() or re.search(r'\[[^]]*\]\s*(?:\[[^]]*\]|:)', code.group()):
            raise ValueError(source + ': link-like syntax inside code is unsupported')
    unsupported = (
        (r'!\[', 'image links'),
        (r'(?m)^\s{0,3}\[[^\]\n]+\]:', 'reference-style link definitions'),
        (r'\[[^\[\]]*\]\s*\[[^\[\]]*\]', 'reference-style links'),
        (r'<[^>\n]*\b(?:href|src)\s*=|<\s*/?\s*(?:a|img|link|source|iframe|object|embed|video|audio)\b',
         'HTML links'),
    )
    for pattern, label in unsupported:
        if re.search(pattern, text, re.IGNORECASE):
            raise ValueError(source + ': unsupported ' + label)
    remainder = INLINE_LINK.sub('', text)
    if '](' in remainder:
        raise ValueError(source + ': unsupported inline link syntax (titles, escapes or nested destinations)')
    if re.search(r'\\\[|\\\]\(', text):
        raise ValueError(source + ': escaped link-like syntax is unsupported')


def render_report(data, source, root, revision, origin, checked):
    text = data.decode('utf-8')
    check_link_syntax(text, source)
    counts = {'relative_files': 0, 'same_document': 0, 'external': 0}

    def replace(match):
        target = match.group('target')
        if target.startswith('#'):
            counts['same_document'] += 1
            return match.group()
        if target.lower().startswith('file:') or re.match(r'^[A-Za-z]:[\\/]', target):
            raise ValueError(source + ': local file URI or drive path is unsupported')
        if target.startswith('//') or SCHEME.match(target):
            counts['external'] += 1
            return match.group()
        path, marker, fragment = target.partition('#')
        if '?' in path or '\\' in path or path.startswith('/'):
            raise ValueError(source + ': unsupported local link destination: ' + target)
        decoded = unquote(path, encoding='utf-8', errors='strict')
        resolved = posixpath.normpath(posixpath.join(str(PurePosixPath(source).parent), decoded))
        if (not decoded or decoded.startswith('/') or '\\' in decoded
                or resolved == '..' or resolved.startswith('../') or '\0' in resolved):
            raise ValueError(source + ': local link escapes the repository or is empty: ' + target)
        require_file(root, revision, resolved, checked)
        if resolved == source and marker:
            destination = '#' + fragment
            counts['same_document'] += 1
        else:
            destination = origin + '/blob/' + revision + '/' + quote(resolved, safe='/')
            destination += marker + fragment
            counts['relative_files'] += 1
        start, end = match.span('target')
        relative_start, relative_end = start - match.start(), end - match.start()
        return match.group()[:relative_start] + destination + match.group()[relative_end:]

    rendered = INLINE_LINK.sub(replace, text)
    header = ('> **Private draft — not submitted.** Source: [' + revision + ']('
              + origin + '/blob/' + revision + '/' + quote(source, safe='/')
              + '). Repository access is required.\n\n')
    return (header + rendered).encode('utf-8'), counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='New directory outside the source checkout; parent must exist')
    parser.add_argument('--revision', default='HEAD', help='Committed Git revision to export (default: HEAD)')
    args = parser.parse_args()
    try:
        root = Path(git(Path(__file__).resolve().parent, 'rev-parse', '--show-toplevel').decode().strip()).resolve()
        if os.path.lexists(args.output.expanduser()):
            raise ValueError('Existing destination is preserved; choose a new output directory')
        output = args.output.expanduser().resolve()
        if output == root or root in output.parents:
            raise ValueError('Output must be outside the source checkout')
        if not output.parent.is_dir():
            raise ValueError('Output parent directory must already exist')
        revision = git(root, 'rev-parse', '--verify', '--end-of-options',
                       args.revision + '^{commit}').decode().strip()
        if not re.fullmatch(r'[0-9a-f]{40,64}', revision):
            raise ValueError('Git did not resolve an unambiguous full commit ID')
        origin = github_origin(root)
        checked, rendered, counts = set(), {}, {}
        for source in (*REPORTS, CSV):
            require_file(root, revision, source, checked)
            data = git(root, 'show', revision + ':' + source)
            name = PurePosixPath(source).name
            if source in REPORTS:
                rendered[name], counts[name] = render_report(data, source, root, revision, origin, checked)
            else:
                rendered[name] = data
        # Nothing is created until all four source files and every link validate.
        output.mkdir(mode=0o700)
        for name, data in rendered.items():
            with (output / name).open('xb') as stream:
                stream.write(data)
        print('Exported four private draft files from commit ' + revision + ' to ' + str(output))
        for name, count in counts.items():
            print(name + ': ' + str(count['relative_files']) + ' file links pinned; '
                  + str(count['same_document']) + ' same-document and '
                  + str(count['external']) + ' external links preserved')
        print(PurePosixPath(CSV).name + ': byte-identical CSV; SHA256 '
              + hashlib.sha256(rendered[PurePosixPath(CSV).name]).hexdigest())
    except (OSError, ValueError, UnicodeError) as error:
        parser.exit(2, 'Export failed: ' + str(error) + '\n')


if __name__ == '__main__':
    main()
