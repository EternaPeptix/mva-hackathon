#!/usr/bin/env python3
"""Compile the pinned DELLY parser verbatim and expose its actual query intervals.

This is a native component probe using the same interval type and installed
Boost/HTSlib. Whole-caller exclusion behavior is tested separately on synthetic reads.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import time

UTIL_SHA = 'ae492aafea144c7746da230f60141aaca6a5a9ff78c27264aa5fb27686a77554'


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def extract(text, start):
    begin = text.index(start); opening = text.index('{', begin); depth = 1; end = opening + 1
    while depth:
        if text[end] == '{': depth += 1
        elif text[end] == '}': depth -= 1
        end += 1
    return text[begin:end]


def build(source, prefix, output, compiler):
    if sha(source / 'src/util.h') != UTIL_SHA:
        raise ValueError('unreviewed parser source')
    typedef = 'typedef boost::icl::interval_set<uint32_t> TChrIntervals;'
    if typedef not in (source / 'src/delly.h').read_text():
        raise ValueError('native interval type changed')
    output.mkdir(parents=True, exist_ok=False)
    text = (source / 'src/util.h').read_text()
    helper = extract(text, '  inline bool\n  is_gz(')
    function = extract(text, '  template<typename TConfig, typename TRegionsGenome>\n  inline int32_t\n  _parseExcludeIntervals(')
    includes = '''#include <cstdint>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <boost/filesystem/path.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <htslib/sam.h>
'''
    main = r'''
struct Config { bool hasExcludeFile; boost::filesystem::path exclude; };
int main(int argc, char** argv) {
    if (argc != 3) return 2;
    Config c; c.hasExcludeFile = std::string(argv[1]) != "NONE"; c.exclude = argv[1];
    std::ifstream fai(argv[2]); std::string line; std::ostringstream header;
    header << "@HD\tVN:1.6\n";
    while (std::getline(fai,line)) {
        std::istringstream fields(line); std::string name; uint32_t length;
        if (!(fields >> name >> length)) return 3;
        header << "@SQ\tSN:" << name << "\tLN:" << length << "\n";
    }
    std::string h = header.str(); bam_hdr_t* hdr = sam_hdr_parse(h.size(), h.c_str());
    if (!hdr) return 4;
    typedef boost::icl::interval_set<uint32_t> TChrIntervals;
    std::vector<TChrIntervals> regions;
    if (!_parseExcludeIntervals(c,hdr,regions)) return 5;
    for (int i = 0; i < hdr->n_targets; ++i)
        for (auto it = regions[i].begin(); it != regions[i].end(); ++it)
            std::cout << hdr->target_name[i] << '\t' << it->lower() << '\t' << it->upper() << '\n';
    bam_hdr_destroy(hdr); return 0;
}
'''
    cpp = output / 'parser.cpp'; cpp.write_text(includes + helper + '\n' + function + '\n' + main)
    binary = output / 'parser'
    args = [str(compiler), '-std=c++17', '-O2', '-I' + str(prefix / 'include'), str(cpp),
            '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib'),
            '-lboost_iostreams', '-lboost_filesystem', '-lboost_system', '-lhts', '-o', str(binary)]
    tick = time.monotonic(); run = subprocess.run(args, capture_output=True, text=True)
    result = {'built_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source_util_sha256': UTIL_SHA, 'source_delly_header_sha256': sha(source / 'src/delly.h'),
              'generator_source_sha256': sha(Path(__file__)), 'extracted_helper_sha256': hashlib.sha256(helper.encode()).hexdigest(),
              'extracted_parser_sha256': hashlib.sha256(function.encode()).hexdigest(),
              'generated_cpp_sha256': sha(cpp), 'command': args, 'exit_code': run.returncode,
              'elapsed_seconds': time.monotonic() - tick, 'stdout': run.stdout, 'stderr': run.stderr,
              'compiler_version': subprocess.run([str(compiler), '--version'], capture_output=True, text=True, check=True).stdout,
              'scope': __doc__}
    if run.returncode == 0: result['binary_sha256'] = sha(binary)
    (output / 'build.json').write_text(json.dumps(result, indent=2) + '\n')
    if run.returncode: raise RuntimeError('native parser probe build failed; see retained build.json')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'prefix', 'output'): p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--compiler', type=Path, default=Path('/usr/bin/clang++'))
    a = p.parse_args(); r = build(a.source, a.prefix, a.output, a.compiler)
    print(json.dumps({key:r[key] for key in ('built_at','exit_code','elapsed_seconds','binary_sha256')}, indent=2))


if __name__ == '__main__': main()
