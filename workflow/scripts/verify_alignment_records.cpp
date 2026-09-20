// Compare every decoded record without writing read names or sequences to logs.
// Strict equality includes all 16 FLAG bits and auxiliary type/value bytes.
// Auxiliary tag order and unused packed-sequence padding are not biological data.
#include <htslib/sam.h>
#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

static void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

static void tags(const bam1_t *record, std::vector<std::string_view> &result) {
    result.clear();
    errno = 0;
    auto current = bam_aux_first(record);
    if (!current) require(errno == ENOENT, "invalid auxiliary data");
    while (current) {
        errno = 0;
        auto next = bam_aux_next(record, current);
        if (!next) require(errno == ENOENT, "invalid auxiliary data");
        const char *start = reinterpret_cast<const char *>(current - 2);
        const char *end = next ? reinterpret_cast<const char *>(next - 2)
                               : reinterpret_cast<const char *>(record->data + record->l_data);
        result.emplace_back(start, end - start);
        current = next;
    }
    std::sort(result.begin(), result.end(), [](auto a, auto b) {
        return a.substr(0, 2) < b.substr(0, 2);
    });
    for (size_t i = 1; i < result.size(); ++i)
        require(result[i - 1].substr(0, 2) != result[i].substr(0, 2), "duplicate auxiliary tag");
}

int main(int argc, char **argv) {
    uint64_t count = 0, duplicate = 0, secondary = 0, supplementary = 0, unmapped = 0;
    samFile *left = nullptr, *right = nullptr;
    sam_hdr_t *lh = nullptr, *rh = nullptr;
    bam1_t *a = nullptr, *b = nullptr;
    int exit_code = 1;
    try {
        require(argc == 6, "usage: verifier reference left right expected_records threads");
        const uint64_t expected = std::stoull(argv[4]);
        const int threads = std::stoi(argv[5]);
        require(threads >= 1 && threads <= 4 && expected > 0, "invalid limits");
        left = sam_open(argv[2], "r"); right = sam_open(argv[3], "r");
        require(left && right, "cannot open alignment");
        require(hts_set_fai_filename(left, argv[1]) == 0 &&
                hts_set_fai_filename(right, argv[1]) == 0, "cannot set reference");
        require(hts_set_threads(left, threads) == 0 && hts_set_threads(right, threads) == 0,
                "cannot configure decoding threads");
        lh = sam_hdr_read(left); rh = sam_hdr_read(right);
        require(lh && rh, "cannot decode header");
        require(sam_hdr_nref(lh) == sam_hdr_nref(rh), "reference count differs");
        for (int i = 0; i < sam_hdr_nref(lh); ++i)
            require(std::strcmp(sam_hdr_tid2name(lh, i), sam_hdr_tid2name(rh, i)) == 0 &&
                    sam_hdr_tid2len(lh, i) == sam_hdr_tid2len(rh, i), "reference dictionary differs");
        // Header metadata is separately checked by the Python preparation gate.
        a = bam_init1(); b = bam_init1();
        require(a && b, "cannot allocate records");
        std::vector<std::string_view> at, bt;
        while (true) {
            const int ar = sam_read1(left, lh, a), br = sam_read1(right, rh, b);
            require(ar >= -1 && br >= -1, "alignment decoding error");
            if (ar == -1 || br == -1) {
                require(ar == br, "alignment record counts differ");
                break;
            }
            const auto &x = a->core; const auto &y = b->core;
            require(std::strcmp(bam_get_qname(a), bam_get_qname(b)) == 0, "read name/order differs");
            require(x.flag == y.flag, "full FLAG field differs");
            require(x.tid == y.tid && x.pos == y.pos && x.qual == y.qual &&
                    x.mtid == y.mtid && x.mpos == y.mpos && x.isize == y.isize,
                    "coordinate, mapping quality or mate field differs");
            require(x.n_cigar == y.n_cigar &&
                    std::memcmp(bam_get_cigar(a), bam_get_cigar(b), x.n_cigar * sizeof(uint32_t)) == 0,
                    "CIGAR differs");
            require(x.l_qseq == y.l_qseq, "read length differs");
            require(std::memcmp(bam_get_seq(a), bam_get_seq(b), x.l_qseq / 2) == 0,
                    "read sequence differs");
            if (x.l_qseq % 2)
                require(bam_seqi(bam_get_seq(a), x.l_qseq - 1) ==
                        bam_seqi(bam_get_seq(b), y.l_qseq - 1), "read sequence differs");
            require(std::memcmp(bam_get_qual(a), bam_get_qual(b), x.l_qseq) == 0, "base quality differs");
            tags(a, at); tags(b, bt);
            require(at == bt, "auxiliary tag type or value differs");
            ++count;
            duplicate += !!(x.flag & BAM_FDUP);
            secondary += !!(x.flag & BAM_FSECONDARY);
            supplementary += !!(x.flag & BAM_FSUPPLEMENTARY);
            unmapped += !!(x.flag & BAM_FUNMAP);
            if (count % 10000000 == 0) std::cerr << "records_compared=" << count << '\n';
        }
        require(count == expected, "total differs from expected source QC");
        std::cout << "{\"records_compared\":" << count
                  << ",\"duplicate_records\":" << duplicate
                  << ",\"secondary_records\":" << secondary
                  << ",\"supplementary_records\":" << supplementary
                  << ",\"unmapped_records\":" << unmapped
                  << ",\"full_flags_equal\":true,\"all_auxiliary_types_values_equal\":true,"
                     "\"absolute_record_order_equal\":true,\"all_compared_fields_equal\":true}\n";
        exit_code = 0;
    } catch (const std::exception &error) {
        std::cerr << "Verification failed after " << count << " records: " << error.what()
                  << "; read identifiers and sequence omitted.\n";
    }
    if (a) bam_destroy1(a);
    if (b) bam_destroy1(b);
    if (lh) sam_hdr_destroy(lh);
    if (rh) sam_hdr_destroy(rh);
    if (left && sam_close(left) != 0) exit_code = 1;
    if (right && sam_close(right) != 0) exit_code = 1;
    return exit_code;
}
