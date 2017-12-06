#ifndef mapping_h
#define mapping_h


#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <string>
#include <utility>
#include "sparsepp/spp.h"


struct Locus {
    // Genomic Locus
    uint8_t chrom;
    uint32_t pos;
    char ref;
    char alt;

    Locus():
        chrom(0), pos(0), ref(0), alt(0) {}
    Locus(uint8_t chrom, uint32_t pos, char ref, char alt):
        chrom(chrom), pos(pos), ref(ref), alt(alt) {}
    Locus(uint8_t chrom, uint32_t pos, char ref):
        chrom(chrom), pos(pos), ref(ref), alt(0) {}

    bool operator==(const Locus& other) const {
        return (chrom == other.chrom &&
                pos == other.pos &&
                ref == other.ref &&
                alt == other.alt);
    }
};


namespace std {
    // inject specialization of std::hash for Locus into namespace std
    template<>
    struct hash<Locus> {
        std::size_t operator()(const Locus& locus) const {
            std::size_t seed = 0;
            spp::hash_combine(seed, locus.chrom);
            spp::hash_combine(seed, locus.pos);
            spp::hash_combine(seed, locus.ref);
            spp::hash_combine(seed, locus.alt);
            return seed;
        }
    };
}


// Each record type is a pair of feature identifier and feature value;
// Since feature IDs are stored as uint8_t, there can be up to 256 unique IDs
typedef std::pair<uint8_t , std::string> StringRec; // string record
typedef std::pair<uint8_t, double> RealRec;         // number record


struct Records {
    // Records associated with a Locus.
    std::vector<StringRec> strings;
    std::vector<RealRec> numbers;

    Records():
        strings(0), numbers(0) {}
    Records(const std::vector<StringRec> & strings,
            const std::vector<RealRec> & numbers):
        strings(strings), numbers(numbers) {}

    explicit operator bool() const {
        return strings.size() || numbers.size();
    }
};


typedef spp::sparse_hash_map<Locus, Records> LocusTable;


//class LocusTable {
//    spp::sparse_hash_map<Locus, Records> table;
//
//    LocusTable() {}
//
//    void setitem(const Locus& locus, const Records& records) {
//        if (table.contains(locus)) {
//            throw std::invalid_argument("Can't replace an existing key");
//        }
//        table[locus] = records;
//    }


#endif
