#ifndef mapping_h
#define mapping_h

#include <cinttypes>
#include <cstddef>
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
typedef std::pair<uint8_t , std::vector<std::string>> StringRecs;  // string records
typedef std::pair<uint8_t, std::vector<float>> FloatRecs;          // float records
typedef std::pair<uint8_t, std::vector<int32_t>> IntRecs;          // int records


struct Records {
    // Records associated with a Locus.
    std::vector<StringRecs> strings;
    std::vector<RealRecs> floats;
    std::vector<FloatRecs> integers;

    Records():
        strings(0), floats(0),  integers(0) {}
    Records(const std::vector<StringRecs> & strings,
            const std::vector<FloatRecs> & floats,
            const std::vector<IntRecs> & integers):
        strings(strings), floats(floats), integers(integers) {}

    explicit operator bool() const {
        return strings.size() || floats.size() || integers.size();
    }
};


typedef spp::sparse_hash_map<Locus, Records> LocusTable;

class StringCache {

private:

    spp::sparse_hash_map<std::string, int32_t> cachemap;
    std::vector<std::string> strings;
    int32_t sizelimit = INT32_MAX;

public:
    StringCache(): cachemap(0), strings(0) {}

    int32_t size() {
        return strings.size();
    }

    int32_t add(const std::string& entry) {
        if (size() == sizelimit) {
            throw std::runtime_error("Exceeded the cache size limit");
        }
        strings.push_back(entry);
        position = strings.size() - 1;
        cachemap[entry] = position;
        return position;
    }

    int32_t find(const std::string& entry) {
        // Return ID for a string
        if (!cachemap.contains(entry)) {
            throw std::invalid_argument("No such entry");
        }
        return cachemap[entry];
    }
    std::string find(int32_t entry_code) {
        // Return string for an ID
        // note: although returning a const reference seems more efficient,
        //       the reference might get invalidated by future calls to `add`,
        //       because the `strings` vector might get reallocated.
        if (!strings.size() > entry_code) {
            throw std::invalid_argument("No such entry");
        }
        return strings[entry_code];
    }
};


#endif
