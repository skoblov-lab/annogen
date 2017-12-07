# distutils: language=c++
# cython: language_level=3

from typing import Dict, Tuple, Union, Iterable, NamedTuple, Mapping
from numbers import Integral, Real

from libc.stdint cimport uint8_t, int32_t, uint32_t, uint64_t
from libcpp cimport bool as cbool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport pair


SUPPORTED_TYPES = [str, Real, Integral]
Value = Union[str, Real, Integral]     # values can be either str, float or int
Annotations = Mapping[str, Value]  # feature ID -> feature value
GenLocus = NamedTuple("Locus", [("chrom", str), ("pos", Integral),
                                ("ref", str), ("alt", str)])


cdef extern from "mapping.hpp":

    cdef cppclass StringRecs:
        uint8_t first
        vector[string] second
        StringRecs()
        StringRecs(uint8_t, vector[string])

    cdef cppclass FloatRecs:
        uint8_t first
        vector[float] second
        FloatRecs()
        FloatRecs(uint8_t, vector[float])

    cdef cppclass IntRecs:
        uint8_t first
        vector[int32_t] second
        IntRecs()
        IntRecs(uint8_t, vector[int32_t])

    cdef cppclass Locus:
        uint8_t chrom
        uint32_t pos
        char ref
        char alt
        Locus() except +
        Locus(uint8_t chrom, uint32_t pos, char ref)
        Locus(uint8_t chrom, uint32_t pos, char ref, char alt)
        cbool operator==(const Locus& other) const

    cdef cppclass Records:
        vector[StringRecs] strings
        vector[FloatRecs] floats
        vector[IntRecs] integers
        Records() except +
        Records(const vector[StringRecs]& strings,
                const vector[FloatRecs]& floats,
                const vector[IntRecs]& integers) except +
        cbool operatorbool() const

    cdef cppclass LocusTable:
        cbool contains(const Locus& key) const
        Records& operator[](const Locus& key)
        uint64_t size()

    cdef cppclass StringCache:
        StringCache()
        int32_t size()
        int32_t add(const string& entry) except +
        int32_t find(const string& entry) except +
        string find(int32_t entry_code) except +

# TODO explicitly ask for data types
# TODO add cached_strings
# TODO rewrite encoder and decoder

cdef class GenomeMapping:

    cdef:
        LocusTable mapping
        dict _key_codes
        dict _key_codes_inv
        dict _chrom_codes
        dict _chrom_codes_inv
        dict _base_codes
        dict _base_codes_inv

    def __init__(self, features: Iterable[str],
                 rows: Iterable[Tuple[GenLocus, Annotations]]):
        self._chrom_codes = {}
        self._base_codes = {"": 0}  # default empty base to 0
        self._key_codes = {key: i for i, key in enumerate(features)}

        cdef:
            str chrom, ref, alt
            int pos, chrom_code, ref_code, alt_code
            Locus locus
            Records records

        for (chrom, pos, ref, alt), annotations in rows:
            # add new chromosome and base codes if needed
            chrom_code = self._chrom_codes.setdefault(chrom, len(self._chrom_codes))
            ref_code = self._base_codes.setdefault(ref, len(self._base_codes))
            alt_code = self._base_codes.setdefault(alt, len(self._base_codes))
            strings, numbers = self._encode_values(annotations)
            locus = Locus(chrom_code, pos, ref_code, alt_code)
            records = convert_records(strings, numbers)
            self.mapping[locus] = records

        self._key_codes_inv = {i: key for key, i in self._key_codes.items()}
        self._chrom_codes_inv = {i: chrom for chrom, i in self._chrom_codes.items()}
        self._base_codes_inv = {i: base for base, i in self._base_codes.items()}

    @property
    def keys(self):
        return self._key_codes

    @property
    def bases(self):
        return self._base_codes

    @property
    def chromosomes(self):
        return self._chrom_codes

    cpdef dict getitem(self, str chrom, int pos, str ref, str alt=""):
        # Return an empty dict if any of (chrom, ref, alt) have not been
        # encountered database construction
        if not (chrom in self._chrom_codes and ref in self._base_codes
                and alt in self._base_codes):
            return {}
        cdef:
            int chrom_code = self._chrom_codes[chrom]
            int ref_code = self._base_codes[ref]
            int alt_code = self._base_codes[alt]
            Locus locus = Locus(chrom_code, pos, ref_code, alt_code)
        if not self.mapping.contains(locus):
            return {}
        return self.decode_values(self.mapping[locus])

    def getitems(self, positions: Iterable):
        return [self.getitem(*position) for position in positions]

    cdef dict decode_values(self, Records records):
        cdef:
            dict decoded = {}
            int feat_id
            str feat
            bytes stringval
            double num_val

        for stringrec in records.strings:
            feat_id = stringrec.first
            feat = self._key_codes_inv[feat_id]
            stringval = stringrec.second
            decoded[feat] = stringval.decode()

        for numrec in records.numbers:
            feat_id = numrec.first
            feat = self._key_codes_inv[feat_id]
            num_val = numrec.second
            decoded[feat] = num_val

        return decoded

    cdef tuple _encode_values(self, annotations: Mapping[str, Value]):
        # remove unsuported types
        cdef:
            list strings = []
            list numbers = []
        try:
            for key, value in annotations.items():
                if isinstance(value, unicode):
                    strings.append((self._key_codes[key], value.encode()))
                elif isinstance(value, Real):
                    numbers.append((self._key_codes[key], float(value)))
                else:
                    raise ValueError("Unsupported type for feature {}".format(key))
            return strings, numbers
        except KeyError:
            raise KeyError("Encountered unspecified feature: {}".format(key))


cdef Records convert_records(list strings, list numbers):
    cdef:
        vector[StringRec] strings_
        vector[RealRec] numbers_
        int feat_id
        string str_val
        double num_val
    for feat_id, str_val in strings:
        strings_.push_back(StringRec(feat_id, str_val))
    for feat_id, num_val in numbers:
        numbers_.push_back(RealRec(feat_id, num_val))

    return Records(strings_, numbers_)
