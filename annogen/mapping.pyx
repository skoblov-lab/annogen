# distutils: language=c++
# cython: language_level=3, c_string_type=unicode, c_string_encoding=utf8

from typing import Dict, Tuple, Union, Iterable, NamedTuple, Mapping, List
from numbers import Integral, Real

from libc.stdint cimport uint8_t, int32_t, uint32_t, uint64_t
from libcpp cimport bool as cbool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.pair cimport pair


cdef:
    frozenset SUPPORTED_TYPES = frozenset([str, int, float])
    char MAXBASES = 127
    char MAXCONTIGS = 127


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
        vector[pair[uint8_t, vector[string]]] strings
        vector[pair[uint8_t, vector[float]]] floats
        vector[pair[uint8_t, vector[int32_t]]] integers
        Records() except +
        Records(const vector[pair[uint8_t, vector[string]]]& strings,
                const vector[pair[uint8_t, vector[float]]]& floats,
                const vector[pair[uint8_t, vector[int32_t]]]) except +
        cbool operatorbool() const

    cdef cppclass LocusTable:
        cbool contains(const Locus& key) const
        Records& operator[](const Locus& key)
        uint64_t size()

    cdef cppclass StringCache:
        StringCache()
        int32_t size()
        int32_t cache(const string& entry) except +
        string cache(int32_t entry_code) except +
        const vector[string]& cache()

 
Site = Tuple[str, int, str, str]  #  chrom, pos, ref, alt
# TODO explicitly ask for data types
# TODO add cached_strings
# TODO rewrite encoder and decoder 

cdef class GenomeMapping:

    cdef:
        LocusTable mapping
        StringCache stringcache
        set _cached
        list _features
        dict _feature_ids
        dict _dtypes
        list _contigs
        dict _contig_ids
        list _bases
        dict _base_ids

    def __init__(self,
                 features: Mapping[str, type],
                 contigs: Iterable[str],
                 alphabet: Iterable[str],
                 cached_strings: Iterable[str],
                 entries: Iterable[Site, Dict[str, List]]):
        """
        :param features: a mapping from feature to type: int, float or
        str; actual data are not required to be of these types, but they must
        be convertible (i.e. calling any of these constructors will not raise
        an error).
        :param entries: entries to insert
        :param cached_strings: cached strings; be warned, caching introduces a
        significant overhead and thus should only be used for really long
        repetitive strings (avoiding caching for strings shorter than 10
        characters and occurring less than 5 times is a rational rule of thumb);
        """
        # initialise features
        self._dtypes = dict(features)
        self._features = list(self._dtypes)
        self._feature_ids = {feat: i for i, feat in enumerate(self._features)}
        if any(dtype not in SUPPORTED_TYPES for dtype in features.values()):
            raise TypeError('dtype not in {}'.format(SUPPORTED_TYPES))
        # initialise contigs
        self._contigs = list(contigs)
        if len(self._contigs) > MAXCONTIGS:
            raise ValueError(f"There can be no more than {MAXCONTIGS} contigs")
        self._contig_ids = {cont: i for i, cont in enumerate(self._contigs)}
        # initialise alphabet
        self._bases = [''] + list(alphabet)  # default empty base
        if len(self._contigs) > MAXCONTIGS:
            raise ValueError(f"There can be no more than {MAXBASES} unique "
                             f"items in the alphabet, including the default "
                             f"empty base")
        self._base_ids = {base: i for i, base in enumerate(self._bases)}
        # make sure all bases and contigs are strings
        if not all(isinstance(val, str) for val in self._bases + self._contigs):
            raise ValueError("All contigs and alphabet items must be strings")
        # make sure all bases can be automatically converted to char
        if not all(len(base) <= 1 for base in self._bases):
            raise ValueError("alphabet can only be comprised of "
                             "single-character strings, except for the default "
                             "empty base string")
        # initialise the list of cached strings
        self._cached = set(cached_strings)
        if any(feature not in self._features for feature in self._cached):
            raise ValueError('`cached_strings` contains features absent in '
                             '`features`')
        if any(self._dtypes[f] is not str for f in self._cached):
            raise ValueError('only string values can be cached')
        for (contig, pos, ref, alt), annotations in entries:
            self.insert(contig, pos, ref, alt, annotations)

    def insert(self, str contig, int pos, str ref, str alt, dict annotations):
        cdef:
            uint8_t contig_code = self.ccode(contig)
            char ref_code = self.bcode(ref)
            char alt_code = self.bcode(alt)
            Locus locus = Locus(contig_code, pos, ref_code, alt_code)
            Records records = self.encode(annotations)
        self.mapping[locus] = records

    cpdef dict getitem(self, str contig, int pos, str ref, str alt):
        # Return an empty dict if any of (contig, ref, alt) have not been
        # indexed
        if (contig not in self._contig_ids or
                    ref not in self._base_ids or alt not in self._base_ids):
            return {}
        cdef:
            uint8_t contig_code = self.ccode(contig)
            char ref_code = self.bcode(ref)
            char alt_code = self.bcode(alt)
            Locus locus = Locus(contig_code, pos, ref_code, alt_code)
            Records records = self.mapping[locus]
        return self.decode(records)

    def getitems(self, positions: Iterable):
        return [self.getitem(*position) for position in positions]

    cdef inline int fcode(self, str feature):
        """
        Return a feature code
        :param feature: 
        :return: 
        """
        if feature not in self._feature_ids:
            raise KeyError(f'encountered unknown feature "{feature}"')
        return self._feature_ids[feature]

    cdef inline str feature(self, int feature_code):
        return self._features[feature_code]

    cdef inline int ccode(self, str contig):
        """
        Return a feature code
        :param contig: 
        :return: 
        """
        if contig not in self._contig_ids:
            raise KeyError(f'encountered unknown contig "{contig}"')
        return self._contig_ids[contig]

    cdef inline str contig(self, int contig_code):
        return self._contigs[contig_code]

    cdef inline int bcode(self, str base):
        """
        Return a base code
        :param base: 
        :return: 
        """
        if base not in self._base_ids:
            raise KeyError(f'encountered unknown feature "{base}"')
        return self._base_ids[base]

    cdef inline str base(self, int base_code):
        return self._bases[base_code]

    @property
    def features(self):
        return self._features

    @property
    def bases(self):
        return self._bases

    @property
    def contigs(self):
        return self._contigs

    cdef dict decode(self, Records records):
        cdef:
            vector[pair[uint8_t, vector[string]]] stringrecs = records.strings
            vector[pair[uint8_t, vector[float]]] floatrecs = records.floats
            vector[pair[uint8_t, vector[int32_t]]] intrecs = records.integers
            list float_features = [
                (self.feature(code), values) for code, values in floatrecs
            ]
            # exclude cached strings from intrecs
            list int_features = [
                (self.feature(code), values)
                for code, values in intrecs
                if self.feature(code) not in self._cached
            ]
            # include cached strings from intrecs
            list string_features = [
                (self.feature(code), self.frombytes(values))
                for code, values in stringrecs
            ] + [
                (self.feature(code), self.fromcache(values))
                for code, values in intrecs
                if self.feature(code) in self._cached
            ]
            dict decoded = {
                feature: values for feature, values in
                float_features + int_features + string_features
            }
        return decoded

    cdef Records encode(self, dict annotations):
        # separate values by type and cast to either int, str or float
        cdef:
            list strings = [
                (self.fcode(f), self.tobytes(values))
                for f, values in annotations.items()
                if self._dtypes[f] is str and f not in self._cached
            ]
            list cached = [
                (self.fcode(f), self.tocache(values))
                for f, values in annotations.items()
                if self._dtypes[f] is str and f in self._cached
            ]
            list floats = [
                (self.fcode(f), self.cast(float, values))
                for f, values in annotations.items()
                if self._dtypes[f] is float
            ]
            # integer records include cached string records
            list ints = [
                (self.fcode(f), self.cast(int, values))
                for f, values in annotations.items()
                if self._dtypes[f] is int
            ] + cached

            vector[pair[uint8_t, vector[string]]] stringrecs = strings
            vector[pair[uint8_t, vector[float]]] floatrecs = floats
            vector[pair[uint8_t, vector[int32_t]]] intrecs = ints

        return Records(stringrecs, floatrecs, intrecs)

    cdef list tocache(self, list strings):
        """
        Cache
        :param strings: 
        :return: 
        """
        cdef:
            list cached = []
            string s
        for s in self.tobytes(self.cast(str, strings)):
            cached.append(self.stringcache.cache(s))
        return cached

    cdef list fromcache(self, list cached_strings):
        cdef:
            list uncached = []
            int cached
        for cached in cached_strings:
            uncached.append(self.stringcache.cache(cached))
        return self.frombytes(uncached)

    cdef inline list tobytes(self, list unicode_strings):
        return [s.encode() if not isinstance(s, bytes) else s
                for s in unicode_strings]

    cdef inline list frombytes(self, list byte_stings):
        return [s.decode() if isinstance(s, bytes) else s
                for s in byte_stings]

    cdef inline list cast(self, type constructor, list values):
        return [constructor(val) for val in values]

# TODO add explicit type conversion
# TODO add annotation for entry type
# TODO improve docs
