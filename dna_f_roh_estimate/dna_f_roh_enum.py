from enum import Enum


class DnaKitFileHeaderAbstract(Enum):
    HEADER_RSID = "RSID"
    HEADER_CHROMOSOME = "CHROMOSOME"
    HEADER_POSITION = "POSITION"


class DnaKitFileHeaderAncestry(Enum):
    HEADER_ALLELE_ONE = "ALLELE1"
    HEADER_ALLELE_TWO = "ALLELE2"


class DnaKitFileHeaderMyHeritage(Enum):
    HEADER_RESULT = "RESULT"
