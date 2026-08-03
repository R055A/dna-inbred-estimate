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


class DnaKitFileHeaderFamilyTreeDna(Enum):
    HEADER_RESULT = "RESULT"


class DnaKitFileHeader23AndMe(Enum):
    HEADER_GENOTYPE = "GENOTYPE"


class DnaKitFileUniqueIdentifier(Enum):
    ANCESTRY = "AncestryDNA"
    MY_HERITAGE = "MyHeritage"
    FAMILY_TREE_DNA = "Family Tree DNA"
    TWENTY_THREE_AND_ME = "23andMe"


class DnaKitFileInvalidAllelesSymbol(Enum):
    ZERO = "0"
    DASH = "--"
