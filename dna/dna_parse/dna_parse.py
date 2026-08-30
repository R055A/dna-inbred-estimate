from dna.dna_parse.dna_parse_csv_family_tree_dna import ParseFileFamilyTreeDNACSV
from dna.dna_parse.dna_parse_csv_my_heritage import ParseFileMyHeritageCSV
from dna.dna_parse.dna_parse_csv_tellme_gen import ParseFileTellmeGenCSV
from dna.dna_parse.dna_parse_txt_23_and_me import ParseFile23AndMeTXT
from dna.dna_parse.dna_parse_txt_ancestry import ParseFileAncestryTXT
from dna.dna_parse.dna_parse_txt_living_dna import ParseFileLivingDNATXT
from dna.dna_parse.dna_parse_abstract import ParseFileDNA
from utils.file_io import read_data_file_first_five_lines
from dna.dna_enum import (
    DnaKitFileHeaderFamilyTreeDNA,
    DnaKitFileHeaderTellmeGen,
    DnaKitFileUniqueIdentifier,
)


def is_format_header_in_file(file_name: str, header: str) -> bool:
    """
    Conditional for if a genotype column header exists in the first rows of a data file.
    :param file_name: The name of the DNA data file.
    :param header: The genotype column header name to search for.
    :return: True if the header is found in the first rows of the data file; otherwise False.
    """
    return (
        header.upper() in read_data_file_first_five_lines(file_name=file_name).upper()
    )


def is_format_identifier_in_file(
    file_name: str, unique_format_identifier: DnaKitFileUniqueIdentifier
) -> bool:
    """
    Conditional for if a unique identifier exists in the first row of a data file to verify DNA kit format.
    :param file_name: The name of the DNA data file.
    :param unique_format_identifier: The unique identifier word(s) in the first row of the DNA data file.
    :return: True if the first row of a data file contains the unique identifier word(s); otherwise False.
    """
    return unique_format_identifier.value in read_data_file_first_five_lines(
        file_name=file_name
    )


def dna_f_parse(file_name: str) -> ParseFileDNA | None:
    """
    Factory function for instantiating an object extending ParseFileDNA with raw DNA kit data file passed as argument.
    :param file_name: The name (and type) of the DNA data file.
    :return: ParseFileDNA object with formatted DNA autosomal SNP data.
    """
    dna_f_parse_obj: ParseFileDNA | None = None
    if file_name.endswith(".csv"):
        if is_format_identifier_in_file(
            file_name=file_name,
            unique_format_identifier=DnaKitFileUniqueIdentifier.MY_HERITAGE,
        ) and not isinstance(dna_f_parse_obj, ParseFileMyHeritageCSV):
            dna_f_parse_obj = ParseFileMyHeritageCSV(file_name=file_name)
        elif is_format_header_in_file(
            file_name=file_name,
            header=DnaKitFileHeaderFamilyTreeDNA.HEADER_RESULT.value,
        ) and not isinstance(dna_f_parse_obj, ParseFileFamilyTreeDNACSV):
            dna_f_parse_obj = ParseFileFamilyTreeDNACSV(file_name=file_name)
        elif is_format_header_in_file(
            file_name=file_name,
            header=DnaKitFileHeaderTellmeGen.HEADER_GENOTYPE.value,
        ) and not isinstance(dna_f_parse_obj, ParseFileTellmeGenCSV):
            dna_f_parse_obj = ParseFileTellmeGenCSV(file_name=file_name)
        elif not isinstance(dna_f_parse_obj, ParseFileMyHeritageCSV):
            raise ValueError("CSV file format is not supported.")
    elif file_name.endswith(".txt"):
        if is_format_identifier_in_file(
            file_name=file_name,
            unique_format_identifier=DnaKitFileUniqueIdentifier.ANCESTRY,
        ) and not isinstance(dna_f_parse_obj, ParseFileAncestryTXT):
            dna_f_parse_obj = ParseFileAncestryTXT(file_name=file_name)
        elif is_format_identifier_in_file(
            file_name=file_name,
            unique_format_identifier=DnaKitFileUniqueIdentifier.LIVING_DNA,
        ) and not isinstance(dna_f_parse_obj, ParseFileLivingDNATXT):
            dna_f_parse_obj = ParseFileLivingDNATXT(file_name=file_name)
        elif is_format_identifier_in_file(
            file_name=file_name,
            unique_format_identifier=DnaKitFileUniqueIdentifier.TWENTY_THREE_AND_ME,
        ) and not isinstance(dna_f_parse_obj, ParseFile23AndMeTXT):
            dna_f_parse_obj = ParseFile23AndMeTXT(file_name=file_name)
        elif (
            not isinstance(dna_f_parse_obj, ParseFileAncestryTXT)
            and not isinstance(dna_f_parse_obj, ParseFile23AndMeTXT)
            and not isinstance(dna_f_parse_obj, ParseFileLivingDNATXT)
        ):
            raise ValueError("Text file format is not supported.")
    else:
        raise ValueError("File type is not supported. Use either .csv or .txt")
    return dna_f_parse_obj if dna_f_parse_obj else None
