from dna.dna_enum import DnaKitFileInvalidAllelesSymbol, DnaKitFileHeaderTellmeGen
from dna.dna_parse.dna_parse_txt import ParseFileTXT


class ParseFileTellmeGenCSV(ParseFileTXT):
    """
    Parse raw DNA data to DataDNA dataclass format from tellmeGen DNA kit .csv file. Extends ParseFileTXT.

    tellmeGen raw data files use a .csv extension but are tab-delimited.
    """

    __INVALID_ALLELES: set[DnaKitFileInvalidAllelesSymbol] = {
        DnaKitFileInvalidAllelesSymbol.DASH
    }

    def __init__(self, file_name: str | None = None) -> None:
        """
        Class constructor with optional file_name parameter.
        :param file_name: Path name to a supported DNA kit data .csv file. If provided, file is read at instantiation.
        """
        super().__init__(
            file_name=file_name,
            headers=[h.value for h in DnaKitFileHeaderTellmeGen],
            invalid_alleles=self.__INVALID_ALLELES,
        )

    def _split_data_row(self, data_row: str) -> list[str]:
        """
        Split a data row read from tellmeGen .csv DNA kit file into column values.
        :param data_row: Row of data from DNA kit file.
        :return: List of column values split from the data row.
        """
        data_row = data_row.lstrip("#").strip()
        if "\t" in data_row:
            return data_row.split("\t")
        return data_row.split()

    def _get_alleles_data(
        self, header: dict[str, int], data_row: list[str]
    ) -> list[str | None]:
        """
        Get allele values from a parsed data row read from .csv DNA kit file.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: List containing allele pair, or None values if the genotype is invalid.
        """
        result: str = (
            data_row[header[DnaKitFileHeaderTellmeGen.HEADER_GENOTYPE.value]]
            .strip()
            .upper()
        )
        if not self._is_valid_alleles(alleles_data=result):
            return [None, None]
        return list[str](result) + [None]
