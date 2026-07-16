from dna_f_roh_estimate.dna_f_roh_enum import (
    DnaKitFileHeader23AndMe,
    DnaKitFileInvalidAllelesSymbol,
)
from dna_f_roh_estimate.dna_f_roh_estimate_abstract_txt import InbredEstimateAbstractTxt


class InbredEstimate23AndMe(InbredEstimateAbstractTxt):
    """
    Estimate F_ROH score, representing inbreeding/shared ancestry coefficient score, from 23AndMe.com DNA kit data.
    """

    __INVALID_ALLELES: set[DnaKitFileInvalidAllelesSymbol] = {
        DnaKitFileInvalidAllelesSymbol.ZERO
    }  # There are likely more invalid entries yet to be discovered in more files

    def __init__(self, file_name: str | None = None) -> None:
        """
        Class constructor for instantiating an F_ROH estimator instance for 23AndMe DNA kits with optional file name.
        :param file_name: Path name to a supported DNA kit data .txt file. If provided, file is read at instantiation.
        """
        super().__init__(
            file_name=file_name,
            headers=[h.value for h in DnaKitFileHeader23AndMe],
            invalid_alleles=self.__INVALID_ALLELES,
        )

    def _get_alleles_data(
        self, header: dict[str, int], data_row: list[str]
    ) -> list[str | None]:
        """
        Get allele values from a parsed data row read from .txt DNA kit file.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: List containing allele pair, or None values if the genotype is invalid.
        """
        result: str = (
            data_row[header[DnaKitFileHeader23AndMe.HEADER_GENOTYPE.value]]
            .strip()
            .upper()
        )
        if not self._is_valid_alleles(alleles_data=result):
            return [None, None]
        return list(result) + [None]
