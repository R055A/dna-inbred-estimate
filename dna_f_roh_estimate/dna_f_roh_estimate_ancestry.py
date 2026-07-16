from dna_f_roh_estimate.dna_f_roh_enum import (
    DnaKitFileHeaderAncestry,
    DnaKitFileInvalidAllelesSymbol,
)
from dna_f_roh_estimate.dna_f_roh_estimate_abstract_txt import InbredEstimateAbstractTxt


class InbredEstimateAncestry(InbredEstimateAbstractTxt):
    """
    Estimate F_ROH score, representing inbreeding/shared ancestry coefficient score, from Ancestry.com DNA kit data.
    """

    __INVALID_ALLELES: set[DnaKitFileInvalidAllelesSymbol] = {
        DnaKitFileInvalidAllelesSymbol.ZERO
    }  # There are likely more invalid entries yet to be discovered in more files

    def __init__(self, file_name: str | None = None) -> None:
        """
        Class constructor for instantiating an F_ROH estimator instance for Ancestry DNA kits with optional file name.
        :param file_name: Path name to a supported DNA kit data .txt file. If provided, file is read at instantiation.
        """
        super().__init__(
            file_name=file_name,
            headers=[h.value for h in DnaKitFileHeaderAncestry],
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
        alleles_one: str = (
            data_row[header[DnaKitFileHeaderAncestry.HEADER_ALLELE_ONE.value]]
            .strip()
            .upper()
        )
        alleles_two: str = (
            data_row[header[DnaKitFileHeaderAncestry.HEADER_ALLELE_TWO.value]]
            .strip()
            .upper()
        )
        if not self._is_valid_alleles(
            alleles_data=alleles_one
        ) or not self._is_valid_alleles(alleles_data=alleles_two):
            return [None, None]
        return [alleles_one, alleles_two]
