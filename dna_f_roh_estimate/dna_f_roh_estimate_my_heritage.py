from dna_f_roh_estimate.dna_f_roh_estimate_abstract_csv import InbredEstimateAbstractCsv
from dna_f_roh_estimate.dna_f_roh_enum import DnaKitFileHeaderMyHeritage


class InbredEstimateMyHeritage(InbredEstimateAbstractCsv):
    """
    Estimate F_ROH score, representing inbreeding/shared ancestry coefficient score, from MyHeritage.com DNA kit data.
    """

    __INVALID_ALLELES: tuple[str] = {"--"}  # There are likely more invalid entries yet to be discovered in more files

    def __init__(self, file_name: str | None = None) -> None:
        """
        Class constructor for instantiating an F_ROH estimator instance for MyHeritage DNA kits with optional file name.
        :param file_name: Path name to a supported DNA kit data .csv file. If provided, file is read at instantiation.
        """
        super().__init__(
            file_name=file_name,
            headers=[h.value for h in DnaKitFileHeaderMyHeritage],
            invalid_alleles=self.__INVALID_ALLELES
        )

    def _get_alleles_data(self, header: dict[str, int], data_row: list[str]) -> list[str | None]:
        """
        Get allele values from a parsed data row read from .csv DNA kit file.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: List containing allele pair, or None values if the genotype is invalid.
        """
        result: str = data_row[header[DnaKitFileHeaderMyHeritage.HEADER_RESULT.value]].strip().upper()
        if not self._is_valid_alleles(alleles_data=result):
            return [None, None]
        return list(result) + [None]
