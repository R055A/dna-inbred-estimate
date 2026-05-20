from dna_f_roh_estimate.dna_f_roh_estimate_abstract import InbredEstimateAbstract


class InbredEstimateAncestry(InbredEstimateAbstract):
    """
    Estimate F_ROH score, representing inbreeding/shared ancestry coefficient score, from Ancestry.com DNA kit data.
    """

    __HEADER_ALLELE_ONE: str = "ALLELE1"
    __HEADER_ALLELE_TWO: str = "ALLELE2"

    def __init__(self, file_name: str | None = None) -> None:
        """
        Class constructor for instantiating an F_ROH estimator instance for Ancestry DNA kits with optional file name.
        :param file_name: Path name to a supported DNA kit data .txt file. If provided, file is read at instantiation.
        """
        super().__init__(file_name=file_name, headers=[self.__HEADER_ALLELE_ONE, self.__HEADER_ALLELE_TWO])

    def _split_data_row(self, data_row: str) -> list[str]:
        """
        Split a data row read from .txt DNA kit file into column values.
        :param data_row: Row of data from DNA kit file.
        :return: List of column values split from the data row.
        """
        return data_row.split("\t")

    def _iter_data_rows(self, header: dict[str, int], data_rows: list[str]):
        """
        Yield valid parsed data rows read from .txt DNA kit file.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_rows: List of data rows from the DNA kit file.
        :return: Iterator of parsed data rows for F_ROH analysis.
        """
        for data_row in data_rows:
            if data_row.strip():
                data_row = self._split_data_row(data_row=data_row)
                if len(data_row) > max(header.values()):
                    yield data_row

    def _get_alleles_data(self, header: dict[str, int], data_row: list[str]) -> list[str | None]:
        """
        Get allele values from a parsed data row read from .txt DNA kit file.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: List containing allele pair, or None values if the genotype is invalid.
        """
        alleles_one: str = data_row[header[self.__HEADER_ALLELE_ONE]].strip().upper()
        alleles_two: str = data_row[header[self.__HEADER_ALLELE_TWO]].strip().upper()
        if not self._is_valid_alleles(alleles_data=alleles_one) or not self._is_valid_alleles(alleles_data=alleles_two):
            return [None, None]
        return [alleles_one, alleles_two]
