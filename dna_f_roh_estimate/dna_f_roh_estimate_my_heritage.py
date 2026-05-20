from dna_f_roh_estimate.dna_f_roh_estimate_abstract import InbredEstimateAbstract
from csv import reader


class InbredEstimateMyHeritage(InbredEstimateAbstract):
    """
    Estimate F_ROH score, representing inbreeding/shared ancestry coefficient score, from MyHeritage.com DNA kit data.
    """

    __HEADER_RESULT: str = "RESULT"

    def __init__(self, file_name: str | None = None) -> None:
        """
        Class constructor for instantiating an F_ROH estimator instance for MyHeritage DNA kits with optional file name.
        :param file_name: Path name to a supported DNA kit data .csv file. If provided, file is read at instantiation.
        """
        super().__init__(file_name=file_name, headers=[self.__HEADER_RESULT])

    def _split_data_row(self, data_row: str) -> list[str]:
        """
        Split a data row read from .csv DNA kit file into column values.
        :param data_row: Row of data from DNA kit file.
        :return: List of column values split from the data row.
        """
        return next(reader([data_row]))

    def _iter_data_rows(self, header: dict[str, int], data_rows: list[str]):
        """
        Yield valid parsed data rows read from .csv DNA kit file.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_rows: List of data rows from the DNA kit file.
        :return: Iterator of parsed data rows for F_ROH analysis.
        """
        for data_row in reader(data_rows):
            if data_row and len(data_row) > max(header.values()):
                yield data_row

    def _get_alleles_data(self, header: dict[str, int], data_row: list[str]) -> list[str | None]:
        """
        Get allele values from a parsed data row read from .csv DNA kit file.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: List containing allele pair, or None values if the genotype is invalid.
        """
        result: str = data_row[header[self.__HEADER_RESULT]].strip().upper()
        if not self._is_valid_alleles(alleles_data=result):
            return [None, None]
        return list(result) + [None]
