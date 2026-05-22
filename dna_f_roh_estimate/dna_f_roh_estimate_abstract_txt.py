from dna_f_roh_estimate.dna_f_roh_estimate_abstract import InbredEstimateAbstract
from abc import abstractmethod


class InbredEstimateAbstractTxt(InbredEstimateAbstract):
    """
    Estimate F_ROH score, representing inbreeding/shared ancestry coefficient score, from .txt format DNA kit data.
    """

    def __init__(
            self,
            file_name: str | None = None,
            headers: list[str] | None = None,
            invalid_alleles: tuple[str] | None = None
    ) -> None:
        """
        Class constructor for instantiating an F_ROH estimator instance for .txt format DNA kit data.
        :param file_name: Path name to a supported DNA kit data .txt file. If provided, file is read at instantiation.
        :param headers: List of format-specific genotype header names in addition to [RSID, CHROMOSOME, POSITION].
        :param invalid_alleles: Allele file entries which are not valid for F_ROH analysis.
        """
        super().__init__(file_name=file_name, headers=headers, invalid_alleles=invalid_alleles)

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

    @abstractmethod
    def _get_alleles_data(self, header: dict[str, int], data_row: list[str]) -> list[str | None]:
        """
        Get allele values from a parsed data row read from .txt DNA kit file.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: List containing allele pair, or None values if the genotype is invalid.
        """
        pass
