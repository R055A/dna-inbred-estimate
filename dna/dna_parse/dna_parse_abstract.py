from dna.dna_enum import (
    DnaKitFileInvalidAllelesSymbol,
    DnaKitFileHeaderAbstract,
)
from utils.file_io import read_data_file_all_lines
from dna.dna_data import DataDNA, DataGenotype
from abc import abstractmethod, ABC


class ParseFileDNA(ABC):
    """
    Abstract class for parsing raw DNA data to DataDNA dataclass format from DNA kit file.
    """

    def __init__(
        self,
        file_name: str | None = None,
        headers: list[str] | None = None,
        invalid_alleles: set[DnaKitFileInvalidAllelesSymbol] | None = None,
    ) -> None:
        """
        Class constructor with optional file_name, headers and invalid_alleles parameters.
        :param file_name: Path name to a supported DNA kit data file. If provided, file is read at instantiation.
        :param headers: List of format-specific genotype header names in addition to [RSID, CHROMOSOME, POSITION].
        :param invalid_alleles: Allele file entries which are not valid for F_ROH analysis.
        """
        self.__invalid_alleles: set[DnaKitFileInvalidAllelesSymbol] = (
            invalid_alleles if invalid_alleles else set()
        )
        self._headers: list[str] = [h.value for h in DnaKitFileHeaderAbstract] + (
            headers or []
        )

        self.__dna_data: DataDNA | None = None
        if file_name:
            self.read_dna_data(file_name=file_name)

    @property
    def dna_data(self) -> DataDNA | None:
        """
        :return: Parsed DNA autosomal SNP data in DataDNA dataclass format, or None if no data is available.
        """
        return self.__dna_data

    @dna_data.setter
    def dna_data(self, dna_data: DataDNA) -> None:
        """
        Set parsed DNA autosomal SNP data to DataDNA dataclass format.
        :param dna_data: parsed SNP data in DataDNA dataclass format.
        """
        self.__dna_data = dna_data

    @abstractmethod
    def _split_data_row(self, data_row: str) -> list[str]:
        """
        Split a data row into column values relevant to a subclass file format.
        :param data_row: Row of data from the DNA kit file.
        :return: List of column values split from the data row.
        """
        pass

    @abstractmethod
    def _iter_data_rows(self, header: dict[str, int], data_rows: list[str]):
        """
        Yield valid parsed data rows relevant to a subclass file format.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_rows: List of data rows from the DNA kit file.
        :return: Iterator of parsed data rows for F_ROH analysis.
        """
        pass

    @abstractmethod
    def _get_alleles_data(
        self, header: dict[str, int], data_row: list[str]
    ) -> list[str | None]:
        """
        Get allele values from a parsed data row relevant to a subclass file format.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: List containing allele pair, or None values if the genotype is invalid.
        """
        pass

    def _is_valid_alleles(self, alleles_data: str) -> bool:
        """
        Conditional for if an allele value is valid for F_ROH analysis.
        :param alleles_data: Allele value from the DNA kit file.
        :return: True if the value is valid; otherwise False.
        """
        return alleles_data not in self.__invalid_alleles

    def __parse_dna_data(self, header: dict[str, int], data_rows: list[str]) -> DataDNA:
        """
        Parse raw DNA kit file data for F_ROH analysis of autosomal SNP genotypes in DataDNA dataclass format.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_rows: List of raw data rows from the DNA kit file.
        :return: Parsed DNA autosomal SNP data in DataDNA dataclass format.
        """
        genotypes: list[DataGenotype] = []
        invalid_cnt: int = 0

        for data_row in self._iter_data_rows(header=header, data_rows=data_rows):
            chromosome = (
                data_row[header[DnaKitFileHeaderAbstract.HEADER_CHROMOSOME.value]]
                .strip()
                .upper()
            )
            if not chromosome or not chromosome.isdigit():
                continue
            chromosome = int(chromosome)
            if not (1 <= int(chromosome) <= 22):
                continue

            pos = data_row[
                header[DnaKitFileHeaderAbstract.HEADER_POSITION.value]
            ].strip()
            if not pos or not pos.isdigit():
                continue
            pos = int(pos)

            allele_one, allele_two = self._get_alleles_data(
                header=header, data_row=data_row
            )[:2]
            if not allele_one or not allele_two:
                invalid_cnt += 1
                continue

            genotypes.append(
                DataGenotype(
                    rs_id=data_row[
                        header[DnaKitFileHeaderAbstract.HEADER_RSID.value]
                    ].strip(),
                    chromosome=chromosome,
                    position=pos,
                    allele_one=allele_one,
                    allele_two=allele_two,
                )
            )
        genotypes.sort(
            key=lambda genotype: (
                genotype.chromosome,
                genotype.position,
            )
        )
        return DataDNA(
            genotypes=tuple[DataGenotype](genotypes),
            invalid_cnt=invalid_cnt,
        )

    def read_dna_data(self, file_name: str) -> DataDNA | None:
        """
        Parse raw DNA kit file data to DataDNA dataclass format for analysis of autosomal SNP genotypes.
        :param file_name: Path name to a supported DNA kit data file.
        :return: Parsed DNA autosomal SNP data in DataDNA dataclass format.
        """
        file_lines: list[str] = read_data_file_all_lines(file_name=file_name)
        for i, line in enumerate[str](file_lines):
            cols = [
                c.strip("#").strip().upper()
                for c in self._split_data_row(data_row=line)
            ]
            if all(col.upper() in cols for col in self._headers):
                self.__dna_data = self.__parse_dna_data(
                    header={
                        col.upper(): cols.index(col.upper()) for col in self._headers
                    },
                    data_rows=file_lines[i + 1 :],
                )
                return self.__dna_data
        raise ValueError("Data format is not supported.")
