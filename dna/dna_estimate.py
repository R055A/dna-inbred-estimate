from dna.dna_f_roh.dna_f_roh_plink import EstimatePlinkFROH
from dna.dna_data import DataDNA, DataFROH, DataConfigFROH
from dna.dna_parse.dna_parse_abstract import ParseFileDNA
from dna.dna_parse.dna_parse import dna_f_parse
from typing import Iterator


class EstimateDNA:
    """
    Estimate F_ROH score from DNA kit data analysis using PLINK.
    """

    def __init__(self, file_name: str | None = None) -> None:
        """
        Class constructor with optional file name parameter.
        :param file_name: Path name to a supported DNA kit data file. If provided, the file is read at instantiation.
        """
        self.__dna_f_roh: EstimatePlinkFROH | None = None
        self.__dna_f_parse: ParseFileDNA | None = None
        if file_name:
            self.__dna_f_parse = dna_f_parse(file_name=file_name)

    def __repr__(self) -> str:
        """
        :return: F_ROH analysis summary data in string format, or an empty string if no data is available.
        """
        if self.__dna_f_roh and self.__dna_f_roh.data:
            return self.__dna_f_roh.data.__repr__()
        return ""

    def __iter__(self) -> Iterator[tuple[str, int | float | None]]:
        """
        :return: F_ROH analysis summary data in dictionary format, or an empty dictionary if no data is available.
        """
        if self.__dna_f_roh and self.__dna_f_roh.data:
            yield from self.__dna_f_roh.data

    @property
    def f_roh_data(self) -> DataFROH | None:
        """
        :return: Data from F_ROH analysis, or None if no data is available.
        """
        if self.__dna_f_roh and self.__dna_f_roh.data:
            return self.__dna_f_roh.data
        return None

    @property
    def dna_data(self) -> DataDNA | None:
        """
        :return: Parsed DNA autosomal SNP data in DataDNA dataclass format, or None if no data is available.
        """
        if self.__dna_f_parse:
            return self.__dna_f_parse.dna_data
        return None

    @dna_data.setter
    def dna_data(self, dna_data: DataDNA) -> None:
        """
        Set parsed DNA autosomal SNP data to DataDNA dataclass format.
        :param dna_data: parsed SNP data in DataDNA dataclass format.
        """
        if self.__dna_f_parse:
            self.__dna_f_parse.dna_data = dna_data

    def read_dna_data(self, file_name: str) -> DataDNA | None:
        """
        Read data from DNA kit file and parse autosomal SNP data to DataDNA dataclass format.
        :param file_name: The name of the DNA data file.
        :return: Parsed DNA autosomal SNP data in DataDNA dataclass format, or None if no data is available.
        """
        dna_f_parse_tmp: ParseFileDNA | None = dna_f_parse(file_name=file_name)
        if dna_f_parse_tmp:
            self.__dna_f_parse = dna_f_parse_tmp
            return self.__dna_f_parse.read_dna_data(file_name=file_name)
        return None

    def compute_f_roh(self, config: DataConfigFROH) -> None:
        """
        Compute F_ROH from parsed autosomal DNA data file for estimating inbreeding/shared-ancestry coefficient score.
        :param config: hyperparameters for computing F_ROH in DataConfigFROH dataclass format.
        """
        if self.dna_data:
            self.__dna_f_roh = EstimatePlinkFROH(dna_data=self.dna_data)
            self.__dna_f_roh.compute_f_roh(config=config)


def dna_estimate_f_roh(args) -> None:
    dna_inbred_estimate: EstimateDNA = EstimateDNA(file_name=args.file)
    dna_inbred_estimate.compute_f_roh(
        config=DataConfigFROH(*list(vars(args).copy().values())[1:])
    )
    print(dna_inbred_estimate)
