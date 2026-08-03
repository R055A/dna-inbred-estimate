from dna_f_roh_estimate.dna_f_roh_estimate_23_and_me import InbredEstimate23AndMe
from dna_f_roh_estimate.dna_f_roh_estimate_my_heritage import InbredEstimateMyHeritage
from dna_f_roh_estimate.dna_f_roh_estimate_family_tree_dna import (
    InbredEstimateFamilyTreeDna,
)
from dna_f_roh_estimate.dna_f_roh_estimate_ancestry import InbredEstimateAncestry
from dna_f_roh_estimate.dna_f_roh_estimate_abstract import InbredEstimateAbstract
from dna_f_roh_estimate.dna_f_roh_enum import DnaKitFileUniqueIdentifier
from dna_f_roh_estimate.dna_f_roh_data import DataFROH
from utils import read_data_file_first_five_lines
from typing import Iterator


class InbredEstimate:
    """
    Estimate F_ROH score, representing inbreeding/shared ancestry coefficient score, from DNA kit data.
    """

    def __init__(self, file_name: str | None = None) -> None:
        """
        Class constructor for instantiating an F_ROH estimator instance with optional file name parameter.
        :param file_name: Path name to a supported DNA kit data file. If provided, the file is read at instantiation.
        """
        self.__dna_inbred_estimate: InbredEstimateAbstract | None = None
        if file_name:
            self.read_dna_data_file(file_name=file_name)

    def __repr__(self) -> str:
        """
        :return: F_ROH analysis summary data in string format, or an empty string if no data is available.
        """
        if self.__dna_inbred_estimate and self.__dna_inbred_estimate.f_roh_data:
            return self.__dna_inbred_estimate.f_roh_data.__repr__()
        return ""

    def __iter__(self) -> Iterator[tuple[str, int | float | None]]:
        """
        :return: F_ROH analysis summary data in dictionary format, or an empty dictionary if no data is available.
        """
        if self.__dna_inbred_estimate and self.__dna_inbred_estimate.f_roh_data:
            yield from self.__dna_inbred_estimate.f_roh_data

    @property
    def f_roh_data(self) -> DataFROH | None:
        """
        :return: Data from F_ROH analysis.
        """
        if self.__dna_inbred_estimate and self.__dna_inbred_estimate.f_roh_data:
            return self.__dna_inbred_estimate.f_roh_data
        return None

    @property
    def snp_data(self) -> list[tuple[int, int, bool]]:
        """
        :return: Parsed DNA autosomal SNP data formatted into tuples of chromosome, base-pair position, homozygosity:
        [
            {
                autosomal_chromosome: int,
                base_pair_position: int,
                is_homozygous: bool
            }
        ]
        """
        if self.__dna_inbred_estimate:
            return self.__dna_inbred_estimate.snp_data
        return []

    @snp_data.setter
    def snp_data(self, snp_data: list[tuple[int, int, bool]]) -> None:
        """
        Set parsed DNA autosomal SNP data for F_ROH analysis.
        :param snp_data: parsed SNP data formatted into tuples of chromosome, base-pair position, homozygosity:
        [
            {
                autosomal_chromosome: int,
                base_pair_position: int,
                is_homozygous: bool
            }
        ]
        """
        if self.__dna_inbred_estimate:
            self.__dna_inbred_estimate.snp_data = snp_data

    def __is_format_identifier_in_file(
        self, file_name: str, unique_format_identifier: DnaKitFileUniqueIdentifier
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

    def read_dna_data_file(self, file_name: str) -> list[tuple[int, int, bool]]:
        """
        Read data from DNA kit file and parse in the following format:
        [
            {
                autosomal_chromosome: int,
                base_pair_position: int,
                is_homozygous: bool
            }
        ]
        :param file_name: The name of the DNA data file.
        :return: Parsed DNA autosomal SNP data formatted into tuples of chromosome, base-pair position, homozygosity.
        """
        dna_inbred_estimate_tmp: InbredEstimateAbstract | None = (
            self.__dna_inbred_estimate
        )
        if file_name.endswith(".csv"):
            if self.__is_format_identifier_in_file(
                file_name=file_name,
                unique_format_identifier=DnaKitFileUniqueIdentifier.MY_HERITAGE,
            ) and not isinstance(dna_inbred_estimate_tmp, InbredEstimateMyHeritage):
                dna_inbred_estimate_tmp = InbredEstimateMyHeritage()
            elif (
                self.__is_format_identifier_in_file(
                    file_name=file_name,
                    unique_format_identifier=DnaKitFileUniqueIdentifier.FAMILY_TREE_DNA,
                )
                or not self.__is_format_identifier_in_file(
                    file_name=file_name,
                    unique_format_identifier=DnaKitFileUniqueIdentifier.MY_HERITAGE,
                )
            ) and not isinstance(dna_inbred_estimate_tmp, InbredEstimateFamilyTreeDna):
                dna_inbred_estimate_tmp = InbredEstimateFamilyTreeDna()
            elif not isinstance(
                dna_inbred_estimate_tmp, InbredEstimateMyHeritage
            ) and not isinstance(dna_inbred_estimate_tmp, InbredEstimateFamilyTreeDna):
                raise ValueError("CSV file format is not supported.")
        elif file_name.endswith(".txt"):
            if self.__is_format_identifier_in_file(
                file_name=file_name,
                unique_format_identifier=DnaKitFileUniqueIdentifier.ANCESTRY,
            ) and not isinstance(dna_inbred_estimate_tmp, InbredEstimateAncestry):
                dna_inbred_estimate_tmp = InbredEstimateAncestry()
            elif self.__is_format_identifier_in_file(
                file_name=file_name,
                unique_format_identifier=DnaKitFileUniqueIdentifier.TWENTY_THREE_AND_ME,
            ) and not isinstance(dna_inbred_estimate_tmp, InbredEstimate23AndMe):
                dna_inbred_estimate_tmp = InbredEstimate23AndMe()
            elif not isinstance(
                dna_inbred_estimate_tmp, InbredEstimateAncestry
            ) and not isinstance(dna_inbred_estimate_tmp, InbredEstimate23AndMe):
                raise ValueError("Text file format is not supported.")
        else:
            raise ValueError("File type is not supported. Use either .csv or .txt")
        self.__dna_inbred_estimate = dna_inbred_estimate_tmp
        return (
            self.__dna_inbred_estimate.read_dna_data_file(file_name=file_name)
            if self.__dna_inbred_estimate
            else []
        )

    def compute_f_roh_estimate(
        self,
        roh_min_snp_cnt: int = 100,
        roh_min_base_pair_len: int = InbredEstimateAbstract.BASE_PAIRS_PER_MB,
        roh_max_base_pair_range: int = InbredEstimateAbstract.BASE_PAIRS_PER_MB,
        snp_max_heterozygous: int = 1,
        snp_max_avg_range_kb: int = 50,
    ) -> None:
        """
        Compute F_ROH from parsed autosomal DNA data file for estimating inbreeding/shared-ancestry coefficient score.
        :param roh_min_snp_cnt: Minimum number of SNP count required for a valid candidate ROH segment.
        :param roh_min_base_pair_len: Minimum base-pair length required for a valid candidate ROH segment.
        :param roh_max_base_pair_range: Maximum base-pair range between adjacent SNP genotypes in a valid candidate ROH.
        :param snp_max_heterozygous: Maximum number of heterozygous SNP genotypes in a valid candidate ROH segment.
        :param snp_max_avg_range_kb: Maximum avg range between adjacent SNP genotypes in a candidate ROH segment in Kb.
        """
        if self.__dna_inbred_estimate:
            self.__dna_inbred_estimate.compute_f_roh_estimate(
                roh_min_snp_cnt=roh_min_snp_cnt,
                roh_min_base_pair_len=roh_min_base_pair_len,
                roh_max_base_pair_range=roh_max_base_pair_range,
                snp_max_heterozygous=snp_max_heterozygous,
                snp_max_avg_range_kb=snp_max_avg_range_kb,
            )
