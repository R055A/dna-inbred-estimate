from dna.dna_data import DataDNA, DataConfigFROH, BASE_PAIRS_PER_MB, DataFROH
from abc import ABC, abstractmethod


class EstimateAbstractFROH(ABC):
    """
    Abstract class for computing F_ROH to DataFROH dataclass from DataDNA-formatted data.
    """

    def __init__(self, dna_data: DataDNA) -> None:
        """
        Class constructor with optional file name and genotype headers.
        :param dna_data: DNA autosomal SNP data in DataDNA dataclass format.
        """
        self._dna_data: DataDNA = dna_data
        self.__f_roh_data: DataFROH = self.__init_f_roh_data()

    def __init_f_roh_data(self) -> DataFROH:
        """
        Store initial SNP data from the DataDNA-formatted data for F_ROH analysis.
        :return: DataFROH dataclass instance with SNP data from DataDNA.
        """
        snp_autosomal_cnt_ttl: int = len(self._dna_data.genotypes)
        snp_homozygous_cnt: int = sum(
            genotype.is_homozygous for genotype in self._dna_data.genotypes
        )
        return DataFROH(
            snp_autosomal_cnt_ttl=snp_autosomal_cnt_ttl,
            snp_homozygous_cnt=snp_homozygous_cnt,
            snp_heterozygous_cnt=snp_autosomal_cnt_ttl - snp_homozygous_cnt,
            snp_invalid_cnt=self._dna_data.invalid_cnt,
        )

    @property
    def data(self) -> DataFROH | None:
        """
        :return: Data from F_ROH analysis in DataFROH dataclass format.
        """
        return self.__f_roh_data

    @staticmethod
    def _compute_chromosome_range(
        snp_pos_first: dict[int, int],
        snp_pos_last: dict[int, int],
    ) -> int:
        """
        Compute base-pair autosomal chromosome range.
        :param snp_pos_first: First valid SNP base-pair position per autosomal chromosome.
        :param snp_pos_last: Last valid SNP base-pair position per autosomal chromosome.
        :return: Total base-pair autosomal chromosome range value.
        """
        chromosome_range: int = sum(
            snp_pos_last[chromosome] - snp_pos_first[chromosome] + 1
            for chromosome in snp_pos_first
        )
        if chromosome_range == 0:
            raise ValueError("Could not estimate F_ROH from DNA data.")
        return chromosome_range

    def _set_roh_data(
        self,
        roh_cnt: int,
        roh_base_pairs_ttl: int,
        chromosome_range: int,
    ) -> None:
        """
        Store computed ROH analysis summary data with F_ROH estimate score and percentage.
        :param roh_cnt: Count of valid ROH segments.
        :param roh_base_pairs_ttl: Total length of valid ROH segments in base-pairs.
        :param chromosome_range: Total valid autosomal chromosome range in base-pairs.
        """
        self.__f_roh_data.roh_cnt = roh_cnt
        self.__f_roh_data.roh_len_ttl_mb = roh_base_pairs_ttl / BASE_PAIRS_PER_MB
        self.__f_roh_data.denominator_mb = chromosome_range / BASE_PAIRS_PER_MB
        self.__f_roh_data.f_roh = roh_base_pairs_ttl / chromosome_range
        self.__f_roh_data.f_roh_percent = roh_base_pairs_ttl / chromosome_range * 100

    def _assert_snp_data(self) -> None:
        if not self._dna_data:
            raise AssertionError("Missing parsed DNA data.")

    @abstractmethod
    def compute_f_roh(
        self,
        config: DataConfigFROH,
    ) -> None:
        """
        Compute F_ROH on autosomal SNP data in DataDNA dataclass format.
        :param config: hyperparameters for computing F_ROH in DataConfigFROH dataclass format.
        """
        pass
