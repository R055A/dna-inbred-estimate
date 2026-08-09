from dataclasses import dataclass, asdict
from typing import Iterator

BASE_PAIRS_PER_MB: int = 1_000_000


@dataclass(frozen=True, slots=True)
class DataConfigFROH:
    """
    Hyperparameters used in computing F_ROH.
    :ivar roh_min_snp_cnt: Minimum number of SNP count required for a valid candidate ROH segment.
    :ivar roh_min_base_pair_len: Minimum base-pair length required for a valid candidate ROH segment.
    :ivar roh_max_base_pair_range: Maximum base-pair range between adjacent SNP genotypes in a valid candidate ROH.
    :ivar snp_max_heterozygous: Maximum number of heterozygous SNP genotypes in a valid candidate ROH segment.
    :ivar snp_max_avg_range_kb: Maximum avg range between adjacent SNP genotypes in a candidate ROH segment in Kb.
    """

    roh_min_snp_cnt: int = 100
    roh_min_base_pair_len: int = BASE_PAIRS_PER_MB
    roh_max_base_pair_range: int = BASE_PAIRS_PER_MB
    snp_max_heterozygous: int = 1
    snp_max_avg_range_kb: int = 50


@dataclass(frozen=True, slots=True)
class DataGenotype:
    """
    A single autosomal SNP genotype.
    :ivar rs_id: SNP identifier, typically an rsID.
    :ivar chromosome: Autosomal chromosome number.
    :ivar position: SNP base-pair position on the chromosome.
    :ivar allele_one: First observed allele.
    :ivar allele_two: Second observed allele.
    """

    rs_id: str
    chromosome: int
    position: int
    allele_one: str
    allele_two: str

    @property
    def is_homozygous(self) -> bool:
        return self.allele_one == self.allele_two


@dataclass(frozen=True, slots=True)
class DataDNA:
    """
    Parsed autosomal SNP data used for F_ROH analysis.
    :ivar genotypes: Valid autosomal SNP genotypes included in the analysis.
    :ivar invalid_cnt: Count of invalid autosomal SNP calls excluded from the analysis.
    """

    genotypes: tuple[DataGenotype, ...]
    invalid_cnt: int


@dataclass
class DataFROH:
    """
    F_ROH analysis summary data, including F_ROH estimate representing inbreeding/shared ancestry coefficient score.
    :ivar snp_autosomal_cnt_ttl: Total count of usable autosomal chromosome SNPs.
    :ivar snp_homozygous_cnt: Count of usable autosomal chromosome SNPs with identical alleles.
    :ivar snp_heterozygous_cnt: Count of usable autosomal chromosome SNPs where alleles differ.
    :ivar snp_invalid_cnt: Count of invalid autosomal chromosome SNP calls excluded from F_ROH analysis.
    :ivar roh_cnt: Count of valid homozygosity runs.
    :ivar roh_len_ttl_mb: Total of all valid ROH segment lengths (Mb).
    :ivar denominator_mb: Total range of autosomal chromosome data used in analysis (Mb).
    :ivar f_roh: Estimated F_ROH value.
    :ivar f_roh_percent: Estimated F_ROH value as a percentage.
    """

    snp_autosomal_cnt_ttl: int | None
    snp_homozygous_cnt: int | None
    snp_heterozygous_cnt: int | None
    snp_invalid_cnt: int | None
    roh_cnt: int | None = None
    roh_len_ttl_mb: float | int | None = None
    denominator_mb: float | int | None = None
    f_roh: float | int | None = None
    f_roh_percent: float | int | None = None

    def __repr__(self) -> str:
        """
        :return: F_ROH analysis summary data in string format.
        """
        return (
            f"Total Autosomal SNPs:\t{self.snp_autosomal_cnt_ttl if self.snp_autosomal_cnt_ttl is not None else 'N/A'}\n"
            f"Homozygous SNPs:\t{self.snp_homozygous_cnt if self.snp_homozygous_cnt is not None else 'N/A'}\n"
            f"Heterozygous SNPs:\t{self.snp_heterozygous_cnt if self.snp_heterozygous_cnt is not None else 'N/A'}\n"
            f"Invalid SNPs:\t\t{self.snp_invalid_cnt if self.snp_invalid_cnt is not None else 'N/A'}\n"
            f"ROH count:\t\t{self.roh_cnt if self.roh_cnt is not None else 'N/A'}\n"
            f"Total ROH length:\t{f"{self.roh_len_ttl_mb:.2f}" if self.roh_len_ttl_mb is not None else 'N/A'} Mb\n"
            f"Denominator:\t\t{f"{self.denominator_mb:.2f}" if self.denominator_mb is not None else 'N/A'} Mb\n"
            f"F_ROH score:\t\t{f"{self.f_roh:.6f}" if self.f_roh is not None else 'N/A'}\n"
            f"F_ROH percent:\t\t{f"{self.f_roh_percent:.2f}" if self.f_roh_percent is not None else 'N/A'}%"
        )

    def __iter__(self) -> Iterator[tuple[str, int | float | None]]:
        """
        :return: F_ROH analysis summary data in dictionary format.
        """
        yield from asdict(self).items()
