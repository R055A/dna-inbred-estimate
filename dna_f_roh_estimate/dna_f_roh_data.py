from dataclasses import dataclass, asdict
from typing import Iterator


@dataclass
class DataFROH:
    """
    F_ROH analysis summary data, including F_ROH estimate representing inbreeding/shared ancestry coefficient score.
    """

    snp_autosomal_cnt_ttl: (
        int | None
    )  # Total count of usable autosomal chromosome SNPs.
    snp_homozygous_cnt: (
        int | None
    )  # Count of usable autosomal chromosome SNPs with identical alleles.
    snp_heterozygous_cnt: (
        int | None
    )  # Count of usable autosomal chromosome SNPs where alleles differ.
    snp_invalid_cnt: (
        int | None
    )  # Count of invalid autosomal chromosome SNP calls excluded from F_ROH analysis.
    roh_cnt: int | None = None  # Count of valid homozygosity runs.
    roh_len_ttl_mb: float | int | None = (
        None  # Total of all valid ROH segment lengths (Mb).
    )
    denominator_mb: float | int | None = (
        None  # Total range of autosomal chromosome data used in analysis (Mb).
    )
    f_roh: float | int | None = None  # Estimated F_ROH value
    f_roh_percent: float | int | None = None  # Estimated F_ROH value as a percentage.

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
