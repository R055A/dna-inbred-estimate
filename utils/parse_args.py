from argparse import ArgumentParser, Namespace
from dna.dna_data import DataConfigFROH


def parse_args() -> Namespace:
    default_config: DataConfigFROH = DataConfigFROH()
    parser = ArgumentParser(description="Estimate F_ROH from DNA data kit file.")
    parser.add_argument("-file", help="Raw DNA data kit file.", required=True)
    parser.add_argument(
        "-min_snp",
        help="Minimum number of SNP count required for a valid candidate ROH segment.",
        type=int,
        default=default_config.roh_min_snp_cnt,
        required=False,
    )
    parser.add_argument(
        "-min_bp_len",
        help="Minimum base-pair length required for a valid candidate ROH segment.",
        type=int,
        default=default_config.roh_min_base_pair_len,
        required=False,
    )
    parser.add_argument(
        "-max_bp_range",
        help="Maximum base-pair range between adjacent SNP genotypes in a valid candidate ROH.",
        type=int,
        default=default_config.roh_max_base_pair_range,
        required=False,
    )
    parser.add_argument(
        "-max_snp_het",
        help="Maximum number of heterozygous SNP genotypes in a valid candidate ROH segment.",
        type=int,
        default=default_config.snp_max_heterozygous,
        required=False,
    )
    parser.add_argument(
        "-max_snp_avg_range",
        help="Maximum avg range between adjacent SNP genotypes in a candidate ROH segment in Kb.",
        type=int,
        default=default_config.snp_max_avg_range_kb,
        required=False,
    )
    return parser.parse_args()
