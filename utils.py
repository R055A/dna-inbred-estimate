from argparse import ArgumentParser, Namespace


def parse_args() -> Namespace:
    parser = ArgumentParser(description="Estimate F_ROH from DNA data kit file.")
    parser.add_argument("-file", help="Raw DNA data kit file.", required=True)
    parser.add_argument(
        "-min_snp",
        help="Minimum number of SNP count required for a valid candidate ROH segment.",
        type=int,
        default=100,
        required=False
    )
    parser.add_argument(
        "-min_bp_len",
        help="Minimum base-pair length required for a valid candidate ROH segment.",
        type=int,
        default=1_000_000,
        required=False
    )
    parser.add_argument(
        "-max_bp_range",
        help="Maximum base-pair range between adjacent SNP genotypes in a valid candidate ROH.",
        type=int,
        default=1_000_000,
        required=False
    )
    parser.add_argument(
        "-max_snp_het",
        help="Maximum number of heterozygous SNP genotypes in a valid candidate ROH segment.",
        type=int,
        default=1,
        required=False
    )
    parser.add_argument(
        "-max_snp_avg_range",
        help="Maximum avg range between adjacent SNP genotypes in a candidate ROH segment in Kb.",
        type=int,
        default=50,
        required=False
    )
    return parser.parse_args()


def open_data_file(file_name: str) -> list[str]:
    with open(file_name, "r", encoding="utf-8-sig") as f:
        return f.read().splitlines()
