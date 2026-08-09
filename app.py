from dna.dna_estimate import dna_estimate_f_roh
from utils.parse_args import parse_args


def main() -> None:
    dna_estimate_f_roh(args=parse_args())


if "__main__" == __name__:
    main()
