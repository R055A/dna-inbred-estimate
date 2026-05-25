from dna_f_roh_estimate.dna_f_roh_estimate import InbredEstimate
from utils import parse_args


def main():
    args = parse_args()
    dna_inbred_estimate: InbredEstimate = InbredEstimate(file_name=args.file)
    dna_inbred_estimate.compute_f_roh_estimate(*list(vars(args).copy().values())[1:])
    print(dna_inbred_estimate)


if "__main__" == __name__:
    main()
