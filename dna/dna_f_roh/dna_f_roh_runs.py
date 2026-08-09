from dna.dna_f_roh.dna_f_roh_abstract import EstimateAbstractFROH
from dna.dna_data import DataDNA, DataConfigFROH


class EstimateRunsFROH(EstimateAbstractFROH):
    """
    Compute F_ROH to DataFROH dataclass from DataDNA-formatted data using consecutive runs. Extends EstimateAbstractFROH
    """

    def __init__(self, dna_data: DataDNA) -> None:
        """
        Class constructor with optional file name and genotype headers.
        :param dna_data: DNA autosomal SNP data in DataDNA dataclass format.
        """
        super().__init__(dna_data=dna_data)

    def compute_f_roh(self, config: DataConfigFROH) -> None:
        """
        Compute F_ROH on autosomal SNP data in DataDNA dataclass format.
        :param config: hyperparameters for computing F_ROH in DataConfigFROH dataclass format.
        """
        self._assert_snp_data()
        (
            candidate_start,
            candidate_end,
            candidate_count,
            candidate_heterozygous_count,
        ) = (None, None, 0, 0)
        roh_cnt, roh_base_pairs_ttl = 0, 0

        def validate_candidate_roh():
            """
            Add candidate ROH segment to base pairs total if valid ROH length, SNP count and density threshold.
            """
            nonlocal candidate_start, candidate_end, candidate_count, candidate_heterozygous_count
            nonlocal roh_cnt, roh_base_pairs_ttl

            if candidate_start is not None and candidate_end is not None:
                length_base_pair = candidate_end - candidate_start + 1
                if (
                    length_base_pair >= config.roh_min_base_pair_len
                    and candidate_count >= config.roh_min_snp_cnt
                    and (length_base_pair / 1_000) / candidate_count
                    <= config.snp_max_avg_range_kb
                ):
                    roh_base_pairs_ttl += length_base_pair
                    roh_cnt += 1
                (
                    candidate_start,
                    candidate_end,
                    candidate_count,
                    candidate_heterozygous_count,
                ) = (None, None, 0, 0)

        snp_pos_first, snp_pos_last = {}, {}
        cur_chromosome, prev_position = None, None
        for genotype in self._dna_data.genotypes:
            if genotype.chromosome not in snp_pos_first:
                snp_pos_first[genotype.chromosome] = genotype.position
            snp_pos_last[genotype.chromosome] = genotype.position

            if genotype.chromosome != cur_chromosome:
                cur_chromosome = genotype.chromosome
                prev_position = None
                validate_candidate_roh()

            if (
                prev_position is not None
                and genotype.position - prev_position > config.roh_max_base_pair_range
            ):
                validate_candidate_roh()
            prev_position = genotype.position

            if not genotype.is_homozygous and (
                candidate_start is None
                or candidate_heterozygous_count >= config.snp_max_heterozygous
            ):
                validate_candidate_roh()
                continue
            elif genotype.is_homozygous and candidate_start is None:
                candidate_start = genotype.position
            elif (
                not genotype.is_homozygous
                and candidate_start is not None
                and candidate_heterozygous_count < config.snp_max_heterozygous
            ):
                candidate_heterozygous_count += 1
            candidate_end = genotype.position
            candidate_count += 1
        validate_candidate_roh()

        chromosome_range = self._compute_chromosome_range(
            snp_pos_first=snp_pos_first,
            snp_pos_last=snp_pos_last,
        )
        self._set_roh_data(
            roh_cnt=roh_cnt,
            roh_base_pairs_ttl=roh_base_pairs_ttl,
            chromosome_range=chromosome_range,
        )
