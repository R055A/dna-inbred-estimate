from dna_f_roh_estimate.dna_f_roh_enum import (
    DnaKitFileHeaderAbstract,
    DnaKitFileInvalidAllelesSymbol,
)
from dna_f_roh_estimate.dna_f_roh_data import DataFROH
from utils import read_data_file_all_lines
from abc import abstractmethod, ABC


class InbredEstimateAbstract(ABC):
    """
    Abstract class for estimating F_ROH, representing inbreeding/shared ancestry coefficient score, from DNA kit data.
    """

    BASE_PAIRS_PER_MB: int = 1_000_000

    def __init__(
        self,
        file_name: str | None = None,
        headers: list[str] | None = None,
        invalid_alleles: set[DnaKitFileInvalidAllelesSymbol] | None = None,
    ) -> None:
        """
        Class constructor for instantiating an F_ROH estimator instance with optional file name and genotype headers.
        :param file_name: Path name to a supported DNA kit data file. If provided, the file is read at instantiation.
        :param headers: List of format-specific genotype header names in addition to [RSID, CHROMOSOME, POSITION].
        :param invalid_alleles: Allele file entries which are not valid for F_ROH analysis.
        """
        self.__invalid_alleles: set[DnaKitFileInvalidAllelesSymbol] = (
            invalid_alleles if invalid_alleles else set()
        )
        self.__snp_data: list[tuple[int, int, bool]] = []
        self.__f_roh_data: DataFROH | None = None
        self._headers: list[str] = [h.value for h in DnaKitFileHeaderAbstract] + (
            headers or []
        )
        if file_name:
            self.read_dna_data_file(file_name=file_name)

    @property
    def f_roh_data(self) -> DataFROH | None:
        """
        :return: Data from F_ROH analysis.
        """
        return self.__f_roh_data

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
        return self.__snp_data

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
        self.__snp_data = snp_data

    @abstractmethod
    def _split_data_row(self, data_row: str) -> list[str]:
        """
        Split a data row into column values relevant to a subclass file format.
        :param data_row: Row of data from the DNA kit file.
        :return: List of column values split from the data row.
        """
        pass

    @abstractmethod
    def _iter_data_rows(self, header: dict[str, int], data_rows: list[str]):
        """
        Yield valid parsed data rows relevant to a subclass file format.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_rows: List of data rows from the DNA kit file.
        :return: Iterator of parsed data rows for F_ROH analysis.
        """
        pass

    @abstractmethod
    def _get_alleles_data(
        self, header: dict[str, int], data_row: list[str]
    ) -> list[str | None]:
        """
        Get allele values from a parsed data row relevant to a subclass file format.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: List containing allele pair, or None values if the genotype is invalid.
        """
        pass

    def __get_autosomal_chromosome_num(
        self, header: dict[str, int], data_row: list[str]
    ) -> int | None:
        """
        Get and validate an autosomal chromosome number from a parsed data row. Autosomal is between 1 and 22 (not sex).
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: Autosomal chromosome number between 1 and 22, or None if invalid.
        """
        chromosome = (
            data_row[header[DnaKitFileHeaderAbstract.HEADER_CHROMOSOME.value]]
            .strip()
            .upper()
        )
        if chromosome.isdigit():
            chromosome = int(chromosome)
            if 1 <= chromosome <= 22:
                return chromosome
        return None

    def __get_base_pair_position_num(
        self, header: dict[str, int], data_row: list[str]
    ) -> int | None:
        """
        Get and validate the base-pair position number from a parsed data row.
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_row: List of parsed data rows from the DNA kit file.
        :return: Base-pair position value, or None if invalid.
        """
        position = data_row[
            header[DnaKitFileHeaderAbstract.HEADER_POSITION.value]
        ].strip()
        if position.isdigit():
            return int(position)
        return None

    def _is_valid_alleles(self, alleles_data: str) -> bool:
        """
        Conditional for if an allele value is valid for F_ROH analysis.
        :param alleles_data: Allele value from the DNA kit file.
        :return: True if the value is valid; otherwise False.
        """
        return alleles_data not in self.__invalid_alleles

    def __init_data(
        self,
        snp_genotypes: list[tuple[int, int, bool]],
        snp_invalid_cnt: int,
    ) -> None:
        """
        Store initial SNP data from the parsed DNA kit file for F_ROH analysis.
        :param snp_genotypes: Parsed autosomal SNP genotypes: {autosomal_chromosome, base_pair_position, is_homozygous}.
        :param snp_invalid_cnt: Count of invalid autosomal SNP genotypes excluded from F_ROH analysis.
        """
        snp_autosomal_cnt_ttl = len(snp_genotypes)
        snp_homozygous_cnt = sum(
            1 for _, _, is_homozygous in snp_genotypes if is_homozygous
        )
        snp_heterozygous_cnt = snp_autosomal_cnt_ttl - snp_homozygous_cnt
        self.__f_roh_data = DataFROH(
            snp_autosomal_cnt_ttl=snp_autosomal_cnt_ttl,
            snp_homozygous_cnt=snp_homozygous_cnt,
            snp_heterozygous_cnt=snp_heterozygous_cnt,
            snp_invalid_cnt=snp_invalid_cnt,
        )

    def __parse_dna_data(
        self, header: dict[str, int], data_rows: list[str]
    ) -> list[tuple[int, int, bool]]:
        """
        Parse DNA kit data for F_ROH analysis of autosomal SNP genotypes in the following format:
        [
            {
                autosomal_chromosome: int,
                base_pair_position: int,
                is_homozygous: bool
            }
        ]
        :param header: Normalized genotype header names mapped to data column indexes.
        :param data_rows: List of raw data rows from the DNA kit file.
        :return: Parsed DNA autosomal SNP data formatted into tuples of chromosome, base-pair position, homozygosity.
        """
        snp_genotypes: list[tuple[int, int, bool]] = []
        snp_invalid_cnt: int = 0
        for data_row in self._iter_data_rows(header=header, data_rows=data_rows):
            autosomal_chromosome: int | None = self.__get_autosomal_chromosome_num(
                header=header, data_row=data_row
            )
            if not autosomal_chromosome:
                continue
            position: int | None = self.__get_base_pair_position_num(
                header=header, data_row=data_row
            )
            if not position:
                continue
            alleles_one, alleles_two = self._get_alleles_data(
                header=header, data_row=data_row
            )[:2]
            if not alleles_one or not alleles_two:
                snp_invalid_cnt += 1
                continue
            snp_genotypes.append(
                (autosomal_chromosome, position, alleles_one == alleles_two)
            )
        snp_genotypes.sort(key=lambda snp: (snp[0], snp[1]))
        self.__init_data(
            snp_genotypes=snp_genotypes,
            snp_invalid_cnt=snp_invalid_cnt,
        )
        return snp_genotypes

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
        :param file_name: Path name to a supported DNA kit data file.
        :return: Parsed DNA autosomal SNP data formatted into tuples of chromosome, base-pair position, homozygosity.
        """
        file_lines: list[str] = read_data_file_all_lines(file_name=file_name)
        for i, line in enumerate(file_lines):
            cols = [
                c.strip("#").strip().upper()
                for c in self._split_data_row(data_row=line)
            ]
            if all(col.upper() in cols for col in self._headers):
                self.__snp_data = self.__parse_dna_data(
                    header={
                        col.upper(): cols.index(col.upper()) for col in self._headers
                    },
                    data_rows=file_lines[i + 1 :],
                )
                return self.__snp_data
        raise ValueError("Data format is not supported.")

    def __is_valid_roh(
        self,
        length_base_pair: int,
        snp_count: int,
        roh_min_base_pair_len: int,
        roh_min_snp_cnt: int,
        snp_max_avg_range_kb: int,
    ) -> bool:
        """
        Conditional for if a candidate homozygous satisfies valid ROH thresholds.
        :param length_base_pair: Base-pair length of candidate ROH segment.
        :param snp_count: Count of SNP genotypes in candidate ROH segment.
        :param roh_min_base_pair_len: Minimum valid base-pair ROH segment length.
        :param roh_min_snp_cnt: Minimum required SNP count in valid ROH segment.
        :param snp_max_avg_range_kb: Maximum valid average range between SNP genotypes in Kb.
        :return: True if candidate segment is valid ROH; otherwise False.
        """
        return (
            length_base_pair >= roh_min_base_pair_len
            and snp_count >= roh_min_snp_cnt
            and (length_base_pair / 1000) / snp_count <= snp_max_avg_range_kb
        )

    def __is_valid_roh_base_pair_range(
        self,
        base_pair_position: int,
        prev_position: int | None,
        roh_max_base_pair_range: int,
    ) -> bool:
        """
        Conditional for if the base-pair range between adjacent SNP genotypes in ROH is valid.
        :param base_pair_position: Base-pair position of the current SNP genotype.
        :param prev_position: Base-pair position of the previous SNP genotype, or None if first SNP in candidate ROH.
        :param roh_max_base_pair_range: Maximum valid base-pair range between adjacent SNP genotypes in candidate ROH.
        :return: True if the base-pair range between adjacent SNP genotypes in ROH is valid; otherwise False.
        """
        return not (
            prev_position is not None
            and base_pair_position - prev_position > roh_max_base_pair_range
        )

    def __compute_autosomal_chromosome_range(
        self,
        snp_pos_first: dict[int, int],
        snp_pos_last: dict[int, int],
    ) -> int:
        """
        Compute base-pair autosomal chromosome range.
        :param snp_pos_first: First valid SNP base-pair position per autosomal chromosome.
        :param snp_pos_last: Last valid SNP base-pair position per autosomal chromosome.
        :return: Total base-pair autosomal chromosome range value.
        """
        return sum(
            snp_pos_last[autosomal_chromosome] - snp_pos_first[autosomal_chromosome] + 1
            for autosomal_chromosome in snp_pos_first
        )

    def __set_roh_data(
        self,
        roh_cnt: int,
        roh_base_pairs_ttl: int,
        autosomal_chromosome_range: int,
    ) -> None:
        """
        Store computed ROH analysis summary data with F_ROH estimate score and percentage.
        :param roh_cnt: Count of valid ROH segments.
        :param roh_base_pairs_ttl: Total length of valid ROH segments in base-pairs.
        :param autosomal_chromosome_range: Total valid autosomal chromosome range in base-pairs.
        """
        self.__f_roh_data.roh_cnt = roh_cnt
        self.__f_roh_data.roh_len_ttl_mb = roh_base_pairs_ttl / self.BASE_PAIRS_PER_MB
        self.__f_roh_data.denominator_mb = (
            autosomal_chromosome_range / self.BASE_PAIRS_PER_MB
        )
        self.__f_roh_data.f_roh = roh_base_pairs_ttl / autosomal_chromosome_range
        self.__f_roh_data.f_roh_percent = (
            roh_base_pairs_ttl / autosomal_chromosome_range
        ) * 100

    def compute_f_roh_estimate(
        self,
        roh_min_snp_cnt: int = 100,
        roh_min_base_pair_len: int = BASE_PAIRS_PER_MB,
        roh_max_base_pair_range: int = BASE_PAIRS_PER_MB,
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
        if not self.__snp_data:
            raise AssertionError("Missing parsed DNA data.")
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
                if self.__is_valid_roh(
                    length_base_pair=length_base_pair,
                    snp_count=candidate_count,
                    roh_min_base_pair_len=roh_min_base_pair_len,
                    roh_min_snp_cnt=roh_min_snp_cnt,
                    snp_max_avg_range_kb=snp_max_avg_range_kb,
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
        for (
            autosomal_chromosome,
            base_pair_position,
            is_homozygous,
        ) in self.__snp_data:
            if autosomal_chromosome not in snp_pos_first:
                snp_pos_first[autosomal_chromosome] = base_pair_position
            snp_pos_last[autosomal_chromosome] = base_pair_position

            if autosomal_chromosome != cur_chromosome:
                cur_chromosome = autosomal_chromosome
                prev_position = None
                validate_candidate_roh()

            if not self.__is_valid_roh_base_pair_range(
                base_pair_position=base_pair_position,
                prev_position=prev_position,
                roh_max_base_pair_range=roh_max_base_pair_range,
            ):
                validate_candidate_roh()
            prev_position = base_pair_position

            if not is_homozygous and (
                candidate_start is None
                or candidate_heterozygous_count >= snp_max_heterozygous
            ):
                validate_candidate_roh()
                continue
            elif is_homozygous and candidate_start is None:
                candidate_start = base_pair_position
            elif (
                not is_homozygous
                and candidate_start is not None
                and candidate_heterozygous_count < snp_max_heterozygous
            ):
                candidate_heterozygous_count += 1
            candidate_end = base_pair_position
            candidate_count += 1
        validate_candidate_roh()

        autosomal_chromosome_range = self.__compute_autosomal_chromosome_range(
            snp_pos_first=snp_pos_first,
            snp_pos_last=snp_pos_last,
        )
        if autosomal_chromosome_range == 0:
            raise ValueError("Could not estimate F_ROH from DNA data.")

        self.__set_roh_data(
            roh_cnt=roh_cnt,
            roh_base_pairs_ttl=roh_base_pairs_ttl,
            autosomal_chromosome_range=autosomal_chromosome_range,
        )
