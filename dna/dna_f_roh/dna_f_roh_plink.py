from dna.dna_f_roh.dna_f_roh_abstract import EstimateAbstractFROH
from dna.dna_data import DataDNA, DataGenotype, DataConfigFROH
from utils.install_plink import install_plink
from tempfile import TemporaryDirectory
from dna.dna_enum import PlinkHeader
from subprocess import run
from pathlib import Path


class EstimatePlinkFROH(EstimateAbstractFROH):
    """
    Compute F_ROH to DataFROH dataclass from DataDNA-formatted data using PLINK. Extends EstimateAbstractFROH.
    """

    def __init__(self, dna_data: DataDNA) -> None:
        """
        Class constructor with optional file name and genotype headers.
        :param dna_data: DNA autosomal SNP data in DataDNA dataclass format.
        """
        super().__init__(dna_data=dna_data)
        self.__plink_exec_path: Path | None = install_plink()
        self.__roh_cnt = 0
        self.__roh_base_pairs_ttl = 0

    def __write_plink_input_files(self, in_prefix: Path) -> None:
        used_ids = set()
        map_rows: list[str] = []
        ped_alleles: list = []
        for i, genotype in enumerate[DataGenotype](
            self._dna_data.genotypes,
            start=1,
        ):
            rs_id: str = genotype.rs_id
            chromosome: int = genotype.chromosome
            pos: int = genotype.position
            if not rs_id or any(char.isspace() for char in rs_id) or rs_id in used_ids:
                rs_id = f"chr{chromosome}:{pos}:{i}"
            used_ids.add(rs_id)
            map_rows.append(f"{chromosome}\t{rs_id}\t0\t{pos}")
            ped_alleles.extend((genotype.allele_one, genotype.allele_two))

        Path(f"{in_prefix}.map").write_text(
            data="\n".join(map_rows) + "\n", encoding="utf-8"
        )
        Path(f"{in_prefix}.ped").write_text(
            data=" ".join(["DNA", "SAMPLE", "0", "0", "0", "-9"] + ped_alleles) + "\n",
            encoding="utf-8",
        )

    def __run_plink(
        self, in_prefix: Path, out_prefix: Path, config: DataConfigFROH
    ) -> None:
        plink_res = run(
            args=[
                self.__plink_exec_path,
                "--file",
                str(in_prefix),
                "--homozyg",
                "--homozyg-snp",
                str(config.roh_min_snp_cnt),
                "--homozyg-kb",
                f"{config.roh_min_base_pair_len / 1000:g}",
                "--homozyg-density",
                str(config.snp_max_avg_range_kb),
                "--homozyg-gap",
                f"{config.roh_max_base_pair_range / 1000:g}",
                "--homozyg-het",
                str(config.snp_max_heterozygous),
                "--homozyg-window-het",
                str(config.snp_max_heterozygous),
                "--allow-no-sex",
                "--out",
                str(out_prefix),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if plink_res.returncode != 0:
            raise RuntimeError(
                f"PLINK fail {plink_res.returncode}: {plink_res.stderr or plink_res.stdout}"
            )

    def __read_plink_output_file(self, out_prefix: Path) -> None:
        self.__roh_cnt = 0
        self.__roh_base_pairs_ttl = 0

        hom_file = Path(f"{out_prefix}.hom")
        if hom_file.exists() and hom_file.stat().st_size != 0:
            lines = hom_file.read_text(encoding="utf-8").splitlines()
            if lines:
                header = lines[0].split()
                if not all(v in header for v in [h.value for h in PlinkHeader]):
                    raise RuntimeError("Invalid PLINK .hom format.")

                for line in lines[1:]:
                    if not line.strip():
                        continue
                    line = line.strip().split()
                    chromosome = int(line[header.index(PlinkHeader.CHR.value)])
                    start = int(line[header.index(PlinkHeader.POS_ONE.value)])
                    end = int(line[header.index(PlinkHeader.POS_TWO.value)])
                    if 1 <= chromosome <= 22 and end >= start:
                        self.__roh_cnt += 1
                        self.__roh_base_pairs_ttl += end - start + 1

    def compute_f_roh(self, config: DataConfigFROH) -> None:
        """
        Compute F_ROH using PLINK on autosomal SNP data in DataDNA dataclass format.
        :param config: hyperparameters for computing F_ROH in DataConfigFROH dataclass format.
        """
        self._assert_snp_data()
        if not self.__plink_exec_path:
            raise AssertionError("Missing PLINK executable.")

        with TemporaryDirectory[str](prefix="plink-") as tmp_dir:
            in_prefix: Path = Path(tmp_dir) / "genotypes"
            out_prefix: Path = Path(tmp_dir) / "roh"

            self.__write_plink_input_files(in_prefix=in_prefix)
            self.__run_plink(in_prefix=in_prefix, out_prefix=out_prefix, config=config)
            self.__read_plink_output_file(out_prefix=out_prefix)

        snp_pos_first = {}
        snp_pos_last = {}
        for genotype in self._dna_data.genotypes:
            if genotype.chromosome not in snp_pos_first:
                snp_pos_first[genotype.chromosome] = genotype.position
            snp_pos_last[genotype.chromosome] = genotype.position

        chromosome_range = self._compute_chromosome_range(
            snp_pos_first=snp_pos_first,
            snp_pos_last=snp_pos_last,
        )
        self._set_roh_data(
            roh_cnt=self.__roh_cnt,
            roh_base_pairs_ttl=self.__roh_base_pairs_ttl,
            chromosome_range=chromosome_range,
        )
