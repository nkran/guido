from typing import List, Dict, Tuple

from guido.mmej import generate_mmej_patterns
from guido.off_targets import run_bowtie
from guido.genome import Genome
from guido.helpers import rev_comp


class Guide:
    def __init__(
        self,
        sequence: str,
        pam_position: int,
        pam_len: int,
        strand: str = "+",
        max_flanking_length: int = 75,
        cut_offset: int = 3,
        chromosome: str = "seq",
        start: int = 0,
    ) -> None:

        length = len(sequence)
        cut_pos = pam_position - cut_offset

        relative_guide_start = cut_pos - 17
        relative_guide_end = cut_pos + cut_offset + pam_len
        guide_seq = sequence[relative_guide_start:relative_guide_end]

        if strand == "-":
            cut_pos = pam_position + pam_len + cut_offset
            relative_guide_start = pam_position
            relative_guide_end = cut_pos + 17
            guide_seq = rev_comp(sequence[relative_guide_start:relative_guide_end])

        absolute_cut_pos = start + cut_pos
        absolute_guide_start = start + relative_guide_start
        absolute_guide_end = start + relative_guide_end

        left_slice = cut_pos - max_flanking_length
        right_slice = cut_pos + max_flanking_length
        relative_cut_pos = max_flanking_length

        if left_slice < 0:
            left_slice = 0
            relative_cut_pos = cut_pos
        if right_slice > length:
            right_slice = length

        # Guide properties
        self.locus_seq: str = sequence[slice(left_slice, right_slice)]
        self.guide_seq: str = guide_seq
        self.guide_pam: str = guide_seq[-pam_len:]
        self.guide_chrom: str = chromosome
        self.guide_start: int = absolute_guide_start
        self.guide_end: int = absolute_guide_end
        self.guide_strand: str = strand
        self.relative_cut_pos: int = relative_cut_pos
        self.absolute_cut_pos: int = absolute_cut_pos

        # MMEJ properties
        self.mmej_patterns: list = []
        self.mmej_sum_score: float = 0
        self.mmej_oof_score: float = 0

        # Off-targets
        self.offtargets_str: str = ""
        self.offtargets_n: int = 0
        self.offtargets_dict: dict = {}

    def __repr__(self) -> str:
        return f"Guide({self.guide_seq}|{self.guide_chrom}:{self.guide_start}-{self.guide_end}|{self.guide_strand}|)"

    @property
    def location(self) -> str:
        return f"{self.guide_chrom}:{self.guide_start}-{self.guide_end}"

    def simulate_end_joining(
        self, n_patterns: int = 5, length_weight: int = 20
    ) -> None:

        if "N" not in self.locus_seq:
            mmej_patterns = generate_mmej_patterns(
                self.relative_cut_pos, self.locus_seq, length_weight
            )

            sorted_mmej_patterns = sorted(
                mmej_patterns, key=lambda x: -x["pattern_score"]
            )[:n_patterns]

            oof_sum = sum(
                [
                    p["pattern_score"]
                    for p in sorted_mmej_patterns
                    if p["frame_shift"] == "+"
                ]
            )

            sum_score = sum([p["pattern_score"] for p in sorted_mmej_patterns])
            oof_score = oof_sum / sum_score * 100
            self.mmej_patterns = sorted_mmej_patterns
            self.mmej_sum_score = sum_score
            self.mmej_oof_score = oof_score

    def find_off_targets(
        self, genome: Genome, **kwargs
    ) -> tuple:

        index_path = genome.bowtie_index

        if index_path:
            guides_offtargets, targets, stdout = run_bowtie(
                guides=[self], genome_index_path=index_path, **kwargs
            )
            print(guides_offtargets)
            # if len(guides_offtargets.items()) > 0:
            #     for ix, g in guides_offtargets.items():
            #         self.offtargets_str = g["offtargets_str"]
            #         self.offtargets_n = g["offtargets_n"]
            #         self.offtargets_dict = g["offtargets_dict"]
            
            return guides_offtargets, targets,stdout 
        else:
            raise ValueError("Bowtie index is not built for the genome / locus.")
