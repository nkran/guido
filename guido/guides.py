from guido.mmej import generate_mmej_patterns
from guido.off_targets import run_bowtie, calculate_ot_sum_score
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
        self.rank_score: float = 0.0
        self.rank: int = 0
        self.id: str = ""

        # MMEJ properties
        self.mmej_patterns: list = []

        # Off-targets
        self.off_targets: list = []

        # Layers
        self._layers: dict = {}

    def __repr__(self) -> str:
        if self.id:
            name = self.id
        else:
            name = "gRNA"
        return f"{name}({self.guide_seq}|{self.guide_chrom}:{self.guide_start}-{self.guide_end}|{self.guide_strand}|)"

    @property
    def location(self) -> str:
        return f"{self.guide_chrom}:{self.guide_start}-{self.guide_end}"

    def simulate_end_joining(
        self, n_patterns: int = 5, length_weight: int = 20
    ) -> None:
        """
        [summary]

        Parameters
        ----------
        n_patterns : int, optional
            [description], by default 5
        length_weight : int, optional
            [description], by default 20
        """

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

            self.add_layer("mmej_sum_score", sum_score)
            self.add_layer("mmej_oof_score", oof_score)

    def find_off_targets(self, genome: Genome, **kwargs) -> None:
        """
        [summary]

        Parameters
        ----------
        genome : Genome
            [description]

        Raises
        ------
        ValueError
            [description]
        ValueError
            [description]
        """

        if genome:
            index_path = genome.bowtie_index
        else:
            raise ValueError("No genome / locus specified.")

        if index_path:
            guides_offtargets = run_bowtie(
                guides=[self], genome_index_path=index_path, **kwargs
            )

            self.off_targets = guides_offtargets[0]
            self.add_layer("ot_sum_score", calculate_ot_sum_score(self.off_targets))

        else:
            raise ValueError("Bowtie index is not built for the genome / locus.")

    @property
    def layers(self):
        return self._layers

    def add_layer(self, name: str, layer_data: float):
        self._layers[name] = layer_data

    def layer(self, key):
        if key in self._layers.keys():
            return self._layers[key]
        else:
            raise ValueError(f"Layer {key} was not added to the gRNA.")
