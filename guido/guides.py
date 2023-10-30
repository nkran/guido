import os
import warnings
from collections import Counter

import numpy as np
from azimuth import model_comparison

from .helpers import load_cfd_scoring_matrix, rev_comp
from .mmej import generate_mmej_patterns
from .off_targets import calculate_cfd_score, calculate_ot_sum_score, run_bowtie

# suppress warnings from azimuth
warnings.filterwarnings("ignore")


class Guide:
    def __init__(
        self,
        sequence,
        pam_position,
        pam_len,
        strand="+",
        max_flanking_length=75,
        cut_offset=3,
        chromosome="seq",
        start=0,
    ):
        """Guide class.

        Defines a gRNA and its properties.

        Parameters
        ----------
        sequence : str
            gRNA sequence
        pam_position : int
            Position of PAM sequence in gRNA sequence
        pam_len : int
            Length of PAM sequence
        strand : str, optional
            Strand of the gRNA, by default "+"
        max_flanking_length : int, optional
            Max flanking length taken into account for MMEJ calculation, by default 75
        cut_offset : int, optional
            Cut offset from the beginning of PAM, by default 3
        chromosome : str, optional
            Chromosome where the gRNA is located, by default "seq"
        start : int, optional
            Start position of the gRNA on the chromosome, by default 0
        """

        length = len(sequence)
        cut_pos = pam_position - cut_offset

        relative_guide_start = cut_pos - 17
        relative_guide_end = cut_pos + cut_offset + pam_len
        guide_seq = sequence[relative_guide_start:relative_guide_end]
        guide_long_seq = sequence[relative_guide_start - 4 : relative_guide_end + 3]

        if strand == "-":
            cut_pos = pam_position + pam_len + cut_offset
            relative_guide_start = pam_position
            relative_guide_end = cut_pos + 17
            guide_seq = rev_comp(sequence[relative_guide_start:relative_guide_end])
            guide_long_seq = rev_comp(
                sequence[relative_guide_start - 3 : relative_guide_end + 4]
            )

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
        self.locus_seq = sequence[slice(left_slice, right_slice)]
        self.guide_seq = guide_seq
        self.guide_long_seq = guide_long_seq
        self.guide_pam = guide_seq[-pam_len:]
        self.guide_chrom = chromosome
        self.guide_start = absolute_guide_start
        self.guide_end = absolute_guide_end - 1  # 1-based (last position inclusive)
        self.guide_strand = strand
        self.relative_cut_pos = relative_cut_pos
        self.absolute_cut_pos = absolute_cut_pos
        self.rank_score = 0.0
        self.rank = 0
        self.id = ""

        # MMEJ properties
        self.mmej_patterns = []

        # Off-targets
        self.off_targets = []
        self.off_target_str = ""

        # Layers
        self._layers = {}

    def __repr__(self):
        if self.id:
            name = self.id
        else:
            name = "gRNA"
        return f"{name}({self.guide_seq}|{self.guide_chrom}:{self.guide_start}-{self.guide_end}|{self.guide_strand}|)"

    def __getattr__(self, attr):
        return self[attr]

    @property
    def location(self):
        """Returns the location of the guide in the format: chr:start-end.

        Returns
        -------
        str
            String representation of gRNA location
        """
        return f"{self.guide_chrom}:{self.guide_start}-{self.guide_end}"

    def _create_mmej_oof_string(self, mmej_patterns):
        return "|".join([p["frame_shift"] for p in mmej_patterns])

    def simulate_end_joining(self, n_patterns=5, length_weight=20):
        """Simulate Microhomology-Mediated End Joining (MMEJ) events for the
        gRNA.

        MMEJ scoring is based on the Bae et al. 2014 paper (https://doi.org/10.1038/nmeth.3015)

        Parameters
        ----------
        n_patterns : int, optional
            Number of top-scoring MMEJ patterns to keep, by default 5
        length_weight : int, optional
            Lengeth weight, by default 20
        """

        if "N" not in self.locus_seq:
            mmej_patterns = generate_mmej_patterns(
                self.relative_cut_pos, self.locus_seq, length_weight
            )

            if mmej_patterns:
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
            else:
                sorted_mmej_patterns = None
                sum_score = None
                oof_score = None

            # TODO: refactor
            self.mmej_patterns = sorted_mmej_patterns
            self.mmej_str = self._create_mmej_oof_string(sorted_mmej_patterns)

            self.add_layer("mmej_sum_score", sum_score)
            self.add_layer("mmej_oof_score", oof_score)
        else:
            self.mmej_patterns = None
            self.mmej_str = None

            self.add_layer("mmej_sum_score", 0)
            self.add_layer("mmej_oof_score", 0)

    def find_off_targets(self, genome, **kwargs):
        """Finds off-targets for the guide. The off-targets are found using
        Bowtie. Bowtie index for the genome must be built before running this
        function.

        Notes
        -----
        The off-targets are stored in the `off_targets` attribute. Based on the
        off-targets, the following layers are added to the guide:

        - ot_sum_score: sum of the off-target scores - the lower the better
        - ot_cfd_score_mean: mean of the CFD scores of the off-targets
        - ot_cfd_score_max: max CFD scores of the off-targets
        - ot_cfd_score_sum: sum CFD scores of the off-targets

        Parameters
        ----------
        genome : Genome
            Genome object with the Bowtie index built
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

            # Load scoring matrices for CFD score calculation
            mm_scores, pam_scores = load_cfd_scoring_matrix()

            # Calculate CFD scores
            cfd_scores = calculate_cfd_score(
                self, self.off_targets, mm_scores, pam_scores
            )

            for ix, cfd in enumerate(cfd_scores.tolist()):
                self.off_targets[ix]["cfd_score"] = cfd

            self.add_layer("ot_cfd_score_mean", cfd_scores.mean())
            self.add_layer("ot_cfd_score_max", cfd_scores.max())
            self.add_layer("ot_cfd_score_sum", cfd_scores.sum())

        else:
            raise ValueError("Bowtie index is not built for the genome / locus.")

    @property
    def off_targets_string(self):
        """Returns a string representation of the off-targets.

        The string representatio captures the number of off-targets with certain number
        of mismatches: n0|n1|n2|n3|n4|n5 (total), where n0 is the number of off-targets
        with 0 mismatches, n1 is the number of off-targets with 1 mismatch, etc.

        For example, if there are 3 off-targets with 0 mismatches, 2 with 1 mismatch,
        1 with 2 mismatches, 0 with 3 mismatches, 5 with 4 mismatches and 1 with 5 mismatches
        the string representation will be "3|2|1|0|5|1 (13)". In the parenthesis, the total
        number of off-targets is given.

        Returns
        -------
        str
            String representation of off-targets.
        """

        # if there are no off-targets, return 0
        if not self.off_targets:
            return "0"

        # count the number of off-targets with certain number of mismatches
        counts = Counter(
            [
                len(ot["mismatches"].keys()) - self.guide_pam.count("N")
                for ot in self.off_targets
            ]
        )

        # create a string representation
        ot_str = "|".join([str(counts[i]) for i in range(max(counts.keys()) + 1)])
        ot_str += f" ({len(self.off_targets)})"

        return ot_str

    @property
    def layers(self):
        return self._layers

    def add_layer(self, name, layer_data):
        """_summary_

        Parameters
        ----------
        name : str
            _description_
        layer_data : float
            _description_
        """
        self._layers[name] = layer_data

    def layer(self, key):
        """_summary_

        Parameters
        ----------
        key : _type_
            _description_

        Returns
        -------
        _type_
            _description_

        Raises
        ------
        ValueError
            _description_
        """
        if key in self._layers.keys():
            return self._layers[key]
        else:
            raise ValueError(f"Layer {key} was not added to the gRNA.")

    def add_azimuth_score(self, model_file="V3_model_nopos.pickle"):
        """Apply Azimuth score to a list of guides.

        Azimuth is a machine learning-based predictive modelling of CRISPR/Cas9 guide efficiency.
        Sometimes its reffered to as Doench 2016 score.

        Described in https://doi.org/10.1038/nbt.3437 (Doench et al., 2016)

        Returns
        -------
        float
            Azimuth score.
        """

        # Azimuth model requires 30bp sequences
        if (
            len(self.guide_long_seq) == 30
            and "N" not in self.guide_long_seq
            and self.guide_pam.endswith("GG")
            and len(self.guide_pam) == 3
        ):
            score_azimuth = float(
                model_comparison.predict(
                    np.array([self.guide_long_seq]),
                    length_audit=True,
                    pam_audit=True,
                    model_file=os.path.dirname(os.path.realpath(__file__))
                    + "/data/"
                    + model_file,
                )[0]
            )
        else:
            score_azimuth = 0.0

        # add layer to gRNA
        self.add_layer("azimuth_score", score_azimuth)

        return score_azimuth
