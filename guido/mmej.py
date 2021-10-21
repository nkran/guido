import re
import math
import itertools
import numpy as np


def generate_mmej_patterns(rel_cut_position, sequence, length_weight):
    """
    Start with predefined k-mer length and extend it until it finds more
    than one match in sequence.
    """

    # split the sequence at the cut site
    left_seq = sequence[:rel_cut_position]
    right_seq = sequence[rel_cut_position:]

    # kmers list
    kmers = []
    # combination of positions
    all_combinations = []
    # k-mer starting length
    min_kmer_length = 2

    # expand k-mer length
    for k in list(range(min_kmer_length, len(left_seq)))[::-1]:

        # iterate through sequence
        for i in range(len(left_seq) - k + 1):
            kmer = left_seq[i : i + k]

            # check if the k-mer is in both sequences and
            # hasn't been found before
            if kmer in right_seq and kmer not in kmers:

                # find positions of kmer
                p = re.compile(kmer)

                # find positions of patterns in each sequence
                left_positions = [(m.start(), m.end()) for m in p.finditer(left_seq)]
                right_positions = [(m.start(), m.end()) for m in p.finditer(right_seq)]

                # get combinations
                pos_combinations = [
                    c[0] + c[1]
                    for c in list(itertools.product(left_positions, right_positions))
                ]

                all_combinations.extend(pos_combinations)

                # save kmer
                kmers.append(kmer)

    if all_combinations:

        combinations_array = np.array(all_combinations)
        m1 = combinations_array[:, 0, None] <= combinations_array[:, 0]
        m2 = combinations_array[:, 1, None] >= combinations_array[:, 1]

        m3 = combinations_array[:, 2, None] <= combinations_array[:, 2]
        m4 = combinations_array[:, 3, None] >= combinations_array[:, 3]

        m12 = m1 & m2
        m34 = m3 & m4

        subpattern_rows = np.triu((m12 & m34), 1).any(0)

        # all mmej patterns
        mmej_patterns = []

        # generate microhomology for every combination
        for combination in combinations_array[~subpattern_rows]:

            # save output to dict
            pattern_dict = {}

            left_start, left_end, right_start, right_end = combination
            pattern_seq = left_seq[left_start:left_end]
            pattern_GC = pattern_seq.count("G") + pattern_seq.count("C")

            # deletion length
            left_deletion_length = len(left_seq) - left_start
            right_deletion_length = right_start
            deletion_length = left_deletion_length + right_deletion_length
            deletion_seq = left_seq[left_start:] + right_seq[:right_start]

            # score pattern
            length_factor = round(1 / math.exp((deletion_length) / (length_weight)), 3)
            pattern_score = (
                100
                * length_factor
                * ((len(pattern_seq) - pattern_GC) + (pattern_GC * 2))
            )

            # frame shift
            if deletion_length % 3 == 0:
                frame_shift = "-"
            else:
                frame_shift = "+"

            # create dictionary
            pattern_dict["left"] = left_seq[:left_start] + "-" * left_deletion_length
            pattern_dict["left_seq"] = left_seq[:left_start]
            pattern_dict["left_seq_position"] = left_start
            pattern_dict["right"] = (
                "+" * right_deletion_length + right_seq[right_start:]
            )
            pattern_dict["right_seq"] = right_seq[right_start:]
            pattern_dict["right_seq_position"] = right_start + len(deletion_seq)
            pattern_dict["pattern"] = pattern_seq
            pattern_dict["pattern_len"] = len(pattern_seq)
            pattern_dict["pattern_score"] = pattern_score
            pattern_dict["deletion_seq"] = deletion_seq
            pattern_dict["frame_shift"] = frame_shift

            # add to list
            mmej_patterns.append(pattern_dict)

        return mmej_patterns
