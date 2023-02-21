import os
import pickle
from shutil import which

import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader

EXCLUDED_ATTRS = [
    "locus_seq",
    "guide_long_seq",
    "relative_cut_pos",
    "mmej_patterns",
    "off_targets",
    "_layers",
]


def rev_comp(seq, rna=False):

    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    if rna:
        complement = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A", "-": "-"}
    return "".join(
        [complement[base] if base in complement.keys() else "N" for base in seq[::-1]]
    )


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    # from whichcraft import which
    return which(name) is not None


def load_cfd_scoring_matrix():
    """Load the cfd scoring matrix."""
    mm_scores = pickle.load(
        open(
            os.path.dirname(os.path.realpath(__file__)) + "/data/mismatch_score.pkl",
            "rb",
        )
    )
    pam_scores = pickle.load(
        open(os.path.dirname(os.path.realpath(__file__)) + "/data/PAM_scores.pkl", "rb")
    )

    return mm_scores, pam_scores


def _guide_instance_to_dict(guide):
    """Convert a class instance to a dictionary."""
    guide_dict = {}
    for attr, value in guide.__dict__.items():
        if attr not in EXCLUDED_ATTRS:
            guide_dict[attr] = value

    for lk, lv in guide.layers.items():
        if isinstance(lv, np.ndarray):
            guide_dict[lk] = lv.mean()
        else:
            guide_dict[lk] = lv

    return guide_dict


def _guides_to_dataframe(guides):
    """Convert a list of guides to a pandas dataframe."""

    guides_ext = [_guide_instance_to_dict(g) for g in guides]
    df = pd.DataFrame(guides_ext)
    df.index = df["id"]

    return df


def _guides_to_csv(guides, file):
    """Write a list of guides to a csv file."""

    guides_df = _guides_to_dataframe(guides)
    guides_df.to_csv(file)

    return None


def _guides_to_bed(guides, file):
    """Write a list of guides to a bed file."""

    guides_df = _guides_to_dataframe(guides).loc[
        :, ["guide_chrom", "guide_start", "guide_end", "id", "rank", "guide_strand"]
    ]
    guides_df.to_csv(file, sep="\t", header=False, index=False)
    return None


def _guides_detailed_table(guides, file):
    """Write a detailed list of guides to a text file."""

    file_loader = FileSystemLoader(
        os.path.dirname(os.path.realpath(__file__)) + "/templates/"
    )
    env = Environment(loader=file_loader)
    env.filters["rev_comp"] = rev_comp
    template = env.get_template("guide_details.j2")
    output_details = template.render(guides=guides)

    with open(file, "w") as f:
        f.write(output_details)


# def render_mmej_table(mmej_patterns):

#     for mmej_pattern in mmej_patterns:
#         left_seq = mmej_pattern["left"][20:]
#         right_seq = mmej_pattern["right"][:-20]
#         pos_mmej = right_seq.index(mmej_pattern["pattern"])
#         right_del = right_seq[:pos_mmej]
#         right_mmej = right_seq[pos_mmej : pos_mmej + len(mmej_pattern["pattern"])]
#         right_rest = right_seq[pos_mmej + len(mmej_pattern["pattern"]) :]
#         mmej_sequence = [
#             left_seq,
#             html.Span([right_del, html.U(right_mmej), right_rest]),
#         ]

#         pattern_layout = html.Div(
#             [
#                 html.Div(
#                     html.Span(mmej_sequence, className="mmej-pattern"),
#                     className="col-7",
#                 ),
#                 html.Div(
#                     html.Span(
#                         round(mmej_pattern["pattern_score"], 2),
#                         className="mmej-pattern",
#                     ),
#                     className="col",
#                 ),
#                 html.Div(
#                     html.Span(mmej_pattern["frame_shift"], className="mmej-pattern"),
#                     className="col",
#                 ),
#                 html.Div(
#                     html.Span(mmej_pattern["deletion_seq"], className="mmej-pattern"),
#                     className="col-2",
#                 ),
#             ],
#             className="row border-bottom border-1",
#         )
