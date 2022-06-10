from shutil import which

import numpy as np
import pandas as pd


def rev_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(
        [complement[base] if base in complement.keys() else "N" for base in seq[::-1]]
    )


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    # from whichcraft import which
    return which(name) is not None


def guides_to_dataframe(guides):
    guides_ext = []
    for g in guides:
        g = g.__dict__
        for lk, lv in g["_layers"].items():
            if isinstance(lv, np.ndarray):
                g[lk] = lv.mean()
            else:
                g[lk] = lv
        guides_ext.append(g)

    return pd.DataFrame(guides_ext)


def guides_to_csv(guides, file):
    guides_df = guides_to_dataframe(guides).drop(
        ["mmej_patterns", "off_targets", "_layers"]
    )
    guides_df.to_csv(file)

    return None
