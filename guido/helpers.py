from shutil import which


def rev_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(
        [complement[base] if base in complement.keys() else "N" for base in seq[::-1]]
    )


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    # from whichcraft import which
    return which(name) is not None
