import tempfile
import subprocess
from io import StringIO
import pandas as pd
from shutil import which

from pathlib import Path
from typing import List, Dict, Tuple, Union, Optional

# TODO typing


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    # from whichcraft import which
    return which(name) is not None


def run_bowtie(
    guides: list,
    pam: str,
    genome_index_path: Union[Path, str],
    max_offtargets: int = 10,
    threads: int = 1,
    bowtie_path: Optional[Union[Path, str]] = "",
) -> tuple:
    mismatches = 3
    # create temporary file for bowtie input
    with tempfile.NamedTemporaryFile(mode="w+t", prefix="guido_") as temp:
        temp.write(
            "\n".join(
                [
                    ">{}|{}|{}|{}|{}\n{}".format(
                        i,
                        G.guide_chrom,
                        G.guide_start,
                        G.guide_end,
                        G.guide_seq,
                        G.guide_seq,
                    )
                    for i, G in enumerate(guides)
                ]
            )
        )
        temp.flush()

        print(
            "\n".join(
                [
                    ">{}|{}|{}|{}|{}\n{}".format(
                        i,
                        G.guide_chrom,
                        G.guide_start,
                        G.guide_end,
                        G.guide_seq,
                        G.guide_seq,
                    )
                    for i, G in enumerate(guides)
                ]
            ))
        
        # run bowtie alignment
        bowtie_command = f"{bowtie_path}bowtie -p {threads} -v {mismatches} --sam --sam-nohead -k {max_offtargets} -x {genome_index_path} -f {temp.name}"
        print(bowtie_command)
        rproc = subprocess.Popen(
            bowtie_command.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = rproc.communicate()

    # output bowtie to pandas dataframe
    targets = pd.read_csv(
        StringIO(stdout),
        sep="\t",
        names=list(range(14)),
        header=None,
        index_col=False,
        converters={3: int},
    )
    
    # split sequence identification
    nsplit = targets[0].str.split("|", n=4, expand=True)
    nsplit.columns = ["id", "chrom", "start", "end", "seq"]

    # add it back to the target df
    targets = pd.concat([targets, nsplit], axis=1)
    targets["id"] = targets["id"].astype(int)
    targets["start"] = targets["start"].astype(int)

    # check if it's on-target
    on_targets_idx = targets[
        (targets[2] == targets["chrom"]) & (targets[3] == targets["start"])
    ].index

    # drop on-targets
    targets = targets.drop(on_targets_idx)

    # extract the number of mismatches in the off-target
    targets["mm"] = targets.apply(lambda x: x[13].split(":")[-1], axis=1)

    # count off-targets for a given guide
    mismatches_count = (
        targets.groupby(["id"])["mm"].value_counts().unstack().fillna(0).reset_index()
    )

    # count mismatched off-targets
    for m in ["0", "1", "2", "3"]:
        if m not in mismatches_count:
            mismatches_count[m] = 0
    mismatches_count["id"] = mismatches_count["id"].astype(int)

    # add mismatch info to the guide dict
    guides_offtargets = {}

    for ix, x in mismatches_count.iterrows():
        i = int(x.id)
        counts = x[["0", "1", "2", "3"]].to_dict()
        guides_offtargets[i] = {}
        guides_offtargets[i]["offtargets_dict"] = counts
        guides_offtargets[i][
            "offtargets_str"
        ] = "{:0.0f}|{:0.0f}|{:0.0f}|{:0.0f}".format(*counts.values())
        guides_offtargets[i]["offtargets_n"] = sum(counts.values())

    return (guides_offtargets, targets, stdout)
