import tempfile
import subprocess

from pathlib import Path
from typing import Union

from guido.helpers import rev_comp

# TODO typing
# TODO assuming PAM has at least one arbitrary nucleotide


def _parse_mismatches(mismatches: str, strand: str, seq_len: int):
    mm_split = mismatches.split(",")
    mm_dict = {}
    for m in mm_split:
        pos, c = m.split(":")
        c = c.split(">")

        if strand == "-":
            mm_dict[seq_len - int(pos)] = rev_comp(c[0])
        else:
            mm_dict[seq_len - int(pos)] = c[0]
    return mm_dict


def _hit_is_valid(mismatches: dict, non_arbitrary_positions: list):
    return all(
        [False if pos in non_arbitrary_positions else True for pos in mismatches.keys()]
    )


def run_bowtie(
    guides: list,
    pam: str = "NGG",
    core_length: int = 10,
    core_mismatches: int = 0,
    total_mismatches: int = 4,
    genome_index_path: Union[str, Path] = None,
    threads: int = 1,
    bowtie_path: str = "bin/bowtie/",
) -> dict:

    pam_mismatches = rev_comp(pam).count("N")
    pam_length = len(pam)
    off_targets = {}

    with tempfile.NamedTemporaryFile(mode="w+t", prefix="guido_") as temp:
        for i, G in enumerate(guides):
            g_seq = rev_comp(G.guide_seq[:-pam_length] + pam)
            temp.write(
                f">{i}|{G.guide_chrom}|{G.guide_start}|{G.guide_end}|{g_seq}|{G.guide_strand}"
            )
            temp.write(f"\n{g_seq}\n")
        temp.seek(0)

        if (core_mismatches + pam_mismatches) > 3:
            raise ValueError(
                f"The value for the parameter core_mismatches is not valid: {core_mismatches}"
            )

        if core_mismatches > total_mismatches:
            raise ValueError(
                f"The value for core_mismatches cannot be greater than total_mismatches: {core_mismatches} > {total_mismatches}"
            )

        bowtie_command = (
            f"{bowtie_path}bowtie -p {threads} --quiet -y "
            f"-n {core_mismatches + pam_mismatches} "
            f"-l {core_length + pam_length} "
            f"-e {(total_mismatches * 30) + 30 } -a "
            f"-x {genome_index_path} "
            f"-f {temp.name}"
        )

        rproc = subprocess.Popen(
            bowtie_command.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        stdout, stderr = rproc.communicate()

        non_arbitrary_positions = [
            21 + ix for ix, nucl in enumerate(pam) if nucl in ["A", "T", "G", "C"]
        ]

        # parse bowtie output
        for line in stdout.split("\n"):
            cols = line.split("\t")

            if len(cols) > 1:
                off_target = {}
                ix, g_chrom, g_start, _, g_seq, _ = cols[0].split("|")
                t_strand, t_chrom, t_start, t_seq, _, _, t_mismatches = cols[1:]
                ix = int(ix)

                if ix not in off_targets.keys():
                    off_targets[ix] = []

                if t_strand == "-":
                    t_strand = "+"

                else:
                    t_seq = rev_comp(t_seq)
                    t_strand = "-"

                mismatches = _parse_mismatches(t_mismatches, t_strand, len(t_seq))

                if (
                    mismatches
                    and _hit_is_valid(mismatches, non_arbitrary_positions)
                    and (g_chrom != t_chrom and int(g_start) != (int(t_start) + 1))
                ):
                    mm_string = "".join(
                        [
                            mismatches[i + 1] if (i + 1) in mismatches.keys() else "."
                            for i, _ in enumerate(g_seq)
                        ]
                    )

                    off_target["ix"] = ix
                    off_target["mismatches"] = mismatches
                    off_target["mismatches_string"] = mm_string
                    off_target["chromosome"] = t_chrom
                    off_target["start"] = int(t_start)
                    off_target["strand"] = t_strand

                    off_targets[ix].append(off_target)

    return off_targets
