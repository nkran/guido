import subprocess
import tempfile

import numpy as np

from .helpers import rev_comp


def calculate_ot_sum_score(
    off_targets, ot_mismatch_weights=[10, 5, 4, 3, 1], pam="NGG"
):
    """Calculate sum score for a list of off-targets."""
    return sum(
        [
            ot_mismatch_weights[len(ot["mismatches"]) - pam.count("N")]
            for ot in off_targets
        ]
    )


def get_off_targets_string(off_targets, pam="NGG"):
    """Get a string representation of a list of off-targets."""
    count_matrix = [0, 0, 0, 0, 0]

    for ot in off_targets:
        count_matrix[len(ot["mismatches"]) - pam.count("N")] += 1

    return "|".join([str(c) for c in count_matrix])


def calculate_cfd_score(guide, offtargets, mm_scores, pam_scores):
    """Calculate CFD score for a guide and a list of off-targets.

    Adapted from Doench et al. 2016, https://doi.org/10.1038/nbt.3437
    """

    scores = []
    wt = guide.guide_seq[: -len(guide.guide_pam)]

    for ot in offtargets:
        pam = guide.guide_pam[1:]
        sg = ot["seq"][:20]

        score = 1
        sg = sg.replace("T", "U")
        wt = wt.replace("T", "U")
        s_list = list(sg)
        wt_list = list(wt)

        for i, sl in enumerate(s_list):
            if wt_list[i] == sl:
                score *= 1
            else:
                key = (
                    "r" + wt_list[i] + ":d" + rev_comp(sl, rna=True) + "," + str(i + 1)
                )
                score *= mm_scores[key]
        score *= pam_scores[pam]
        scores.append(score)

    return np.array(scores)


def _parse_mismatches(mismatches, strand, seq_len):
    """Parse mismatches from bowtie output.

    Returns a dictionary with mismatches positions as keys and
    mismatched nucleotides as values.
    """
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


def _hit_is_valid(mismatches, non_arbitrary_positions):
    return all(
        [False if pos in non_arbitrary_positions else True for pos in mismatches.keys()]
    )


def run_bowtie(
    guides,
    pam="NGG",
    core_length=10,
    core_mismatches=2,
    total_mismatches=4,
    genome_index_path=None,
    threads=1,
    bowtie_path="bin/bowtie/",
):
    """Run bowtie to find off-targets for a list of guides.

    # TODO: add description of mismatch settings

    Parameters
    ----------
    guides : Guide[]
        List of guides to find off-targets for.
    pam : str, optional
        PAM sequence, by default "NGG"
    core_length : int, optional
        gRNA core length, by default 10
    core_mismatches : int, optional
        Allowed mismatches in core gRNA sequence, by default 2
    total_mismatches : int, optional
        Total allowed mismatches. Max 4, by default 4
    genome_index_path : Genome, optional
        Genome object to search off-targets in, by default None
    threads : int, optional
        Number of threads to run bowtie with, by default 1
    bowtie_path : str, optional
        Path to bowtie binary, by default "bin/bowtie/"

    Returns
    -------
    list
        List of off-targets for each guide.
    """

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

        stdout, _ = rproc.communicate()

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
                    g_seq = rev_comp(g_seq)
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

                    mm_seq = "".join(
                        [
                            mismatches[i + 1]
                            if (i + 1) in mismatches.keys()
                            else t_seq[i]
                            for i, _ in enumerate(g_seq)
                        ]
                    )

                    off_target["ix"] = ix
                    off_target["mismatches"] = mismatches
                    off_target["mismatches_string"] = mm_string
                    off_target["chromosome"] = t_chrom
                    off_target["start"] = int(t_start)
                    off_target["strand"] = t_strand
                    off_target["seq"] = mm_seq

                    off_targets[ix].append(off_target)

    return off_targets
