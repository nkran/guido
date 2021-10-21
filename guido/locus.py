from __future__ import annotations
from os import name
import re
from numpy.lib.function_base import iterable
import pandas as pd
import pyranges

from typing import Union, List, Set
from pyfaidx import Sequence, Fasta

from guido.off_targets import run_bowtie
from guido.guides import Guide
from guido.genome import Genome


class Locus:
    def __init__(
        self,
        sequence: Union[Sequence, str],
        name: str = None,
        start: int = 1,
        end: int = None,
        genome: Genome = None,
        annotation: pd.DataFrame = None,
        **kwargs,
    ) -> None:

        self.sequence = sequence
        self.chromosome = name
        self.start = start
        self.genome = genome
        self.annotation = annotation
        self.pam = ''
        self.guides = []

        if isinstance(sequence, Sequence):
            self.sequence = sequence
            self.start = sequence.start
            self.end = sequence.end
            self.chromosome = sequence.name
        else:
            self.sequence = Sequence(
                seq=sequence, start=start, end=end, name=name, **kwargs
            )
            if end:
                self.end = end
            else:
                self.end = self.start + len(
                    self.sequence.seq
                )  # TODO: test if length corresponds to given end

    def __iter__(self):
        return iterable(self.guides)

    def __len__(self):
        return len(self.guides)

    def _flatten_intervals(self, intervals):
        fi = []
        for start, end in intervals:
            if fi and fi[-1][1] >= start - 1:
                fi[-1][1] = max(fi[-1][1], end)
            else:
                if start < self.start:
                    start = self.start
                if end > self.end:
                    end = self.end
                fi.append([start, end])
        return fi

    def _find_guides_in_interval(self, sequence, start, pam, min_flanking_length):

        iupac_dict = {
            "A": ("A", "T"),
            "C": ("C", "G"),
            "G": ("G", "C"),
            "T": ("T", "A"),
            "R": ("[AG]", "[CT]"),
            "Y": ("[CT]", "[AG]"),
            "S": ("[GC]", "[GC]"),
            "W": ("[AT]", "[AT]"),
            "K": ("[GT]", "[AC]"),
            "M": ("[AC]", "[GT]"),
            "B": ("[CGT]", "[ACG]"),
            "D": ("[AGT]", "[ACT]"),
            "H": ("[ACT]", "[AGT]"),
            "V": ("[ACG]", "[CGT]"),
            "N": ("[ACGT]", "[ACGT]"),
        }

        iupac_pam = "".join([iupac_dict[letter][0] for letter in pam])
        rev_iupac_pam = "".join([iupac_dict[letter][1] for letter in pam[::-1]])
        fwd_seq = sequence.upper()

        pams_fw = [
            m.start()
            for m in re.finditer(rf"(?=({iupac_pam}))", fwd_seq)
            if m.start(0) - min_flanking_length > 0
            and m.end(0) + min_flanking_length < len(fwd_seq)
        ]
        pams_rv = [
            m.start()
            for m in re.finditer(rf"(?=({rev_iupac_pam}))", fwd_seq)
            if m.start(0) - min_flanking_length > 0
            and m.end(0) + min_flanking_length < len(fwd_seq)
        ]

        guides_fw = [
            Guide(
                sequence=fwd_seq,
                pam_position=pam_position,
                pam_len=len(pam),
                strand="+",
                chromosome=self.chromosome,
                start=start,
            )
            for pam_position in pams_fw
        ]
        guides_rv = [
            Guide(
                sequence=fwd_seq,
                pam_position=pam_position,
                pam_len=len(pam),
                strand="-",
                chromosome=self.chromosome,
                start=start,
            )
            for pam_position in pams_rv
        ]
        return guides_fw + guides_rv

    def find_guides(
        self,
        pam: str = "NGG",
        min_flanking_length: int = 75,
        selected_features: Set[str] = {"all"},
    ) -> None:

        # save searched PAM
        self.pam = pam
        
        # take default locus bounds
        loci_intervals = [[self.start, self.end]]
        
        if "all" not in selected_features:
            try:
                locus_annotation = self.annotation.query(
                    f"(Feature in @selected_features) & \
                    (Chromosome == @self.chromosome) &  \
                    (((Start >= @self.start) & (Start <= @self.end)) | \
                    ((End >= @self.start) & (End <= @self.end)))"
                ).sort_values("Start")

                # split the locus into smaller loci defined by features
                if len(locus_annotation) > 0:
                    loci_intervals = self._flatten_intervals(
                        [
                            [f_start, f_end]
                            for f_start, f_end in locus_annotation[
                                ["Start", "End"]
                            ].values
                        ]
                    )
            except:
                raise ValueError("Annotation not accessible.")

        # search for guides in each locus
        for locus_start, locus_end in loci_intervals:
            if locus_start != self.start and locus_end != self.end:
                locus_start -= self.start
                locus_end -= self.start
                locus_sequence = self.sequence.seq[locus_start - 1 : locus_end]
            else:
                locus_sequence = self.sequence.seq

            locus_guides = self._find_guides_in_interval(
                locus_sequence, locus_start, pam, min_flanking_length
            )
            self.guides.extend(locus_guides)

    def simulate_end_joining(
        self, n_patterns: int = 5, length_weight: int = 20
    ) -> None:

        if len(self.guides) == 0:
            raise ValueError("No gRNAs saved yet.")

        for G in self.guides:
            G.simulate_end_joining(n_patterns, length_weight)

    def find_off_targets(
        self, external_genome: Union[Genome, None] = None, **kwargs
    ) -> None:

        if len(self.guides) == 0:
            raise ValueError("No gRNAs saved yet.")

        if external_genome:
            index_path = external_genome.bowtie_index
        elif self.genome:
            index_path = self.genome.bowtie_index
        else:
            raise ValueError("No genome / locus specified.")

        if index_path:
            guides_offtargets, targets = run_bowtie(
                guides=self.guides, pam=self.pam, genome_index_path=index_path, **kwargs
            )

            if len(guides_offtargets.items()) > 0:
                for ix, g in guides_offtargets.items():
                    self.guides[ix].offtargets_str = g["offtargets_str"]
                    self.guides[ix].offtargets_n = g["offtargets_n"]
                    self.guides[ix].offtargets_dict = g["offtargets_dict"]
        else:
            raise ValueError("Bowtie index is not built for the genome / locus.")


def _prepare_annotation(annotation_file_abspath):
    if annotation_file_abspath.suffix in [".gff3", ".gff"]:
        ann_db = pyranges.read_gff3(str(annotation_file_abspath), as_df=True)
        ann_db.rename(columns={"Name": "Exon"}, inplace=True)
    elif annotation_file_abspath.suffix in [".gtf"]:
        ann_db = pyranges.read_gtf(str(annotation_file_abspath), as_df=True)
        ann_db.rename(columns={"gene_id": "ID", "exon_number": "Exon"}, inplace=True)
    else:
        raise ValueError("Annotation file not recognised")  # TODO: fix message
    return ann_db


def load_from_coordinates(
    genome: Genome, chromosome: str, start: int, end: int
) -> Locus:

    locus_sequence = Fasta(str(genome.genome_file_abspath)).get_seq(
        chromosome, start, end
    )

    if genome.annotation_file_abspath.exists():
        ann_db = _prepare_annotation(genome.annotation_file_abspath)
        locus_annotation = ann_db.query(
            f"(Chromosome == @chromosome) &  \
                (((Start >= {start - 1}) & (Start <= {end + 1})) | \
                ((End >= {start - 1}) & (End <= {end + 1})))"
        )

        if len(locus_annotation) > 0:
            return Locus(
                genome=genome, sequence=locus_sequence, annotation=locus_annotation
            )
        else:
            return Locus(
                genome=genome, sequence=locus_sequence, annotation=None
            )  # TODO: refactor
    else:
        return Locus(genome=genome, sequence=locus_sequence, annotation=None)


def locus_from_sequence(sequence: str, sequence_name: str = None) -> Locus:
    return Locus(sequence=sequence, name=sequence_name)


def load_from_gene(genome: Genome, gene_name: str) -> Locus:
    if genome.annotation_file_abspath.exists():
        try:
            ann_db = _prepare_annotation(genome.annotation_file_abspath)
            locus_annotation = ann_db.query(f'ID == @gene_name & Feature == "gene"')
            chromosome = locus_annotation.Chromosome.values[0]
            start = int(locus_annotation.Start)
            end = int(locus_annotation.End)

            locus_sequence = Fasta(str(genome.genome_file_abspath)).get_seq(
                chromosome, start, end
            )

            return Locus(
                genome=genome, sequence=locus_sequence, annotation=locus_annotation
            )

        except Exception:
            raise ValueError(f"Gene {gene_name} not found.")
    else:
        raise ValueError("Annotation file not recognised")  # TODO: fix message
