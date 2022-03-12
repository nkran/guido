import re
import numpy as np
import pandas as pd
import pyranges
import mcdm

from typing import Union, Set
from pyfaidx import Sequence, Fasta

from guido.off_targets import run_bowtie, calculate_ot_sum_score
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
        """
        [summary]

        Parameters
        ----------
        sequence : Union[Sequence, str]
            [description]
        name : str, optional
            [description], by default None
        start : int, optional
            [description], by default 1
        end : int, optional
            [description], by default None
        genome : Genome, optional
            [description], by default None
        annotation : pd.DataFrame, optional
            [description], by default None
        """

        self.sequence = sequence
        self.chromosome = name
        self.start = start
        self.intervals = []
        self.genome = genome
        self.annotation = annotation
        self.pam = ""
        self.guides = []
        self._layers = {}

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

    # TODO: length of the locus
    def __len__(self):
        return self.end - self.start

    @property
    def lenght(self) -> int:
        return self.end - self.start

    def guide(self, ix: Union[str, int]) -> Union[Guide,None]:
        if isinstance(ix, str):
            for g in self.guides:
                if g.id == ix:
                    return g
        elif isinstance(ix, int):
            return self.guides[ix]
        else:
            raise ValueError('Provided gRNA index is not valid.')

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

    def _find_guides_in_interval(self, sequence, start, pam):

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

        pams_fw = [m.start() for m in re.finditer(rf"(?=({iupac_pam}))", fwd_seq)]
        pams_rv = [m.start() for m in re.finditer(rf"(?=({rev_iupac_pam}))", fwd_seq)]

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
        """
        [summary]

        Parameters
        ----------
        pam : str, optional
            [description], by default "NGG"
        min_flanking_length : int, optional
            [description], by default 75
        selected_features : Set[str], optional
            [description], by default {"all"}
        """

        # save searched PAM
        self.pam = pam

        # take default locus bounds
        self.intervals = [[self.start, self.end]]

        if "all" not in selected_features and self.annotation:
            locus_annotation = self.annotation.query(
                f"(Feature in @selected_features) & \
                (Chromosome == @self.chromosome) &  \
                (((Start >= @self.start) & (Start <= @self.end)) | \
                ((End >= @self.start) & (End <= @self.end)))"
            ).sort_values("Start")

            # split the locus into smaller loci defined by features
            if len(locus_annotation) > 0:
                self.intervals = self._flatten_intervals(
                    [
                        [f_start, f_end]
                        for f_start, f_end in locus_annotation[["Start", "End"]].values
                    ]
                )

        # search for guides in each locus
        for interval_start, interval_end in self.intervals:

            # check if there is sequence 30 bp upstrean
            rel_interval_start = interval_start - min_flanking_length - self.start
            rel_interval_end = interval_end + min_flanking_length - self.start

            if rel_interval_start < 0:
                rel_interval_start = 0

            if rel_interval_end > len(self.sequence.seq):
                rel_interval_end = len(self.sequence.seq)

            locus_sequence = self.sequence.seq[rel_interval_start:rel_interval_end]
            locus_guides = self._find_guides_in_interval(
                locus_sequence, interval_start - min_flanking_length, pam
            )

            # dont add shorter sequences, check that cut position is in the interval
            locus_guides = [
                g
                for g in locus_guides
                if len(g.guide_seq) == 23
                and g.absolute_cut_pos >= interval_start
                and g.absolute_cut_pos <= interval_end
            ]

            self.guides.extend(locus_guides)

        sorted_guides = []
        for ix, g in enumerate(sorted(self.guides, key=lambda g: g.absolute_cut_pos)):
            g.id = f"gRNA-{ix+1}"
            sorted_guides.append(g)

        self.guides = sorted_guides

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

        Raises
        ------
        ValueError
            [description]
        """
        if len(self.guides) == 0:
            raise ValueError("No gRNAs saved yet.")

        for G in self.guides:
            G.simulate_end_joining(n_patterns, length_weight)

    def find_off_targets(
        self, external_genome: Union[Genome, None] = None, **kwargs
    ) -> None:
        """
        [summary]

        Parameters
        ----------
        external_genome : Union[Genome, None], optional
            [description], by default None

        Raises
        ------
        ValueError
            [description]
        ValueError
            [description]
        ValueError
            [description]
        """
        if len(self.guides) == 0:
            raise ValueError("No gRNAs saved yet.")

        if external_genome:
            index_path = external_genome.bowtie_index
        elif self.genome:
            index_path = self.genome.bowtie_index
        else:
            raise ValueError("No genome / locus specified.")

        if index_path:
            guides_bowtie_offtargets = run_bowtie(
                guides=self.guides, pam=self.pam, genome_index_path=index_path, **kwargs
            )

            for ix, G in enumerate(self.guides):
                if ix in guides_bowtie_offtargets.keys():
                    G.off_targets = guides_bowtie_offtargets[ix]
                    G.add_layer("ot_sum_score", calculate_ot_sum_score(G.off_targets))

        else:
            raise ValueError("Bowtie index is not built for the genome / locus.")

    def _apply_clipped_layer_data(self, guides, layer_name, layer_data) -> None:
        """
        _summary_

        Parameters
        ----------
        guides : _type_
            _description_
        layer_name : _type_
            _description_
        layer_data : _type_
            _description_
        """
        if len(guides) > 0:
            for g in self.guides:
                ix = g.guide_start - self.start
                g.add_layer(layer_name, layer_data[ix : ix + 23])

    @property
    def layers(self):
        """
        _summary_

        Returns
        -------
        _type_
            _description_
        """
        return self._layers

    def add_layer(
        self, name: str, layer_data: np.ndarray, apply_to_guides: bool = True
    ):
        """
        _summary_

        Parameters
        ----------
        name : str
            _description_
        layer_data : np.ndarray
            _description_
        apply_to_guides : bool, optional
            _description_, by default True

        Raises
        ------
        ValueError
            _description_
        ValueError
            _description_
        """
        if layer_data.shape[0] == self.lenght:
            self._layers[name] = layer_data

            if apply_to_guides:
                if len(self.guides) > 0:
                    self._apply_clipped_layer_data(self.guides, name, layer_data)
                else:
                    raise ValueError("No guides to apply the data to.")
        else:
            raise ValueError("Layer and locus not the same lenght.")

    def layer(self, key: str) -> np.ndarray:
        """
        _summary_

        Parameters
        ----------
        key : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        return self._layers[key]

    def _guide_layers(self):
        layers = []
        for g in self.guides:
            layers.extend(g.layers.keys())
        return set(layers)

    def _prepare_alt_matrix(self, rank_layer_names: list, method=np.mean) -> np.ndarray:
        """
        [summary]

        Parameters
        ----------
        rank_layer_names : list
            [description]
        method : [type], optional
            [description], by default np.mean

        Returns
        -------
        np.ndarray
            [description]
        """
        locus_data = []
        for g in self.guides:
            guide_data = []
            for layer_name in rank_layer_names:
                if (
                    layer_name not in g.layers.keys()
                    and layer_name not in self.layers.keys()
                ):
                    raise ValueError(
                        f"Layer {layer_name} does not exist on `Locus` or `Guide` object."
                    )

                if (
                    layer_name not in g.layers.keys()
                    and layer_name in self.layers.keys()
                ):
                    self._apply_clipped_layer_data(
                        [g], layer_name=layer_name, layer_data=self.layer(layer_name)
                    )

                layer_data = g.layer(layer_name)

                if (
                    isinstance(layer_data, np.ndarray)
                    and layer_data.squeeze().ndim == 1
                ):
                    layer_data = method(layer_data)
                elif isinstance(layer_data, float) or isinstance(layer_data, int):
                    layer_data = float(layer_data)
                else:
                    layer_data = np.nan

                guide_data.append(layer_data)
            locus_data.append(guide_data)

        return np.array(locus_data)

    def rank_guides(
        self,
        rank_layer_names: list = None,
        is_benefit_layer: list = None,
        weight_vector: list = None,
        ranking_method: str = "TOPSIS",
        norm_method: str = "Vector",
        **kwargs,
    ) -> list:
        """
        [summary]

        Returns
        -------
        [type]
            [description]

        Raises
        ------
        ValueError
            [description]
        ValueError
            [description]
        """

        if len(self.guides) == 0:
            raise ValueError(
                f"No gRNAs to rank. Try running `find_guides()` method first."
            )

        # get layer names registered on individual guides
        guide_layer_names = sorted(
            list(
                set(
                    [
                        guide_layer_name
                        for g in self.guides
                        for guide_layer_name in g.layers.keys()
                    ]
                )
            )
        )

        # check if requested layer exists on guides or on locus
        if rank_layer_names:
            for layer in guide_layer_names:
                if layer not in self._layers.keys() and layer not in guide_layer_names:
                    print(self._layers.keys(), guide_layer_names)
                    raise ValueError(
                        f"Layer {layer} is not added to `Locus` or `Guide`."
                    )
        else:
            rank_layer_names = list(guide_layer_names) + list(self._layers.keys())

        x_matrix = self._prepare_alt_matrix(rank_layer_names=rank_layer_names)
        rank_scores = mcdm.rank(
            x_matrix,
            n_method=norm_method,
            w_vector=weight_vector,
            is_benefit_x=is_benefit_layer,
            s_method=ranking_method,
            alt_names=[g.id for g in self.guides],
            **kwargs,
        )

        for i, (g_id, rank_score) in enumerate(rank_scores):
            self.guide(g_id).rank_score = rank_score
            self.guide(g_id).rank = i+1
        
        return rank_scores


"""
Locus creation ------------------
"""


def _prepare_annotation(annotation_file_abspath, as_df=True):
    if annotation_file_abspath.suffix in [".gff3", ".gff"]:
        ann_db = pyranges.read_gff3(str(annotation_file_abspath), as_df=as_df)
        if as_df:
            ann_db =ann_db.rename(columns={"Name": "Exon"})  # type: ignore
    elif annotation_file_abspath.suffix in [".gtf"]:
        ann_db = pyranges.read_gtf(str(annotation_file_abspath), as_df=as_df)
        if as_df:
            ann_db = ann_db.rename(columns={"gene_id": "ID", "exon_number": "Exon"})  # type: ignore
    else:
        raise ValueError("Annotation file not recognised. Annotation file needs to be GFF3 or GTF formnat.")
    return ann_db


def locus_from_coordinates(
    genome: Genome, chromosome: str, start: int, end: int
) -> Locus:
    """
    [summary]

    Parameters
    ----------
    genome : Genome
        [description]
    chromosome : str
        [description]
    start : int
        [description]
    end : int
        [description]

    Returns
    -------
    Locus
        [description]
    """

    locus_sequence = Fasta(str(genome.genome_file_abspath)).get_seq(
        chromosome, start, end
    )

    if genome.annotation_file_abspath and genome.annotation_file_abspath.exists():
        ann_db = _prepare_annotation(genome.annotation_file_abspath, as_df=False)
        locus_annotation = ann_db.intersect(pyranges.PyRanges(chromosomes=[chromosome], starts=[start], ends=[end])).df
        
        if genome.annotation_file_abspath.suffix in [".gff3", ".gff"]:
            locus_annotation = locus_annotation.rename(columns={"Name": "Exon"})  # type: ignore
        elif genome.annotation_file_abspath.suffix in [".gtf"]:
            locus_annotation = locus_annotation.rename(columns={"gene_id": "ID", "exon_number": "Exon"})  # type: ignore
            
        
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
    """
    [summary]

    Parameters
    ----------
    sequence : str
        [description]
    sequence_name : str, optional
        [description], by default None

    Returns
    -------
    Locus
        [description]
    """
    return Locus(sequence=sequence, name=sequence_name)


def locus_from_gene(genome: Genome, gene_name: str) -> Locus:
    """
    [summary]

    Parameters
    ----------
    genome : Genome
        [description]
    gene_name : str
        [description]

    Returns
    -------
    Locus
        [description]

    Raises
    ------
    ValueError
        [description]
    ValueError
        [description]
    """
    if genome.annotation_file_abspath.exists():
        try:
            ann_db = _prepare_annotation(genome.annotation_file_abspath, as_df=True)
            gene_annotation = ann_db.query(f'ID == @gene_name & Feature == "gene"')
            chromosome = gene_annotation.Chromosome.values[0]
            start = int(gene_annotation.Start)
            end = int(gene_annotation.End)

            locus_annotation = ann_db.query(
                f"(ID == @gene_name) & (Chromosome == @chromosome) &  \
                  (((Start >= {start - 1}) & (Start <= {end + 1})) | \
                  ((End >= {start - 1}) & (End <= {end + 1})))"
            )

            locus_sequence = Fasta(str(genome.genome_file_abspath)).get_seq(
                chromosome, start, end
            )

            return Locus(
                genome=genome, sequence=locus_sequence, annotation=locus_annotation
            )

        except Exception:
            raise ValueError(f"Gene {gene_name} not found.")
    else:
        raise ValueError("Annotation file not valid.")  # TODO: fix message
