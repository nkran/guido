import re

import allel
import mcdm
import numpy as np
import pyranges
from pyfaidx import Fasta, Sequence

from guido.guides import Guide
from guido.helpers import _guides_to_bed, _guides_to_csv, _guides_to_dataframe
from guido.off_targets import calculate_ot_sum_score, run_bowtie


class Locus:
    def __init__(
        self,
        sequence,
        name=None,
        start=1,
        end=None,
        genome=None,
        annotation=None,
        **kwargs,
    ):
        """Locus object that represents a genomic locus in which gRNAs are
        searched.

        Parameters
        ----------
        sequence : Sequence, str
            Nucleotide sequence. Can be defined as `pyfaidx.Sequence` or `str`.
        name : str, optional
            Name of the sequence, contig or chromosome, by default None.
        start : int, optional
            Starting position of the locus. It should be 1-based. By default 1.
        end : int, optional
            Ending position of the locus.
        genome : Genome, optional
            Genome object that also includes the genomic locus defined here. It is
            used by default to search for gRNA off-targets.
        annotation : pd.DataFrame, optional
            Annotation table that provides genomic features and is used to annotate
            gRNAs that are found in the locus.

        Examples
        ----------

        ```
        >>> import guido
        >>> seq = "TTATCATCCACTCTGACGGGTGGTATTGCGCAACTCCACGCCATCAAACATGTTCAGATTATGCAATCGTGAGTATTCGTTGACCACCGCTTGACCTGTGT"
        >>> l = guido.Locus(
        ...     sequence=seq, name="AgamP4_2R", start=48714554, end=48714654
        ... )
        >>> l.find_guides()
        >>> l.guides
        [gRNA-1(CGCAATACCACCCGTCAGAGTGG|AgamP4_2R:48714561-48714584|-|),
         gRNA-2(TTATCATCCACTCTGACGGGTGG|AgamP4_2R:48714554-48714577|+|),
         gRNA-3(TCTGAACATGTTTGATGGCGTGG|AgamP4_2R:48714589-48714612|-|),
         gRNA-4(CATAATCTGAACATGTTTGATGG|AgamP4_2R:48714594-48714617|-|)]
        ```
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
            if end and isinstance(end, int):
                self.end = end
            elif not end:
                self.end = self.start + len(
                    self.sequence.seq
                )  # TODO: test if length corresponds to given end

        if self.end:
            self.length = self.end - self.start
        else:
            self.length = len(self.sequence)

    def guide(self, ix):
        """Fetch a guide from the locus by it's index or name.

        Parameters
        ----------
        ix : str or int
            Index of the gRNA.

        Returns
        -------
        g: Guide
            Guide object representing a gRNA

        Example
        -------
        ```
        >>> import guido
        >>> seq = "TTATCATCCACTCTGACGGGTGGTATTGCGCAACTCCACGCCATCAAACATGTTCAGATTATGCAATCGTGAGTATTCGTTGACCACCGCTTGACCTGTGT"
        >>> l = guido.Locus(
        ...     sequence=seq, name="AgamP4_2R", start=48714554, end=48714654
        ... )
        >>> l.find_guides()
        >>> l.guide("gRNA-1")
        gRNA-1(CGCAATACCACCCGTCAGAGTGG|AgamP4_2R:48714561-48714584|-|)
        >>> l.guide(0)
        gRNA-1(CGCAATACCACCCGTCAGAGTGG|AgamP4_2R:48714561-48714584|-|)
        ```
        """

        if isinstance(ix, str):
            for g in self.guides:
                if g.id == ix:
                    return g
            else:
                raise ValueError("Provided gRNA index is not valid.")
        elif isinstance(ix, int):
            return self.guides[ix]
        else:
            raise ValueError("Provided gRNA index is not valid.")

    def _flatten_intervals(self, intervals):
        """Flattens overlapping intervals into an union."""
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
        """Finds guides in interval."""
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
        pam="NGG",
        min_flanking_length=0,
        selected_features={"all"},
    ):
        """Find gRNAs in the locus.

        Parameters
        ----------
        pam : str, optional
            gRNA PAM sequence, by default "NGG"
        min_flanking_length : int, optional
            Defines flanking region from the locus where gRNAs are ignored. By default
            0, however `simulate_end_joining()` requires flanking region of 75 bp to
            simulate MMEJ.
        selected_features : str, optional
            Limit gRNA search on only specified genomic features. Features are defined
            in the provided genome annotation file. By default {"all"}
        """

        # save searched PAM
        self.pam = pam

        # take default locus bounds
        self.intervals = [[self.start, self.end]]

        if "all" not in selected_features and self.annotation:
            locus_annotation = self.annotation.query(
                "(Feature in @selected_features) & \
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

        return sorted_guides

    def simulate_end_joining(self, n_patterns=5, length_weight=20):
        """Simulate end-joining and find MMEJ deletion patterns for each gRNA.

        Microhomology scores are calculated based on proposed scoring model descibed by
        Bae et al. 2014.

        Parameters
        ----------
        n_patterns : int, optional
            Number of top scored MMEJ deletion patterns reported. By default 5.
        length_weight : int, optional
            Length weight parameter used in MMEJ scoring as defined by Bae et al. 2015.
            By default 20.
        """

        if len(self.guides) == 0:
            raise ValueError("No gRNAs saved yet.")

        for G in self.guides:
            G.simulate_end_joining(n_patterns, length_weight)

    def find_off_targets(self, external_genome=None, **kwargs):
        """Find off-targets in the genome for each gRNA.

        Parameters
        ----------
        external_genome : Genome, optional
            If provided, off-target search is performed in the external genome rather
            than in the genome which Locus is a part of. By default None.
        """

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

        return guides_bowtie_offtargets

    def _apply_clipped_layer_data(self, guides, layer_name, layer_data):
        """_summary_"""
        if len(guides) > 0:
            for g in self.guides:
                ix = g.guide_start - self.start
                g.add_layer(layer_name, layer_data[ix : ix + 23])

    # TODO: move to util functions
    def _get_guide_regions(self, guide):
        if guide.guide_strand == "-":
            pam_pos = (guide.guide_start, guide.guide_start + 1)
            seed_pos = (guide.guide_start + 3, guide.guide_start + 13)
            seed_small_pos = (guide.guide_start + 3, guide.guide_start + 7)
        else:
            pam_pos = (guide.guide_end - 1, guide.guide_end)
            seed_pos = (guide.guide_end - 13, guide.guide_end - 3)
            seed_small_pos = (guide.guide_end - 7, guide.guide_end - 3)

        guide_pos = (guide.guide_start, guide.guide_end)
        guide_regions = [pam_pos, seed_small_pos, seed_pos, guide_pos]

        return guide_regions

    # TODO: move to util functions
    def _guide_sequence_diversity(self, guide, g, pos):
        guide_regions = self._get_guide_regions(guide)
        regions_vals = []
        for r in guide_regions:
            try:
                region_loc = pos.locate_range(r[0], r[1])
                region_pos = pos[region_loc]
                region_ac = g[region_loc].count_alleles()
                pi = allel.sequence_diversity(region_pos, region_ac)
                regions_vals.append(pi)
            except Exception:
                regions_vals.append(0)

        return regions_vals

    def _guide_alt_ac(self, guide, g, pos):
        guide_regions = self._get_guide_regions(guide)
        regions_vals = []
        for r in guide_regions:
            try:
                region_loc = pos.locate_range(r[0], r[1])
                regions_vals.append(g[region_loc].count_alleles()[:, 1:].sum())
            except Exception:
                regions_vals.append(0)

        return regions_vals

    def _guide_n_variants(self, guide, g, pos):
        guide_regions = self._get_guide_regions(guide)
        regions_vals = []
        for r in guide_regions:
            try:
                region_loc = pos.locate_range(r[0], r[1])
                regions_vals.append(g[region_loc].n_variants)
            except Exception:
                regions_vals.append(0)

        return regions_vals

    def _apply_variation_layer_data(
        self, guides, layer_name, layer_genotype_data, layer_pos
    ):
        """_summary_"""
        guide_regions = ["pam", "seed", "small_seed", "guide"]

        if len(guides) > 0:
            for g in self.guides:
                guide_sequence_diversity = self._guide_sequence_diversity(
                    g, layer_genotype_data, layer_pos
                )
                guide_alt_ac = self._guide_alt_ac(g, layer_genotype_data, layer_pos)
                guide_n_variants = self._guide_n_variants(
                    g, layer_genotype_data, layer_pos
                )
                for i in range(len(guide_regions)):
                    g.add_layer(
                        f"var_{layer_name}_{guide_regions[i]}_pi",
                        guide_sequence_diversity[i],
                    )
                    g.add_layer(
                        f"var_{layer_name}_{guide_regions[i]}_alt_ac", guide_alt_ac[i]
                    )
                    g.add_layer(
                        f"var_{layer_name}_{guide_regions[i]}_variants",
                        guide_n_variants[i],
                    )

    @property
    def layers(self):
        """_summary_

        Returns
        -------
        _type_
            _description_
        """
        return self._layers

    def add_layer(
        self,
        name,
        layer_data,
        layer_pos=None,
        apply_to_guides=True,
        is_variation=False,
    ):
        """Adds a layer with the data to the locus.

        TODO: add examples

        Parameters
        ----------
        name : str
            Name of the layer
        layer_data : np.ndarray
            Layer data. Needs to be the same shape as the locus.
        apply_to_guides : bool, optional
            Apply layer data to gRNAs when adding it to the locus. By default True.

        Raises
        ------
        ValueError
            Layer data and locus are not the same size
        ValueError
            No gRNAs in the locus to apply the data to.
        """
        self._layers[name] = layer_data

        if apply_to_guides:
            if len(self.guides) > 0:
                if is_variation and layer_pos:
                    self._apply_variation_layer_data(
                        self.guides, name, layer_data, layer_pos
                    )
                else:
                    if layer_data.shape[0] == self.length:
                        self._apply_clipped_layer_data(self.guides, name, layer_data)
                    else:
                        raise ValueError("Layer and locus not the same lenght.")
        else:
            raise ValueError("No guides to apply the data to.")

    def _guide_layers(self):
        layers = []
        for g in self.guides:
            layers.extend(g.layers.keys())
        return set(layers)

    def _prepare_alt_matrix(self, rank_layer_names, method=np.mean):
        """Prepares numerical matrix with the gRNA layer data to be used later
        in the ranking.

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
                        [g], layer_name=layer_name, layer_data=self._layers[layer_name]
                    )

                layer_data = g._layers[layer_name]

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
        layer_names=None,
        layer_is_benefit=None,
        weight_vector=None,
        ranking_method="TOPSIS",
        norm_method="Vector",
    ):
        """[summary]

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

        TODO: - How do I handle instances where a layer exists only on some guides?
              - nan vs 0 when no data on var layers
        """

        if len(self.guides) == 0:
            raise ValueError(
                "No gRNAs to rank. Try running `find_guides()` method first."
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
        if layer_names:
            for layer in guide_layer_names:
                if layer not in self._layers.keys() and layer not in guide_layer_names:
                    print(self._layers.keys(), guide_layer_names)
                    raise ValueError(
                        f"Layer {layer} is not added to `Locus` or `Guide`."
                    )
        else:
            layer_names = list(guide_layer_names) + list(self._layers.keys())

        x_matrix = self._prepare_alt_matrix(rank_layer_names=layer_names)
        rank_scores = mcdm.rank(
            x_matrix,
            n_method=norm_method,
            w_vector=weight_vector,
            is_benefit_x=layer_is_benefit,
            s_method=ranking_method,
            alt_names=[g.id for g in self.guides],
        )

        rank_asc = np.argsort([-r[1] for r in rank_scores])

        for i, (g_id, rank_score) in enumerate(rank_scores):
            self.guide(g_id).rank_score = rank_score
            self.guide(g_id).rank = rank_asc[i] + 1

        return rank_scores

    def guides_to_dataframe(self):
        """Returns gRNAs in Pandas dataframe."""
        return _guides_to_dataframe(self.guides)

    def guides_to_csv(self, filename):
        """Save gRNAs in CSV file."""
        if filename:
            return _guides_to_csv(self.guides, filename)
        else:
            raise ValueError("Filename required to save CSV with gRNAs.")

    def guides_to_bed(self, filename):
        """Save gRNAs in BED file."""
        if filename:
            return _guides_to_bed(self.guides, filename)
        else:
            raise ValueError("Filename required to save BED with gRNAs.")


"""
Locus creation ------------------
"""


def _prepare_annotation(annotation_file_abspath, as_df=True):
    if annotation_file_abspath.suffix in [".gff3", ".gff"]:
        ann_db = pyranges.read_gff3(str(annotation_file_abspath), as_df=as_df)
        if as_df:
            ann_db = ann_db.rename(columns={"Name": "Exon"})  # type: ignore
    elif annotation_file_abspath.suffix in [".gtf"]:
        ann_db = pyranges.read_gtf(str(annotation_file_abspath), as_df=as_df)
        if as_df:
            ann_db = ann_db.rename(columns={"gene_id": "ID", "exon_number": "Exon"})  # type: ignore
    else:
        raise ValueError(
            "Annotation file not recognised. Annotation file needs to be GFF3 or GTF formnat."
        )
    return ann_db


def locus_from_coordinates(genome, chromosome, start, end):
    """[summary]

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
        locus_annotation = ann_db.intersect(
            pyranges.PyRanges(chromosomes=[chromosome], starts=[start], ends=[end])
        ).df

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


def locus_from_sequence(sequence, sequence_name=None):
    """[summary]

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


def locus_from_gene(genome, gene_name):
    """[summary]

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
            gene_annotation = ann_db.query('ID == @gene_name & Feature == "gene"')
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
