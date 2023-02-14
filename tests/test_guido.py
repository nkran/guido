import pytest
from pyfaidx import Fasta

from guido.genome import load_genome_from_file
from guido.locus import locus_from_coordinates, locus_from_gene

# prepare files
# download to temp from vectorbase

genome = load_genome_from_file(
    guido_file="/Users/nkranjc/imperial/ref/new/AgamP4.guido"
)


@pytest.mark.parametrize(
    "genome, gene_name, expected_locus",
    [
        (
            genome,
            "AGAP005958",
            ("AgamP4_2L", 24049834, 24051517, "TCCAGTCCAAGGTAGTCAGTATCA"),
        ),
        (
            genome,
            "AGAP007280",
            ("AgamP4_2L", 44989066, 44998118, "CGGTCGTGGCTGAACACACAGTCCA"),
        ),
    ],
)
def test_locus_from_gene(genome, gene_name, expected_locus):
    loc = locus_from_gene(genome, gene_name)

    assert loc.chromosome == expected_locus[0]
    assert loc.start == expected_locus[1]
    assert loc.end == expected_locus[2]
    assert loc.sequence.seq[0 : len(expected_locus[3])] == expected_locus[3]


@pytest.mark.parametrize(
    "genome, chromosome, start, end, expected_locus",
    [
        (
            genome,
            "AgamP4_2R",
            48714541,
            48714666,
            ("AgamP4_2R", 48714541, 48714666, "TGGTGCGGAAAGTTTAT", "TGTGTTAAACATAAATG"),
        ),
    ],
)
def test_locus_from_coordinates(genome, chromosome, start, end, expected_locus):
    loc = locus_from_coordinates(genome, chromosome, start, end)

    assert loc.chromosome == expected_locus[0]
    assert loc.start == expected_locus[1]
    assert loc.end == expected_locus[2]
    assert loc.sequence.seq[0 : len(expected_locus[3])] == expected_locus[3]
    assert loc.sequence.seq[-len(expected_locus[4]) :] == expected_locus[4]

    assert "AGAP004050" in loc.annotation["ID"].values


@pytest.mark.parametrize(
    "genome, chromosome, start, end",
    [
        (genome, "AgamP4_2R", 48714541, 48714666),
    ],
)
def test_find_guides(genome, chromosome, start, end):
    loc = locus_from_coordinates(genome, chromosome, start, end)

    guides_all = loc.find_guides()
    guides_exon = loc.find_guides(selected_features="exon")

    # check if we have the correct number of guides
    assert len(guides_all) == 8
    assert len(guides_exon) == 6

    # check if guide sequences are correct
    for g in guides_all:
        rc = False
        if g.guide_strand == "-":
            rc = True

        seq = Fasta(str(genome.genome_file_abspath), one_based_attributes=True).get_seq(
            g.guide_chrom, g.guide_start, g.guide_end, rc
        )

        assert g.guide_seq == seq.seq
