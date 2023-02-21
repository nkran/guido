import os
from distutils import dir_util

import pytest
from pyfaidx import Fasta
from pytest import fixture

from guido.genome import load_genome_from_file
from guido.locus import locus_from_coordinates, locus_from_gene
from guido.off_targets import calculate_ot_sum_score

# prepare files - fixtures
# move genome to separate function to load when needed
# download to temp from vectorbase


@fixture
def datadir(tmpdir, request):
    """Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely."""
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, str(tmpdir))

    return tmpdir


def load_genome():
    return load_genome_from_file(
        guido_file="/Users/nkranjc/imperial/ref/new/AgamP4.guido"
    )


@pytest.mark.parametrize(
    "gene_name, expected_locus",
    [
        (
            "AGAP005958",
            ("AgamP4_2L", 24049834, 24051517, "TCCAGTCCAAGGTAGTCAGTATCA"),
        ),
        (
            "AGAP007280",
            ("AgamP4_2L", 44989066, 44998118, "CGGTCGTGGCTGAACACACAGTCCA"),
        ),
    ],
)
def test_locus_from_gene(gene_name, expected_locus):
    genome = load_genome()
    loc = locus_from_gene(genome, gene_name)

    assert loc.chromosome == expected_locus[0]
    assert loc.start == expected_locus[1]
    assert loc.end == expected_locus[2]
    assert loc.sequence.seq[0 : len(expected_locus[3])] == expected_locus[3]


@pytest.mark.parametrize(
    "chromosome, start, end, expected_locus",
    [
        (
            "AgamP4_2R",
            48714541,
            48714666,
            ("AgamP4_2R", 48714541, 48714666, "TGGTGCGGAAAGTTTAT", "TGTGTTAAACATAAATG"),
        ),
    ],
)
def test_locus_from_coordinates(chromosome, start, end, expected_locus):
    genome = load_genome()
    loc = locus_from_coordinates(genome, chromosome, start, end)

    assert loc.chromosome == expected_locus[0]
    assert loc.start == expected_locus[1]
    assert loc.end == expected_locus[2]
    assert loc.sequence.seq[0 : len(expected_locus[3])] == expected_locus[3]
    assert loc.sequence.seq[-len(expected_locus[4]) :] == expected_locus[4]

    assert "AGAP004050" in loc.annotation["ID"].values


@pytest.mark.parametrize(
    "chromosome, start, end",
    [
        ("AgamP4_2R", 48714541, 48714666),
    ],
)
def test_find_guides(chromosome, start, end):
    genome = load_genome()
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


@pytest.mark.parametrize(
    "chromosome, start, end",
    [
        ("AgamP4_2R", 48714541, 48714666),
    ],
)
def test_find_offtargets(chromosome, start, end):
    genome = load_genome()
    loc = locus_from_coordinates(genome, chromosome, start, end)

    guides_all = loc.find_guides()

    # check if we have the correct number of guides
    assert len(guides_all) == 8

    offtargets = loc.find_off_targets()
    assert len(offtargets.keys()) == 8

    loc.guide(1).off_targets = offtargets[1]
    loc.guide(0).off_targets = offtargets[0]

    assert calculate_ot_sum_score(loc.guide(0).off_targets) == 8
    assert calculate_ot_sum_score(loc.guide(1).off_targets) == 7

    assert loc.guide(0).off_targets_string == "0|0|0|0|1|5 (6)"
    assert loc.guide(1).off_targets_string == "0|0|0|0|0|7 (7)"


@pytest.mark.parametrize(
    "chromosome, start, end",
    [
        ("AgamP4_2R", 48714541, 48714666),
    ],
)
def test_azimuth_score(chromosome, start, end):
    genome = load_genome()
    loc = locus_from_coordinates(genome, chromosome, start, end)

    # find guides
    loc.find_guides()

    # add azimuth scores
    scores = loc.add_azimuth_score()

    assert len(scores) == len(loc.guides)
    assert loc.guide(7).layer("azimuth_score") == 0.0
    assert len(loc.guide(7).guide_long_seq) != 30

    # test if method works the same when called on the locus or on the guide
    # test if compare_model works the same when called on the locus or on the sequence independently in the test


@pytest.mark.parametrize(
    "chromosome, start, end",
    [
        ("AgamP4_2R", 48714541, 48714666),
    ],
)
def test_cfd_score(chromosome, start, end):

    genome = load_genome()
    loc = locus_from_coordinates(genome, chromosome, start, end)

    # find guides
    loc.find_guides()
    loc.find_off_targets()
    loc.add_azimuth_score()
    loc.simulate_end_joining()

    guide = loc.guide(0)

    assert "ot_cfd_score_mean" in guide.layers and "ot_cfd_score_max" in guide.layers

    # test the expected values from cfd score


@pytest.mark.parametrize(
    "chromosome, start, end",
    [
        ("AgamP4_2R", 48714541, 48714666),
    ],
)
def test_ranking(chromosome, start, end):

    genome = load_genome()
    loc = locus_from_coordinates(genome, chromosome, start, end)

    # find guides
    loc.find_guides()
    loc.find_off_targets()
    loc.add_azimuth_score()
    loc.simulate_end_joining()
    loc.rank_guides()

    assert loc.guide(0).rank == 3
    assert loc.guide(0).rank_score == pytest.approx(0.4596, rel=1e-3)

    assert loc.guide(1).rank_score == pytest.approx(0.6323, rel=1e-3)
    assert loc.guide(1).rank == 2

    assert loc.guide(0).rank > loc.guide(1).rank
    assert loc.guide(1).rank_score > loc.guide(0).rank_score

    # test the expected values from cfd score
