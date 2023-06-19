# Guido
[![Documentation Status](https://readthedocs.org/projects/guido/badge/?version=latest)](https://guido.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/419776638.svg)](https://zenodo.org/badge/latestdoi/419776638)
<!-- [![PyPI version](https://badge.fury.io/py/guido.svg)](https://badge.fury.io/py/guido) -->

Guido is a Python package developed to search for gRNA targets in any reference genome or DNA sequence. It integrates MMEJ prediction and scoring, off-target search, and allows users to define their own data layers that can be used in the gRNA evaluation and ranking.

## Installation
Install `guido` and azimuth (dependency) via PyPi using pip:

```bash
$ pip install guido
$ pip install git+https://github.com/Biomatters/Azimuth

```

Please note that `guido` requires `bowtie` to be installed on your system. Please follow the instructions on the [bowtie website](http://bowtie-bio.sourceforge.net/index.shtml) to install it.

## Usage
You can use `guido` to search for gRNAs in a reference genome or DNA sequence. The following example shows how to search for gRNAs in the malaria mosquito (*Anopheles gambiae*) reference genome (AgamP4):

### Create a `Genome` instance

First we need to create a `Genome` instance. This instance will be used to search for gRNAs in the genome. The `Genome` class will link the genome sequence to the annotation file and will create a bowtie index for the genome. The `Genome` class takes the following arguments: `genome_name` (name of the genome), `genome_file_abspath` (FASTA sequence), and `annotation_file_abspath` (GTF annotation file).

```python
>>> import guido

>>> genome = guido.Genome(genome_name='AgamP4',
                      genome_file_abspath='data/AgamP4.fa',
                      annotation_file_abspath='data/AgamP4.12.gtf')
```

Build the FASTA index and bowtie index files and create a `AgamP4.guido` file in `genome_file_abspath` that contains the genome information and can be used to create a `Genome` instance next time without needing to build the genome indices again.

```python
>>> genome.build(n_threads=2)
```

The `Genome` class also has a `bowtie_index_abspath` argument that can be used to specify the directory and bowtie index name. If this argument is not specified, the bowtie index will be created when running `genome.build()`.

To load a `Genome` instance from a `AgamP4.guido` file that was created using `genome.build()`, use the `Genome.load()` method:

```python
>>> import guido
>>> genome = guido.load_genome_from_file(guido_file='data/AgamP4.guido')
```

### Create a `Locus` instance and search for gRNAs
`Locus` instances are used to search for gRNAs in a specific genomic region. The `Locus` class takes the following arguments: `genome` (a `Genome` instance), `chromosome` (chromosome name), `start` (start position), and `end` (end position). The `start` and `end` positions are 1-based and inclusive.

```python
>>> import guido

>>> genome = guido.load_genome_from_file(guido_file='data/AgamP4.guido')
>>> loc = guido.locus_from_coordinates(genome, 'AgamP4_2R', 48714541, 48714666)

>>> loc.find_guides()
    ['gRNA-1(AAGTTTATCATCCACTCTGACGG|AgamP4_2R:48714550-48714572|+|)',
     'gRNA-2(CGCAATACCACCCGTCAGAGTGG|AgamP4_2R:48714561-48714583|-|)',
     ...
     'gRNA-7(GTTTAACACAGGTCAAGCGGTGG|AgamP4_2R:48714637-48714659|-|)',
     'gRNA-8(TATGTTTAACACAGGTCAAGCGG|AgamP4_2R:48714640-48714662|-|)']
```

Alternatively, you can use the `locus_from_gene()` function to create a `Locus` instance and limit the search for gRNAs to a specific feature that is defined in the genome annotation file:

```python
>>> import guido
>>> loc = guido.locus_from_gene(genome, 'AGAP005958')

>>> loc.find_guides(feature_type='exon')
```

After running `loc.find_guides()`, the `Locus` instance will contain a list of `gRNA` instances that can be accessed using the `loc.guides` attribute.

```python
>>> loc.guides

    ['gRNA-1(AAGTTTATCATCCACTCTGACGG|AgamP4_2R:48714550-48714572|+|)',
     'gRNA-2(CGCAATACCACCCGTCAGAGTGG|AgamP4_2R:48714561-48714583|-|)',
     ...
     'gRNA-7(GTTTAACACAGGTCAAGCGGTGG|AgamP4_2R:48714637-48714659|-|'),
     'gRNA-8(TATGTTTAACACAGGTCAAGCGG|AgamP4_2R:48714640-48714662|-|)']

```

You can access a gRNA by its index or a name:

```python
>>> loc.guides[0]
    'gRNA-1(AAGTTTATCATCCACTCTGACGG|AgamP4_2R:48714550-48714572|+|)'

>>> loc.guides['gRNA-1']
    'gRNA-1(AAGTTTATCATCCACTCTGACGG|AgamP4_2R:48714550-48714572|+|)'
```

## Docs
More extensive documentation can be found on [Read the Docs](https://guido.readthedocs.io/en/latest/).

## Cite as
Nace Kranjc, & Courty Thomas. (2023). nkran/guido: v0.1.3 (v0.1.3). Zenodo. https://doi.org/10.5281/zenodo.8056051


## Developer setup
Install [poetry](https://python-poetry.org/docs/#installation):

```bash
$ pip install poetry
```

Create development environment:

```bash
$ cd guido
$ poetry install
```

Activate development environment:

```bash
$ poetry shell
```

Install pre-commit hooks:

```bash
$ pre-commit install
```

Run pre-commit checks (isort, black, blackdoc, flake8, ...) manually:

```bash
$ pre-commit run --all-files
```

Bump version, build and publish to PyPI:

```bash
$ poetry version prerelease
$ poetry build
$ poetry publish
```
