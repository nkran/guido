import pickle
import subprocess
from pathlib import Path

from pyfaidx import Faidx, Fasta


def check_file(filename, supported_ext):
    """Check if file exists and has the right extension.

    Parameters
    ----------
    filename : str
        Path to the file.
    supported_ext : str
        Supported file extensions.

    Returns
    -------
    bool
        True if file exists and has the right extension, False otherwise.
    """
    filename = Path(filename)
    if filename.exists() and filename.suffix in supported_ext:
        return True
    else:
        return False


class Genome:
    def __init__(
        self,
        genome_name,
        genome_file_abspath=None,
        annotation_file_abspath=None,
        bowtie_index_abspath=None,
    ):
        """Genome class.

        Parameters
        ----------
        genome_name : str
            Name of the genome.
        genome_file_abspath : str, optional
            Path to the genome Fasta file, by default None
        annotation_file_abspath : str, optional
            Path to the genome annotation file (GTF or GFF3), by default None
        bowtie_index_abspath : str, optional
            Path to the bowtie index files, by default None
        """
        self.genome_name = genome_name
        self._bowtie_ignore = False
        self.bowtie_index = None

        if genome_file_abspath:
            if check_file(genome_file_abspath, [".fa", ".fasta", ".fna"]):
                self.genome_file_abspath = Path(genome_file_abspath).resolve()
            else:
                raise ValueError(
                    f"Genome file {genome_file_abspath} does noth exist or is not in the right format."
                )

        if bowtie_index_abspath:
            self.bowtie_index = Path(bowtie_index_abspath)
            self._bowtie_ignore = True

        if annotation_file_abspath:
            if check_file(annotation_file_abspath, [".gtf", ".gff3"]):
                self.annotation_file_abspath = Path(annotation_file_abspath).resolve()
            else:
                raise ValueError(
                    f"Annotation file file {genome_file_abspath} does noth exist or is not in the right format."
                )
        else:
            self.annotation_file_abspath = None

    @property
    def is_built(self):
        """Check if genome index is built.

        Returns
        -------
        Bool
            True if genome index is built, False otherwise.
        """
        ready = []
        if self.genome_file_abspath and self.fai_file:
            ready.append(True)
        else:
            return False

        if self.annotation_file_abspath and not self.tbi_file:
            ready.append(False)
        else:
            ready.append(True)

        bowtie_index_files = [
            i
            for i in self.genome_file_abspath.parent.glob("*ebwt")
            if f"{self.genome_name}.rev." in str(i) and f"{self.genome_name}." in str(i)
        ]
        if not self._bowtie_ignore and any(bowtie_index_files):
            ready.append(True)
        elif not self._bowtie_ignore and not any(bowtie_index_files):
            ready.append(False)

        return all(ready)

    def build(
        self,
        n_threads=1,
        save_pickle=True,
        bowtie_ignore=False,
        bowtie_path="",
    ):
        """Build genome index.

        This method creates index files for the genome and annotation files and saves them
        in a `.guido` file in the same directory as the genome Fasta file. This file can be
        later loaded using :meth:`guido.genome.load_genome_from_file` without having to
        re-build the index.

        Parameters
        ----------
        n_threads : int, optional
            Number of threads, by default 1
        save_pickle : bool, optional
            Pickle the dictionary into `.guido` file, by default True
        bowtie_ignore : bool, optional
            Ignore building bowtie index. Use if you already have bowtie index built, by default False
        bowtie_path : str, optional
            Path to bowtie binary if it's not in the path, by default ""
        """

        self._bowtie_ignore = bowtie_ignore
        if self.annotation_file_abspath:
            print("Indexing genome annotation.")
            annotation_sorted_bgz = self.annotation_file_abspath.with_suffix(
                str(self.annotation_file_abspath.suffix) + ".gz"
            )
            annotation_sort_bgz_cmd = f"sort -k1,1 -k4,4n {str(self.annotation_file_abspath)} | bgzip > {annotation_sorted_bgz}"
            self.annotation_sorted_gz_file = annotation_sorted_bgz

            _ = subprocess.run(annotation_sort_bgz_cmd, shell=True)
            tabix_cmd = f"tabix {str(self.annotation_sorted_gz_file)}"

            _ = subprocess.run(tabix_cmd, stderr=subprocess.PIPE, shell=True)
            self.tbi_file = self.annotation_sorted_gz_file.with_suffix(".tbi")

        if self.genome_file_abspath.exists():
            Faidx(str(self.genome_file_abspath))
            self.fai_file = self.genome_file_abspath.with_suffix(
                str(self.genome_file_abspath.suffix) + ".fai"
            )

            if not self._bowtie_ignore:
                bowtie_index_command = f"{bowtie_path}bowtie-build {str(self.genome_file_abspath)} {str(self.genome_file_abspath.parent / self.genome_name)} --threads {n_threads}"
                print("Building Bowtie index")
                bowtie_index_proc = subprocess.run(
                    bowtie_index_command.strip(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    shell=True,
                )
                if bowtie_index_proc.stderr:
                    print(bowtie_index_proc.stderr)
                else:
                    self.bowtie_index = Path(
                        f"{str(self.genome_file_abspath.parent / self.genome_name)}"
                    )
                    print(f"Done: {self.bowtie_index }")

            if save_pickle:
                guido_file = Path(
                    f"{str(self.genome_file_abspath.parent / self.genome_name)}.guido"
                )
                with open(guido_file, "wb") as f:
                    pickle.dump(self.__dict__, f, protocol=pickle.HIGHEST_PROTOCOL)
                print(
                    f"{self.genome_name} genome data can now be used by Guido: {str(guido_file)}"
                )

    @property
    def sequence(self):
        if self.is_built:
            return Fasta(str(self.genome_file_abspath))


def load_genome_from_file(guido_file):
    """Load a genome from a .guido file. This file is created when a genome is
    built. It contains all the information needed to use the genome. Guido
    files are saved in the same directory as the genome FASTA file. They are
    named after the genome name and have the .guido extension. Genome object
    can be created by using the :method:`build` method.

    Parameters
    ----------
    guido_file : str
        Path to the .guido file.

    Returns
    -------
    Genome
        :class:`Genome` object.
    """

    if Path(guido_file).exists():
        with open(guido_file, "rb") as f:
            guido_dict = pickle.load(f)
        G = Genome(genome_name=guido_dict["genome_name"])
        for key, val in guido_dict.items():
            setattr(G, key, val)
        return G
    else:
        raise ValueError(f"Provided path to .guido ({guido_file}) file does not exist.")
