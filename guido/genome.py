import pickle
import subprocess
from pathlib import Path

from pyfaidx import Faidx


def check_file(filename, supported_ext):
    """
    _summary_

    Parameters
    ----------
    filename : _type_
        _description_
    supported_ext : _type_
        _description_

    Returns
    -------
    bool
        _description_
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
    ):
        """
        _summary_

        Parameters
        ----------
        genome_name : _type_
            _description_
        genome_file_abspath : _type_, optional
            _description_, by default None
        annotation_file_abspath : _type_, optional
            _description_, by default None

        Raises
        ------
        ValueError
            _description_
        ValueError
            _description_
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
        """
        _summary_

        Returns
        -------
        _type_
            _description_
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
        """
        _summary_

        Parameters
        ----------
        n_threads : int, optional
            _description_, by default 1
        save_pickle : bool, optional
            _description_, by default True
        bowtie_ignore : bool, optional
            _description_, by default False
        bowtie_path : str, optional
            _description_, by default ""
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


def load_genome_from_file(guido_file):
    """
    _summary_

    Parameters
    ----------
    guido_file : _type_
        _description_

    Returns
    -------
    _type_
        _description_

    Raises
    ------
    ValueError
        _description_
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
