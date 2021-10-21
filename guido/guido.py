import os
import re
import sys
import pickle
import itertools
import time
import tempfile
import subprocess

# import pyranges
import pandas as pd
import multiprocessing as mp

# from pyfaidx import Fasta

from typing import List, Dict, Tuple

from guido.guides import Guide
from guido.mmej import generate_mmej_patterns
from guido.off_targets import run_bowtie

"""
- check if N in flanks / guide sequence
- remove partial gRNAs                                  done
- refactor load_from_coordinates() return
- multithreading
- create Genome from FASTA file
- refactor guide.find_off_targets()
- check if bowtie is installed
"""


def main():

    seq = "GCCGACCCATTCTGCTGCCCTTCTGTACCGTGGTGCGGCTCTCTCGCTCCACTCCTTAAACACTAGTTTGAACTTATCGGCATCAGTTGCGCACGCGGCTTGATTTAAAATAGCACAGAACTATTGAATTCGTTTCACCAAacacacatacacacacccacatacaAAGATACGGACAGTTACAGTGGTGCGGAAAGTTTATCATCCACTCTGACGGGTGGTATTGCGCAACTCCACGCCATCAAACATGTTCAGATTATGCAATCGTGAGTATTCGTTGACCACCGCTTGACCTGTGTTAAACATAAATGAATGGAAAGGTAAGGCTTTGAAGGTCACTGCTGCTGGCTGACGGAATTCACAAtttggtttttgatgttttggttttttttttGTATCGAATTTTGAAGTCAGTGAACGTGGCATAACACCATATGCCGCTACCTTCAAGATGCAGATACTCCTAACTTCTCGTGTCTGAGCTAGCTAA"
    locus = Locus(sequence=seq)

    locus.find_guides("NGG")
