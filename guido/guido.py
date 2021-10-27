from guido.locus import Locus

"""
- check if N in flanks / guide sequence
- remove partial gRNAs                                  done
- refactor load_from_coordinates() return
- multithreading
- create Genome from FASTA file
- refactor guide.find_off_targets()                     done
- check if bowtie is installed
"""


def main():

    seq = "GCCGACCCATTCTGCTGCCCTTCTGTACCGTGGTGCGGCTCTCTCGCTCCACTCCTTAAACACTAGTTTGAACTTATCGGCATCAGTTGCGCACGCGGCTTGATTTAAAATAGCACAGAACTATTGAATTCGTTTCACCAAacacacatacacacacccacatacaAAGATACGGACAGTTACAGTGGTGCGGAAAGTTTATCATCCACTCTGACGGGTGGTATTGCGCAACTCCACGCCATCAAACATGTTCAGATTATGCAATCGTGAGTATTCGTTGACCACCGCTTGACCTGTGTTAAACATAAATGAATGGAAAGGTAAGGCTTTGAAGGTCACTGCTGCTGGCTGACGGAATTCACAAtttggtttttgatgttttggttttttttttGTATCGAATTTTGAAGTCAGTGAACGTGGCATAACACCATATGCCGCTACCTTCAAGATGCAGATACTCCTAACTTCTCGTGTCTGAGCTAGCTAA"
    locus = Locus(sequence=seq)

    locus.find_guides("NGG")
