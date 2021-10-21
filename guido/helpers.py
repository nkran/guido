def rev_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join([complement[base] if base in complement.keys() else "N" for base in seq[::-1]])
