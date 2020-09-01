import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def convert_to_paac(dataset):
    # Separate sequences
    sequences = [str(ro.seq) for ro in dataset]
    # Map amino acids
    amino_map = {amino: i for i, amino in enumerate(list(set(''.join(sequences))))}
    list_pair_dict = {}
    for first in amino_map.keys():
        for second in amino_map.keys():
            if first + second not in list_pair_dict:
                list_pair_dict[first + second] = len(list_pair_dict)

    row = 0
    X_paac = np.zeros((len(dataset), len(list_pair_dict)))
    for seq in dataset:
        pair_amino_dict = {}
        i = 0
        pair_str = ''
        for char in str(seq.seq):
            pair_str += char
            i += 1
            if i == 2:
                if pair_str not in pair_amino_dict:
                    pair_amino_dict[pair_str] = 1
                    i = 1
                    pair_str = pair_str[1]
                else:
                    pair_amino_dict[pair_str] += 1
                    i = 1
                    pair_str = pair_str[1]
        for pair in pair_amino_dict:
            X_paac[row][list_pair_dict[pair]] = pair_amino_dict[pair] / len(seq)
        row += 1
    return X_paac

data = list(SeqIO.parse("../Dataset/Selected30_families_tcdb2.fasta", "fasta"))
paac = convert_to_paac(data)
print(paac)