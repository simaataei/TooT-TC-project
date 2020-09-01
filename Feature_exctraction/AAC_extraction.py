import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def convert_to_aac(dataset):
    i = 0
    j = 0
    X_aac = np.zeros((len(dataset), 20))
    for seq in dataset:
        analysed_seq = ProteinAnalysis(str(seq.seq))
        for val in analysed_seq.count_amino_acids().values():
            X_aac[i][j] = val / len(seq)
            j += 1
        i += 1
        j = 0
    return X_aac

data = list(SeqIO.parse("../Dataset/Selected30_families_tcdb2.fasta", "fasta"))
aac = convert_to_aac(data)
