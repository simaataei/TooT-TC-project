import numpy as np
from scipy.signal import resample
from Bio import SeqIO

sample_length = 20
def one_hot_encode(sequence, amino_map):
    rep = np.zeros((len(amino_map), len(sequence)))
    for i, amino in enumerate(sequence):
        rep[amino_map[amino], i] = 1
    return rep


def resample_sequence(sample):
    resampled_sample = np.zeros((sample.shape[0], sample_length))
    for i in range(sample.shape[0]):
        resampled_sample[i, :] = resample(sample[i, :], sample_length)
    return resampled_sample
def extract_one_hot(dataset):
    sequences = [str(ro.seq) for ro in dataset]
    amino_map = {amino: i for i, amino in enumerate(list(set(''.join(sequences))))}
    features = [one_hot_encode(seq, amino_map) for seq in sequences]
    resampled_features = []
    for sample in features:
        resampled_features.append(resample_sequence(sample).ravel())
    resampled_features = np.array(resampled_features)
    return resampled_features
#data = list(SeqIO.parse("../Dataset/Selected30_families_tcdb2.fasta", "fasta"))
#one_hot = extract_one_hot(data)
#print(one_hot)