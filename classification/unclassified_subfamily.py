

from Bio import SeqIO
from classification.Extract_lables import extract_subfamilies, extract_families
from Feature_exctraction.AAC_extraction import extract_aac
from Feature_exctraction.PAAC_exraction import extract_paac
from Feature_exctraction.One_hot_encoding import extract_one_hot
import numpy as np
import sklearn.svm
import matplotlib.pyplot as plt
from sklearn import metrics

def count_subfamilies(data):
    subfamily_count ={}
    for item in data:
        subfamil = item.id.split('|')[3].split('.')[0:4]
        subfamily = '.'.join(subfamil)
        if subfamily not in subfamily_count.keys():
            subfamily_count[subfamily] = 1
        else:
            subfamily_count[subfamily] += 1
    return subfamily_count

subfamily_data = list(SeqIO.parse("../Dataset/Selected30_subfamilies_tcdb2.fasta", "fasta"))
subfamily_test_data = list(SeqIO.parse("../Results/BLAST/Unclassified/subfamily_unclassified_evalue_0.0001.fasta", "fasta"))

subfamily_unclassified= {}
subfamily_train_data = []
subfamily_test_list = []

for t in subfamily_test_data:
    subfamily_test_list.append(t.id)
for seq in subfamily_data:
    if seq.id not in subfamily_test_list:
        subfamily_train_data.append(seq)

family_count ={}
family_count = count_subfamilies(subfamily_test_data)
plt.bar(family_count.keys(), family_count.values(),color=(0.1, 0.6, 0.8, 0.8))
plt.grid()
plt.xticks(rotation= 90, fontsize= 6)
plt.xlabel("Subfamily", fontsize = 8)
plt.ylabel("Number of unclassified sequences")
plt.savefig('../Results/unclassified_subfamilies_blast_0.0001.png', dpi=300)
plt.show()
