from Bio import SeqIO
def extract_families(data):
    family = []
    for item in data:
        famil = item.id.split('|')[3].split('.')[0:3]
        family.append('.'.join(famil))
    return family


def extract_subfamilies(data):
    subfamily = []

    for item in data:
        subfamil = item.id.split('|')[3].split('.')[0:4]
        subfamily.append('.'.join(subfamil))
    return subfamily
#data = list(SeqIO.parse("../Dataset/normal_len_seq2.fasta", "fasta"))
#print(extract_families(data))