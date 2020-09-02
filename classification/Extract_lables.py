from Bio import SeqIO
def extract_families(data):
    for item in data:
        famil = item.id.split('|')[3].split('.')[0:3]
        family = '.'.join(famil)


def extract_subfamilies(data):
    for item in data:
        subfamil = item.id.split('|')[3].split('.')[0:4]
        subamily = '.'.join(subfamil)
        print(subamily)

data = list(SeqIO.parse("../Dataset/normal_len_seq2.fasta", "fasta"))
extract_subfamilies(data)