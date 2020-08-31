from Bio import SeqIO
original_data = list(SeqIO.parse("../Dataset/tcdb17feb2020.fasta", "fasta"))


# remove unacceptable amino acids
'''
unacceptable_amino=['B', 'Z', 'X', 'J', 'O', 'U']
delete_these=[]
for od in original_data:
    for amino_acid in od.seq:
        if amino_acid in unacceptable_amino:
           # original_data.remove(od)
           if od.id not in delete_these:
               delete_these.append(od.id)


acceptable_records =[]
for item in original_data:
    if item.id not in delete_these:
        acceptable_records.append(item)
print(len(original_data))
print(len(acceptable_records))
with open("../Dataset/acceptable_aminoacid_seq.fasta", "w") as handle:
    for rec in acceptable_records:
        SeqIO.write(rec, handle, "fasta")


'''
'''
#Sequence Length refinment
data = list(SeqIO.parse("../Dataset/acceptable_aminoacid_seq.fasta", "fasta"))
Max_len = 1000
Min_len = 50
normal_len_dataset=[]
for s in data :
    if Min_len<len(s.seq)<Max_len:
        normal_len_dataset.append(s)

with open("../Dataset/normal_len_seq.fasta", "w") as handle:
    for rec in normal_len_dataset:
        SeqIO.write(rec, handle, "fasta")
'''
