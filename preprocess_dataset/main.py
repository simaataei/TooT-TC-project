from Bio import SeqIO

original_data = list(SeqIO.parse("../Dataset/tcdb17feb2020.fasta", "fasta"))
import numpy as np

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


# Making a family map for data extraction
'''
family_map = {} #key:family number, value: id
i = 0
for item in original_data:
     famil = item.id.split('|')[3].split('.')[0:3]
     family = '.'.join(famil)

     if family not in family_map:
         family_map[family] = i
         i += 1


'''


'''
#extracting 30 random seqs for each family
data = list(SeqIO.parse("../Dataset/normal_len_seq.fasta", "fasta"))
selected_list_seqReq = []
removed_list_seqReq = []
for fm in family_map:
     list_temp = []
     for ro in data:
         if '.'.join(ro.id.split('|')[-1].split('.')[0:3]) == fm:
             list_temp.append(ro)
     if len(list_temp) > 29:
         new_random_list = list(np.random.permutation(list_temp))[0:30]
         selected_list_seqReq.append(new_random_list)
     else:
         removed_list_seqReq.append(fm)

print(removed_list_seqReq)
print(len(removed_list_seqReq))
print(selected_list_seqReq)

with open("../Dataset/Selected30_families_tcdb.fasta", "w") as handle:
     for rec in selected_list_seqReq:
         SeqIO.write(rec, handle, "fasta")

'''

# making a subfamily map for data extraction
subfamily_map = {} #key:subfamily number, value: id
i = 0
for item in original_data:
     subfamil = item.id.split('|')[-1].split('.')[0:4]
     subfamily = '.'.join(subfamil)

     if subfamily not in subfamily_map:
         subfamily_map[subfamily] = i
         i += 1
print(subfamily_map)

#extracting 30 random seqs for each subfamily
data = list(SeqIO.parse("../Dataset/normal_len_seq.fasta", "fasta"))
selected_list_seqReq = []
removed_list_seqReq = []
for sfm in subfamily_map:
     list_temp = []
     for ro in data:
         if '.'.join(ro.id.split('|')[-1].split('.')[0:4]) == sfm:
             list_temp.append(ro)
     if len(list_temp) > 29:
         new_random_list = list(np.random.permutation(list_temp))[0:30]
         selected_list_seqReq.append(new_random_list)
     else:
         removed_list_seqReq.append(sfm)

print(removed_list_seqReq)
print(len(removed_list_seqReq))
print(len(selected_list_seqReq))

with open("../Dataset/Selected30_subfamilies_tcdb.fasta", "w") as handle:
     for rec in selected_list_seqReq:
         SeqIO.write(rec, handle, "fasta")
data = list(SeqIO.parse("../Dataset/Selected30_subfamilies_tcdb.fasta", "fasta"))
print(len(data))

