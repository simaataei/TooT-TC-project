# %%

from Bio import SeqIO
import random

# extract a balanced dataset from family 3.a.7 and all other families
def extract_balanced_data():
    data_id = []
    data = list(SeqIO.parse("../Dataset/tcdb17feb2020.fasta", "fasta"))
    for item in data:
        data_id.append(item.id)

    not_family_3A7 = [item for item in data_id if not item.startswith('3.A.7')]
    selected_not_id = random.sample(not_family_3A7, 30)

    family_3A7 = list(SeqIO.parse("../Dataset/family_3A7/Random_family_3A7.fasta", "fasta"))

    # train_not_id = random.sample(selected_not_id, 6)
    # train_3A7_id =

    new_dataset = []
    for item in data:
        if item.id in selected_not_id:
            new_dataset.append(item)

    #train_set = random.sample(new_dataset, 24) + random.sample(family_3A7, 24)
    new_dataset = new_dataset + family_3A7


    with open("../Dataset/family_3A7/Balanced_30family3A7_30other.fasta", "w") as f:
        for rec in new_dataset:
            SeqIO.write(rec, f, "fasta")


def extract_test_train():
    dataset = list(SeqIO.parse("../Dataset/family_3A7/Balanced_30family3A7_30other.fasta", "fasta"))
    data_id = []
    for item in dataset:
        data_id.append(item.id)

    not_family_3A7 = [item for item in data_id if not item.split("|")[3].startswith('3.A.7')]
    family_3A7 = [item for item in data_id if item.split("|")[3].startswith('3.A.7')]

    train_not_family = random.sample(not_family_3A7, 24)
    train_family = random.sample(family_3A7, 24)
    test_family = [item for item in family_3A7 if item not in train_family]
    test_not_family = [item for item in not_family_3A7 if item not in train_not_family ]

    train_id = train_family + train_not_family
    test_id =test_family + test_not_family


    train =[]
    test =[]
    for item in dataset:
        if item.id in train_id:
            train.append(item)
    for item in dataset:
        if item.id in test_id:
            test.append(item)

    with open("../Dataset/family_3A7/Balanced30_trainingset.fasta", "w") as f:
        for rec in train:
            SeqIO.write(rec, f, "fasta")
    with open("../Dataset/family_3A7/Balanced30_testset.fasta", "w") as f:
        for rec in test:
            SeqIO.write(rec, f, "fasta")


#extract_balanced_data()
#extract_test_train()


