from Bio import SeqIO
import numpy as np
import os


def blast(evalue, cv, data):
    test_data = []
    train_data = []
    test_list = []
    random_data = list(np.random.permutation(data))
    for i in range(1, cv+1):
        test_data = random_data[int((i-1)*1/(cv)*len(random_data)):int(i*1/(cv)*len(random_data))]
        for t in test_data:
            test_list.append(t.id)
        for seq in random_data:
            if seq.id not in test_list:
                train_data.append(seq)

        with open("../Results/BLAST/evalue-"+str(evalue)+"/testsetE"+str(evalue)+"Fold"+str(i)+".fasta", "w") as handle:
            for rec in test_data:
                SeqIO.write(rec, handle, "fasta")
        with open("../Results/BLAST/evalue-"+str(evalue)+"/trainsetE"+str(evalue)+"Fold"+str(i)+".fasta", "w") as handle:
            for rec in train_data:
                SeqIO.write(rec, handle, "fasta")
        cmd1 = f"makeblastdb -in ../Results/BLAST/evalue-{evalue}/trainsetE{evalue}Fold{i}.fasta -dbtype prot -out ../Results/BLAST/evalue-{evalue}/trainsetE{evalue}Fold{i}DB"
        os.system(cmd1)
        cmd2 = f"blastp -query ../Results/BLAST/evalue-{evalue}//testsetE{evalue}Fold{i}.fasta -db ../Results/BLAST/evalue-{evalue}/trainsetE{evalue}Fold{i}DB -outfmt 6 -out ../Results/BLAST/evalue-{evalue}/BLASTE{evalue}Fold{i}.list -evalue {evalue}"
        os.system(cmd2)
        test_data = []
        train_data = []
        test_list = []
    s=2

def evaluation(hits):


data = list(SeqIO.parse("../Dataset/Selected30_families_tcdb2.fasta", "fasta"))
blast(10,5,data)