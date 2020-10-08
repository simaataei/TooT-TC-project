from Bio import SeqIO
import numpy as np
import os
from sklearn import metrics

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

def evaluation(neighbours,evalue, cv):

    accuracy_list =[]
    precision_list= []
    f1_list = []
    recall_list = []
    for i in range(1, cv + 1):
        true_family ={}
        pred_family = {}
        unclassified =[]
        y_true = []
        y_pred = []
        test_sequences = list(
            SeqIO.parse('../Results/BLAST/evalue-'+str(evalue)+'/testsetE'+str(evalue)+'Fold'+str(i)+'.fasta',
                        'fasta'))
        for item in test_sequences:
            family = '.'.join(item.id.split('|')[3].split('.')[0:3])
            true_family[item.id] = family

        with open(f'../Results/BLAST/evalue-{evalue}/BLASTE{evalue}Fold{i}.list') as file:
            results = file.readlines()
            for result in results:
                hit = result.split('\t')
                if hit[0] not in pred_family.keys():
                    hit_family = '.'.join(hit[1].split('|')[3].split('.')[0:3])
                    pred_family[hit[0]] = hit_family
                    a=1
        for key in true_family.keys():
            if key not in pred_family.keys():
                unclassified.append(key)
            else:
                y_true.append(true_family[key])
                y_pred.append(pred_family[key])
        accuracy_list.append(metrics.accuracy_score(y_true, y_pred))
        precision_list.append(metrics.precision_score(y_true, y_pred, average="macro"))
        recall_list.append(metrics.recall_score(y_true, y_pred, average="macro"))
        f1_list.append(metrics.f1_score(y_true, y_pred, average="macro"))
    avg_accuracy = sum(accuracy_list)/cv
    avg_precision = sum(precision_list) / cv
    avg_recall = sum(recall_list) / cv
    avg_f1 = sum(f1_list)/cv
    return avg_accuracy,avg_precision,avg_recall,avg_f1








data = list(SeqIO.parse("../Dataset/Selected30_families_tcdb2.fasta", "fasta"))

#11, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001,
evalues= [0.00000000001,0.000000000001,0.000000000001]
for e in evalues:
    blast(e,5,data)


accuracy_list = []
precision_list = []
recall_list = []
f1_list = []
for e in evalues:
    acc,precision,recall,f1 = evaluation(1,e,5)
    accuracy_list.append(acc)
    precision_list.append(precision)
    recall_list.append(recall)
    f1_list.append(f1)
a=1