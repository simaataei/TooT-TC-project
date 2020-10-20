from Bio import SeqIO
import numpy as np
import os
from sklearn import metrics
import matplotlib.pyplot as plt
from math import log
import csv


def psiblast(evalue, cv, data):

    random_data = list(np.random.permutation(data))
    for i in range(1, cv+1):
        cmd1 = f"makeblastdb -in ../Results/BLAST/evalue-{evalue}/subfamily_trainsetE{evalue}Fold{i}.fasta -dbtype prot -out ../Results/BLAST/evalue-{evalue}/subfamily_trainsetE{evalue}Fold{i}DB"
        os.system(cmd1)
        cmd2 = f"psiblast -query ../Results/BLAST/evalue-{evalue}/subfamily_testsetE{evalue}Fold{i}.fasta -db ../Results/BLAST/evalue-{evalue}/subfamily_trainsetE{evalue}Fold{i}DB -num_iterations 5 -outfmt 6 -out ../Results/PsiBLAST/evalue-{evalue}/subfamily_PsiBLASTE{evalue}Fold{i}.list -evalue {evalue}"
        os.system(cmd2)

    s=2

def evaluation(neighbours,evalue, cv):

    accuracy_list =[]
    precision_list= []
    f1_list = []
    recall_list = []

    total_unclassified =[]
    for i in range(1, cv + 1):
        unclassified = []
        true_family ={}
        pred_family = {}

        y_true = []
        y_pred = []
        test_sequences = list(
            SeqIO.parse('../Results/BLAST/evalue-'+str(evalue)+'/subfamily_testsetE'+str(evalue)+'Fold'+str(i)+'.fasta',
                        'fasta'))
        for item in test_sequences:
            family = '.'.join(item.id.split('|')[3].split('.')[0:4])
            true_family[item.id] = family
        with open(f'../Results/PsiBLAST/evalue-{evalue}/subfamily_PsiBLASTE{evalue}Fold{i}.list') as file:
            results = file.readlines()
            for result in results:
                hit = result.split('\t')
                if hit[0] not in pred_family.keys() and result != '\n' and result != 'Search has CONVERGED!\n':
                    hit_family = '.'.join(hit[1].split('|')[3].split('.')[0:4])
                    pred_family[hit[0]] = hit_family

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
        total_unclassified.append(len(unclassified)/len(test_sequences))

    with open("../Results/PsiBLAST/Unclassified/unclassified_subfamily_evalue_" + str(evalue) + ".fasta",
                  "w") as handle:
        for rec in test_sequences:
            if rec.id in unclassified:
                SeqIO.write(rec, handle, "fasta")
    avg_accuracy = sum(accuracy_list)/cv
    avg_precision = sum(precision_list) / cv
    avg_recall = sum(recall_list) / cv
    avg_f1 = sum(f1_list)/cv
    avg_unclassified = sum(total_unclassified)/5
    return avg_accuracy,avg_precision,avg_recall,avg_f1,avg_unclassified


def plotting(log_evalue, accuracy_list, f1_list, precision_list, recall_list):

    plt.ylim(0.8, 1.0)
    plt.rcParams['xtick.color'] = '#333F4B'
    plt.bar(log_evalues, accuracy_list, width=0.25)
    plt.plot(log_evalues, accuracy_list)
    plt.savefig('../Results/PsiBlast_subfamily_accuracy.png', dpi=300)

    plt.show()

    plt.bar(log_evalues, f1_list)

    plt.savefig('../Results/PsiBlast_subfamily_f1.png', dpi=300)
    plt.show()

    plt.bar(log_evalues, precision_list)
    plt.savefig('../Results/PsiBlast_subfamily_result.png', dpi=300)
    plt.show()

    plt.bar(log_evalues, recall_list)
    plt.savefig('../Results/PsiBlast_subfamily_recall.png', dpi=300)
    plt.show()

    plt.bar(log_evalues, unclassified_list)
    plt.savefig('../Results/PsiBlast_subfamily_unclassified.png', dpi=300)
    plt.show()


data = list(SeqIO.parse("../Dataset/Selected30_subfamilies_tcdb2.fasta", "fasta"))

evalues= [10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001, 0.0000000001,0.00000000001,0.000000000001,
      1e-13,1e-14, 1e-15, 1e-16, 1e-17, 1e-18,1e-19,1e-20, 1e-21, 1e-22, 1e-23, 1e-24, 1e-25, 1e-26, 1e-27,1e-28, 1e-29, 1e-30,
            1e-31, 1e-32]

log_evalues = [-1*log(y,10) for y in evalues]
#evalues = np.asarray(evalues)
#log_evalues = np.log10(evalues)


#for e in evalues:
 #   psiblast(e,5,data)

#evalues = [0.0001]
accuracy_list = []
precision_list = []
recall_list = []
f1_list = []
total_unclassified = []
unclassified_list = []
for e in evalues:
    acc, precision, recall, f1, unclassified = evaluation(1, e, 5)
    accuracy_list.append(acc)
    precision_list.append(precision)
    recall_list.append(recall)
    f1_list.append(f1)
    unclassified_list.append(unclassified)
#log_evalues = [log(y,10) for y in evalues]
rows = zip(evalues,accuracy_list, precision_list, recall_list,f1_list,unclassified_list)

with open('../Results/PsiBLAST/PsiBlast_result-subfamily.csv', "w", newline='') as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)



plt.xticks(log_evalues, rotation=90)
plt.plot(log_evalues, accuracy_list)
plt.plot(log_evalues, f1_list)
plt.plot(log_evalues, precision_list)
plt.plot(log_evalues, recall_list)
# plt.xlim(-20, 2)
plt.legend(['Accuracy', 'F1-Score', 'Precision', 'Recall'])
plt.xlabel('-Log(E-value)')
plt.ylabel('Percentage')
plt.grid()
plt.savefig('../Results/PsiBlast_subfamily_result.png', dpi=300)
plt.show()
m =max(f1_list)