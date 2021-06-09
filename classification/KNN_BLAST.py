from Bio import SeqIO
import numpy as np
import os
from sklearn import metrics
import matplotlib.pyplot as plt
from math import log
import csv



def evaluation_family(neighbours,evalue, cv):

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
            SeqIO.parse('../Results/BLAST/evalue-'+str(evalue)+'/testsetE'+str(evalue)+'Fold'+str(i)+'.fasta',
                        'fasta'))
        for item in test_sequences:
            family = '.'.join(item.id.split('|')[3].split('.')[0:3])
            true_family[item.id] = family

        with open(f'../Results/BLAST/evalue-{evalue}/BLASTE{evalue}Fold{i}.list') as file:
            results = file.readlines()
            hit = results[0].split('\t')
            seq = hit[0]
            hit_list = []
            for result in results:
                hit = result.split('\t')
                if seq == hit[0] and len(hit_list) <= neighbours:
                    hit_list.append(hit[1])
                else:
                    counter = 0
                    max = hit_list[0]
                    for i in hit_list:
                        curr_frequency = hit_list.count(i)
                        if (curr_frequency > counter):
                            counter = curr_frequency
                            max = i
                    family = '.'.join(max.split('|')[3].split('.')[0:3])
                    pred_family[seq] = family
                    hit_list = []
                    hit_list.append(hit[1])
                    seq = hit[0]



#                    hit_family = '.'.join(hit[1].split('|')[3].split('.')[0:3])
 #                   pred_family[hit[0]] = hit_family

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
    #return accuracy_list, precision_list, recall_list, f1_list
    with open("../Results/BLAST/Unclassified/KNN_unclassified_evalue_" + str(evalue) + ".fasta",
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

def evaluation_subfamily(neighbours,evalue, cv):

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
            family = '.'.join(item.id.split('|')[3].split('.')[0:3])
            true_family[item.id] = family

        with open(f'../Results/BLAST/evalue-{evalue}/subfamily_BLASTE{evalue}Fold{i}.list') as file:
            results = file.readlines()
            hit = results[0].split('\t')
            seq = hit[0]
            hit_list = []
            for result in results:
                hit = result.split('\t')
                if seq == hit[0] and len(hit_list) <= neighbours:
                    hit_list.append(hit[1])
                else:
                    counter = 0
                    max = hit_list[0]
                    for i in hit_list:
                        curr_frequency = hit_list.count(i)
                        if (curr_frequency > counter):
                            counter = curr_frequency
                            max = i
                    family = '.'.join(max.split('|')[3].split('.')[0:3])
                    pred_family[seq] = family
                    hit_list = []
                    hit_list.append(hit[1])
                    seq = hit[0]



#                    hit_family = '.'.join(hit[1].split('|')[3].split('.')[0:3])
 #                   pred_family[hit[0]] = hit_family

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
    #return accuracy_list, precision_list, recall_list, f1_list
    with open("../Results/BLAST/Unclassified/KNN_subfamily_unclassified_evalue_" + str(evalue) + ".fasta",
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


accuracy_list = []
precision_list = []
recall_list = []
f1_list = []
total_unclassified = []
unclassified_list = []
#acc, precision, recall, f1 = evaluation(1, e, 5)

for k in range (1,31):
    acc, precision, recall, f1, unclassified = evaluation_family(k, 1e-30, 5)
  #  acc, precision, recall, f1, unclassified = evaluation_subfamily(k, 1e-30, 5)
    accuracy_list.append(acc)
    precision_list.append(precision)
    recall_list.append(recall)
    f1_list.append(f1)
    unclassified_list.append(unclassified)
k = [*range(1, 31 , 1)]
rows = zip(k,accuracy_list, precision_list, recall_list,f1_list,unclassified_list)

with open('../Results/BLAST/KNN_Blast_result_family_evalue_1e-30.csv', "w", newline='') as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)

#log_evalues = [log(y,10) for y in evalues]
'''folds = ['1', '2','3','4','5']
fig, axs = plt.subplots(2,2)
axs[0,0].bar(folds, acc, color=(0.2, 0.4, 0.6, 0.6))
plt.setp(axs, xlabel = 'Fold Number')
axs[0,0].set(ylabel = 'Accuracy')
axs[0,0].set_ylim([0.9, 1])
axs[1,0].set(ylabel = 'Precision')
axs[1,0].bar(folds, precision, color=(0.2, 0.4, 0.6, 0.6))
axs[0,1].set_ylim([0.9, 1])
axs[0,1].set(ylabel = 'Recall')
axs[0,1].bar(folds, recall, color=(0.2, 0.4, 0.6, 0.6))
axs[1,0].set_ylim([0.9, 1])
axs[1,1].set(ylabel = 'F1-score')
axs[1,1].bar(folds, f1, color=(0.2, 0.4, 0.6, 0.6))
axs[1,1].set_ylim([0.9, 1])
plt.show()
fig.savefig('../Results/Blast_5-fold_family.png')
'''

plt.xticks(k, rotation=90)
plt.plot(k, accuracy_list)
plt.plot(k, f1_list)
plt.plot(k, precision_list)
plt.plot(k, recall_list)
# plt.xlim(-20, 2)
plt.legend(['Accuracy', 'F1-Score', 'Precision', 'Recall'])
plt.xlabel('Number of Neighbors')
plt.ylabel('Percentage')
plt.grid()
plt.savefig('../Results/KNN_Blast_family_result_evalue_1e-30.png', dpi=300)
plt.show()
m =max(f1_list)
