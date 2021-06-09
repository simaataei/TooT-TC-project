from Bio import SeqIO
from sklearn import metrics
def BLAST_evaluation():


    total_unclassified =[]
    unclassified = []
    true_family ={}
    pred_family = {}

    y_true = []
    y_pred = []
    test_sequences = list(
    SeqIO.parse('..\Dataset\\family_3A7\Balanced30_testset.fasta','fasta'))
    for item in test_sequences:
        family = '.'.join(item.id.split('|')[3].split('.')[0:3])
       # true_family[item.id] = family
        if family.startswith("3.A.7"):
            true_family[item.id] = 1
        else :
            true_family[item.id] = 0


    with open('..\Results\\family_analysis\\family_3A7\Balanced_test_on_train_diff_scov_without_selfhit.list') as file:
        results = file.readlines()
        for result in results:
            hit = result.split('\t')
            if hit[0] not in pred_family.keys():
                hit_family = '.'.join(hit[1].split('|')[3].split('.')[0:3])

                if hit_family.startswith("3.A.7"):
                    pred_family[hit[0]] = 1
                else:
                    pred_family[hit[0]] = 0
  #              pred_family[hit[0]] = hit_family

        for key in true_family.keys():
            if key not in pred_family.keys():
                unclassified.append(key)

            else:
                y_true.append(true_family[key])
                y_pred.append(pred_family[key])

    return [metrics.accuracy_score(y_true, y_pred), metrics.precision_score(y_true, y_pred, average="macro"),metrics.recall_score(y_true, y_pred, average="macro"),metrics.f1_score(y_true, y_pred, average="macro")]

a = BLAST_evaluation()
#tn, fp, fn, tp = a.ravel()

s=2

