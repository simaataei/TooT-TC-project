from sklearn import svm, metrics
True_class = [0,0,0,0,0,0,1,1,1,1,1,1]
BLAST_pred = [0,0,0,0,0,0,0,0,0,0,0,0]
HMM_pred =[0,0,0,0,0,0,0,0,0,1,0,0]
SVM_pred = [0,0,0,0,0,0,1,0,0,0,1,1]

accuracy_BLAST = metrics.accuracy_score(True_class, BLAST_pred)
precision_BLAST = metrics.precision_score(True_class , BLAST_pred, average="macro")
recall_BLAST = metrics.recall_score(True_class,  BLAST_pred, average="macro")
fscore_BLAST = metrics.f1_score(True_class, BLAST_pred, average="macro")
MCC_BLAST = metrics.matthews_corrcoef(True_class, BLAST_pred)


accuracy_HMM = metrics.accuracy_score(True_class, HMM_pred)
precision_HMM = metrics.precision_score(True_class , HMM_pred, average="macro")
recall_HMM = metrics.recall_score(True_class,  HMM_pred, average="macro")
fscore_HMM= metrics.f1_score(True_class, HMM_pred, average="macro")
MCC_HMM = metrics.matthews_corrcoef(True_class, HMM_pred)


accuracy_SVM = metrics.accuracy_score(True_class, SVM_pred)
precision_SVM = metrics.precision_score(True_class , SVM_pred, average="macro")
recall_SVM = metrics.recall_score(True_class,  SVM_pred, average="macro")
fscore_SVM= metrics.f1_score(True_class, SVM_pred, average="macro")
MCC_SVM = metrics.matthews_corrcoef(True_class, SVM_pred)
a=1