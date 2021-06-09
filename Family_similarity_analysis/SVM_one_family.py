from Bio import SeqIO
import random
from Feature_exctraction.AAC_extraction import extract_aac
from Feature_exctraction.PAAC_exraction import extract_paac
from Feature_exctraction.One_hot_encoding import extract_one_hot
from classification.Extract_lables import  extract_families
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn.metrics import make_scorer
from sklearn import svm, metrics


train_data  = list(SeqIO.parse("Dataset/family_3A7/Balanced30_trainingset.fasta", 'fasta'))

test_data = list( SeqIO.parse('Dataset/family_3A7/Balanced30_testset.fasta',
                        'fasta'))

y_train = np.zeros(len(train_data))
y_test = np.zeros(len(test_data))
i=0
for item in train_data:
    hit_family = '.'.join(item.id.split('|')[3].split('.')[0:3])
    if hit_family.startswith("3.A.7"):
        y_train[i] = 1
    i+=1

i=0
for item in test_data:
    hit_family = '.'.join(item.id.split('|')[3].split('.')[0:3])
    if hit_family.startswith("3.A.7"):
        y_test[i] = 1
    i+=1


#extract features
train_aac_x = extract_aac(train_data)
train_paac_x = extract_paac(train_data)
train_onehot_x = extract_one_hot(train_data)

test_aac_x = extract_aac(test_data)
test_paac_x = extract_paac(test_data)
test_onehot_x = extract_one_hot(test_data)




#5-fold cross validation
#linear kernel
#kf = KFold(n_splits=5, random_state=0, shuffle=True)
#for train_index, test_index in kf.split(family_onehot_x):
#    print("TRAIN:", len(train_index), "TEST:", len(test_index))

clf = svm.SVC(kernel='linear',probability=True) # Linear Kernel

#Train the model using the training sets
clf.fit(train_paac_x , y_train)

#Predict the response for test dataset
y_pred = clf.predict(test_paac_x)

accuracy = metrics.accuracy_score(y_test, y_pred)
precision = metrics.precision_score(y_test, y_pred, average="macro")
recall = metrics.recall_score(y_test, y_pred, average="macro")
fscore = metrics.f1_score(y_test, y_pred, average="macro")
MCC = metrics.matthews_corrcoef(y_test, y_pred)
Prob = clf.predict_proba(test_paac_x)
pred = clf.predict(test_paac_x)
decision =clf.decision_function(test_paac_x)

a=1

f1 = open("Results/family_analysis/family_3A7/SVM/BestResults.txt","w+")
f1.write('kernel=linear, C=1\n')
f1.write("Precision: "+str(precision)+'\n')
f1.write("recall: "+str(recall)+'\n')
f1.write("F-score: "+str(fscore)+'\n')
f1.write("Accuracy: "+str(accuracy)+'\n')
f1.write("MCC:"+str(MCC)+'\n')
f1.write("results:\n")
i=0
a =clf.predict(test_onehot_x)
b= clf.predict_proba(test_onehot_x)
for item in test_data:
    hit_family = '.'.join(item.id.split('|')[3].split('.')[0:3])

    f1.write(str(hit_family)+', '+str(y_pred[i])+ ', ' + str(Prob[i])+', '+str(clf)+'\n')

    i += 1



a=1
#RBF


clf = svm.SVC(kernel='rbf', probability =True)
clf.fit(train_paac_x , y_train)

#Predict the response for test dataset
y_pred = clf.predict(test_paac_x)
accuracy = metrics.accuracy_score(y_test, y_pred)
precision = metrics.precision_score(y_test, y_pred, average="macro")
recall = metrics.recall_score(y_test, y_pred, average="macro")
fscore = metrics.f1_score(y_test, y_pred, average="macro")
MCC =metrics.matthews_corrcoef(y_test, y_pred)
Prob = clf.predict_proba(test_paac_x)
i=0

f1.write('kernel=RBF\n')
f1.write("Precision: "+str(precision)+'\n')
f1.write("recall: "+str(recall)+'\n')
f1.write("F-score: "+str(fscore)+'\n')
f1.write("Accuracy: "+str(accuracy)+'\n')
f1.write("MCC:"+str(MCC)+'\n')
f1.write("results:\n")
for item in test_data:
    hit_family = '.'.join(item.id.split('|')[3].split('.')[0:3])
    f1.write(str(hit_family) + ', ' + str(y_pred[i]) + ', ' + str(Prob[i]) + ', ' + str(clf) + '\n')
    i += 1


#poly

clf = svm.SVC(kernel='poly', probability=True)
clf.fit(train_aac_x , y_train)

#Predict the response for test dataset
y_pred = clf.predict(test_aac_x)

accuracy = metrics.accuracy_score(y_test, y_pred)
precision = metrics.precision_score(y_test, y_pred, average="macro")
recall = metrics.recall_score(y_test, y_pred, average="macro")
fscore = metrics.f1_score(y_test, y_pred, average="macro")
Prob = clf.predict_proba(test_aac_x)

i=0
f1.write('kernel=poly\n')
f1.write("Precision: "+str(precision)+'\n')
f1.write("recall: "+str(recall)+'\n')
f1.write("F-score: "+str(fscore)+'\n')
f1.write("Accuracy: "+str(accuracy)+'\n')
f1.write("MCC:"+str(MCC)+'\n')
f1.write("results:\n")
for item in test_data:
    hit_family = '.'.join(item.id.split('|')[3].split('.')[0:3])
    f1.write(str(hit_family) + ', ' + str(y_pred[i]) + ', ' + str(Prob[i]) + ', ' + str(clf) + '\n')
    i += 1

#sigmoid
clf = svm.SVC(kernel='sigmoid', probability=True)
clf.fit(train_onehot_x , y_train)

#Predict the response for test dataset
y_pred = clf.predict(test_onehot_x)

accuracy = metrics.accuracy_score(y_test, y_pred)
precision = metrics.precision_score(y_test, y_pred, average="macro")
recall = metrics.recall_score(y_test, y_pred, average="macro")
fscore = metrics.f1_score(y_test, y_pred, average="macro")
MCC =metrics.matthews_corrcoef(y_test, y_pred)
Prob = clf.predict_proba(test_onehot_x)
i=0
f1.write('kernel=sigmoid\n')
f1.write("Precision: "+str(precision)+'\n')
f1.write("recall: "+str(recall)+'\n')
f1.write("F-score: "+str(fscore)+'\n')
f1.write("Accuracy: "+str(accuracy)+'\n')
f1.write("MCC:"+str(MCC)+'\n')
f1.write("results:\n")
for item in test_data:
    hit_family = '.'.join(item.id.split('|')[3].split('.')[0:3])
    f1.write(str(hit_family) + ', ' + str(y_pred[i]) + ', ' + str(Prob[i]) + ', ' + str(clf) + '\n')
    i += 1
f1.close()
a=1