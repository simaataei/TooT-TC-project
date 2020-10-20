from Bio import SeqIO
from classification.Extract_lables import extract_subfamilies, extract_families
from Feature_exctraction.AAC_extraction import extract_aac
from Feature_exctraction.PAAC_exraction import extract_paac
from Feature_exctraction.One_hot_encoding import extract_one_hot
import numpy as np
import sklearn.svm
import matplotlib.pyplot as plt
from sklearn import metrics

def count_families(data):
    family_count ={}
    for item in data:
        famil = item.id.split('|')[3].split('.')[0:3]
        family = '.'.join(famil)
        if family not in family_count.keys():
            family_count[family] = 1
        else:
            family_count[family] += 1
    return family_count





family_data = list(SeqIO.parse("../Dataset/Selected30_families_tcdb2.fasta", "fasta"))
test_data = list(SeqIO.parse("../Results/BLAST/Unclassified/unclassified_evalue_0.0001.fasta", "fasta"))


unclassified = {}
train_data = []
test_list = []

for t in test_data:
    test_list.append(t.id)
for seq in family_data:
    if seq.id not in test_list:
        train_data.append(seq)

'''family_count ={}
family_count = count_families(test_data)
plt.bar(family_count.keys(), family_count.values(),color=(0.1, 0.6, 0.8, 0.8))
plt.grid()
plt.xticks(rotation= 90, fontsize= 6)
plt.xlabel("family", fontsize = 8)
plt.ylabel("Number of unclassified sequences")
plt.savefig('../Results/unclassified_blast_0.0001.png', dpi=300)
plt.show()
'''
#extract lables
y_test = []
y_train = []

family_train =  extract_families(train_data)
test_families =  extract_families(test_data)
j=0
for item in test_families:
    for i in item:
        y_test.append(i)
y_test = np.array(y_test)



for item in family_train:
    for i in item:
        y_train.append(i)
y_train = np.array(y_train)
#family_y = extract_families(family_data)
y = np.zeros(2940)



#extract features_train
family_aac_x_train = extract_aac(train_data)
family_paac_x_train = extract_paac(train_data)
family_onehot_x_train = extract_one_hot(train_data)


#extract features_test
family_aac_x_test = extract_aac(test_data)
family_paac_x_test = extract_paac(test_data)
family_onehot_x_test = extract_one_hot(test_data)


#SVMs
#RBF
#AAC
rbf_svc = sklearn.svm.SVC(kernel='rbf')
rbf_svc.fit(family_aac_x_train, y_train )
predicted = rbf_svc.predict(family_aac_x_test)

f1_rbf = metrics.f1_score(y_test, predicted, average="macro")


#PAAC
rbf_svc = sklearn.svm.SVC(kernel='rbf')
rbf_svc.fit(family_paac_x_train, y_train )
predicted = rbf_svc.predict(family_paac_x_test)

f1_rbf = metrics.f1_score(y_test, predicted, average="macro")

#one_hot
rbf_svc = sklearn.svm.SVC(kernel='rbf')
rbf_svc.fit(family_onehot_x_train, y_train )
predicted = rbf_svc.predict(family_onehot_x_test)

f1_rbf = metrics.f1_score(y_test, predicted, average="macro")

linear_svc = sklearn.svm.SVC(kernel='linear')
linear_svc.fit(family_onehot_x_train, y_train )
predicted = linear_svc.predict(family_onehot_x_test)

f1_linear = metrics.f1_score(y_test, predicted, average="macro")
a=1



