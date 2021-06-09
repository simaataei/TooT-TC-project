from Bio import SeqIO
from classification.Extract_lables import extract_subfamilies, extract_families
from Feature_exctraction.AAC_extraction import extract_aac
from Feature_exctraction.PAAC_exraction import extract_paac
from Feature_exctraction.One_hot_encoding import extract_one_hot
import numpy as np
import sklearn.svm
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
import csv
from sklearn.neighbors import KNeighborsClassifier


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
test_data = list(SeqIO.parse("../Results/BLAST/Unclassified/unclassified_evalue_1e-30.fasta", "fasta"))


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

with open('../Results/BLAST/unclassified_classification_1e-30.csv', "w", newline='') as f:
    writer = csv.writer(f, delimiter=',')
    result =[]
    result.append("SVM")
    result.append("Kernel = RBF")
    rbf_svc = sklearn.svm.SVC(kernel='rbf')
    rbf_svc.fit(family_aac_x_train, y_train )
    predicted = rbf_svc.predict(family_aac_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))



    #PAAC
    rbf_svc = sklearn.svm.SVC(kernel='rbf')
    rbf_svc.fit(family_paac_x_train, y_train )
    predicted = rbf_svc.predict(family_paac_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))




    #one_hot
    rbf_svc = sklearn.svm.SVC(kernel='rbf')
    rbf_svc.fit(family_onehot_x_train, y_train )
    predicted = rbf_svc.predict(family_onehot_x_test)
    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    writer.writerow(result)


#Sigmoid

    result = []
    result.append("SVM")
    result.append("Kernel = Sigmoid")
    linear_svc = sklearn.svm.SVC(kernel='sigmoid')
    linear_svc.fit(family_aac_x_train, y_train)
    predicted = linear_svc.predict(family_aac_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    # PAAC
    linear_svc = sklearn.svm.SVC(kernel='sigmoid')
    linear_svc.fit(family_paac_x_train, y_train)
    predicted = linear_svc.predict(family_paac_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    # one_hot
    linear_svc = sklearn.svm.SVC(kernel='sigmoid')
    linear_svc.fit(family_onehot_x_train, y_train)
    predicted = linear_svc.predict(family_onehot_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    writer.writerow(result)

#Poly

    result = []
    result.append("SVM")
    result.append("Kernel = Polynomial")
    linear_svc = sklearn.svm.SVC(kernel='poly')
    linear_svc.fit(family_aac_x_train, y_train)
    predicted = linear_svc.predict(family_aac_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    # PAAC
    linear_svc = sklearn.svm.SVC(kernel='poly')
    linear_svc.fit(family_paac_x_train, y_train)
    predicted = linear_svc.predict(family_paac_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    # one_hot
    linear_svc = sklearn.svm.SVC(kernel='poly')
    linear_svc.fit(family_onehot_x_train, y_train)
    predicted = linear_svc.predict(family_onehot_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    writer.writerow(result)

    #Linear
#AAC
    result= []
    result.append("SVM")
    result.append("Kernel = Linear")
    linear_svc = sklearn.svm.SVC(kernel='linear')
    linear_svc.fit(family_aac_x_train, y_train )
    predicted = linear_svc.predict(family_aac_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    #PAAC
    linear_svc = sklearn.svm.SVC(kernel='linear')
    linear_svc.fit(family_paac_x_train, y_train )
    predicted = linear_svc.predict(family_paac_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    #one_hot
    linear_svc = sklearn.svm.SVC(kernel='linear')
    linear_svc.fit(family_onehot_x_train, y_train )
    predicted = linear_svc.predict(family_onehot_x_test)

    result.append(metrics.f1_score(y_test, predicted, average="macro"))
    writer.writerow(result)

#Random forest
    result = []
    result.append("Random Forest")
    result.append("Max Depth =30")
#AAC
    clf = RandomForestClassifier(max_depth=30, random_state=0)
    clf.fit(family_aac_x_train, y_train )
    predicted = clf.predict(family_aac_x_test)
    result.append(metrics.f1_score(y_test, predicted, average="macro"))

#PAAC
    clf = RandomForestClassifier(max_depth=30, random_state=0)
    clf.fit(family_paac_x_train, y_train)
    predicted = clf.predict(family_paac_x_test)
    result.append(metrics.f1_score(y_test, predicted, average="macro"))

#one_hot
    clf = RandomForestClassifier(max_depth=30, random_state=0)
    clf.fit(family_onehot_x_train, y_train)
    predicted = clf.predict(family_onehot_x_test)
    result.append(metrics.f1_score(y_test, predicted, average="macro"))
    writer.writerow(result)

#KNN
    result = []
    result.append("KNN")
    result.append("Number of Neighbors =5")
#AAC
    neigh = KNeighborsClassifier(n_neighbors=3)
    neigh.fit(family_aac_x_train, y_train )
    predicted = neigh.predict(family_aac_x_test)
    result.append(metrics.f1_score(y_test, predicted, average="macro"))

#PAAC
    neigh = KNeighborsClassifier(n_neighbors=3)
    neigh.fit(family_paac_x_train, y_train)
    predicted = neigh.predict(family_paac_x_test)
    result.append(metrics.f1_score(y_test, predicted, average="macro"))

#one_hot
    neigh = KNeighborsClassifier(n_neighbors=3)
    neigh.fit(family_onehot_x_train, y_train)
    predicted = neigh.predict(family_onehot_x_test)
    result.append(metrics.f1_score(y_test, predicted, average="macro"))

    writer.writerow(result)

