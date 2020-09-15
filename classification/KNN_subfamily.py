
from sklearn import svm, metrics
from Bio import SeqIO
import numpy as np
from sklearn.metrics import make_scorer
from classification.Extract_lables import extract_families, extract_subfamilies
from Feature_exctraction.AAC_extraction import extract_aac
from Feature_exctraction.PAAC_exraction import extract_paac
from Feature_exctraction.One_hot_encoding import extract_one_hot
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.metrics import recall_score,f1_score
##################subfamily classification##############
#read the data
subfamily_data = list(SeqIO.parse("../Dataset/Selected30_subfamilies_tcdb2.fasta", "fasta"))


#extract lables
subfamily_y = extract_subfamilies(subfamily_data)
y = np.zeros(2940)



#extract features
subfamily_aac_x = extract_aac(subfamily_data)
subfamily_paac_x = extract_paac(subfamily_data)
subfamily_onehot_x = extract_one_hot(subfamily_data)


#5-fold cross validation
#linear kernel
f = open("../Results/KNN_subfamilies_aac.txt", "w+")
Neighbor_n =[1, 5, 10, 20]
for n in Neighbor_n:
    clf = KNeighborsClassifier(n_neighbors=n)

    scoring = {'precision_macro','recall_macro','f1_macro'}
    Linear_scores = cross_validate(clf, subfamily_aac_x, subfamily_y, scoring=scoring, cv=5)

    f.write('n_neighbors=' + str(n) + '\n')
    f.write("Average Precision: "+str(sum(Linear_scores['test_precision_macro'])/5)+'\n')
    f.write("Average Recall: "+str(sum(Linear_scores['test_recall_macro'])/5)+'\n')
    f.write("Average f1-Measure: "+str(sum(Linear_scores['test_f1_macro'])/5)+'\n')
    f.write(str(Linear_scores)+'\n')


###PAAC
#5-fold cross validation
f1 = open("../Results/KNN_subfamilies_PAAC.txt", "w+")
for n in Neighbor_n:
    clf = KNeighborsClassifier(n_neighbors=n)
    scoring = {'precision_macro','recall_macro','f1_macro'}
    Linear_scores = cross_validate(clf, subfamily_paac_x, subfamily_y, scoring=scoring, cv=5)

    f1.write('n_neighbors='+str(n)+'\n')
    f1.write("Average Precision: "+str(sum(Linear_scores['test_precision_macro'])/5)+'\n')
    f1.write("Average Recall: "+str(sum(Linear_scores['test_recall_macro'])/5)+'\n')
    f1.write("Average f1-Measure: "+str(sum(Linear_scores['test_f1_macro'])/5)+'\n')
    f1.write(str(Linear_scores)+'\n')

#################one_hot
#5-fold cross validation
#linear kernel
f2 = open("../Results/KNN_subfamilies_onehot.txt", "w+")
for n in Neighbor_n:
    clf = KNeighborsClassifier(n_neighbors=n)
    scoring = {'precision_macro','recall_macro','f1_macro'}
    Linear_scores = cross_validate(clf, subfamily_onehot_x , subfamily_y, scoring=scoring, cv=5)

    f2.write('n_neighbors=' + str(n) + '\n')
    f2.write("Average Precision: "+str(sum(Linear_scores['test_precision_macro'])/5)+'\n')
    f2.write("Average Recall: "+str(sum(Linear_scores['test_recall_macro'])/5)+'\n')
    f2.write("Average f1-Measure: "+str(sum(Linear_scores['test_f1_macro'])/5)+'\n')
    f2.write(str(Linear_scores)+'\n')
