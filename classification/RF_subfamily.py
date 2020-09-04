from sklearn import svm, metrics
from Bio import SeqIO
import numpy as np
from sklearn.metrics import make_scorer
from classification.Extract_lables import extract_subfamilies, extract_families
from Feature_exctraction.AAC_extraction import extract_aac
from Feature_exctraction.PAAC_exraction import extract_paac
from Feature_exctraction.One_hot_encoding import extract_one_hot
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score, cross_validate
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.metrics import recall_score,f1_score
##################Subfamily classification##############
#read the data
Subfamily_data = list(SeqIO.parse("../Dataset/Selected30_subfamilies_tcdb2.fasta", "fasta"))


#extract lables
Subfamily_y = extract_subfamilies(Subfamily_data)
y = np.zeros(2940)



#extract features
Subfamily_aac_x = extract_aac(Subfamily_data)
Subfamily_paac_x = extract_paac(Subfamily_data)
Subfamily_onehot_x = extract_one_hot(Subfamily_data)


#5-fold cross validation

kf = KFold(n_splits=5, random_state=0, shuffle=True)
for train_index, test_index in kf.split(Subfamily_aac_x):
    print("TRAIN:", len(train_index), "TEST:", len(test_index))
    X_train, X_test = Subfamily_aac_x[train_index], Subfamily_aac_x[test_index]
    y_train, y_test = Subfamily_y[train_index], Subfamily_y[test_index]

clf = RandomForestClassifier(max_depth=20, random_state=0)
scoring = {'precision_macro','recall_macro','f1_macro'}
Linear_scores = cross_validate(clf, Subfamily_aac_x, Subfamily_y, scoring=scoring, cv=5)
f = open("RF_Subfamilies_aac.txt","w+")
f.write('max_depth=20\n')
f.write("Average Precision: "+str(sum(Linear_scores['test_precision_macro'])/5)+'\n')
f.write("Average Recall: "+str(sum(Linear_scores['test_recall_macro'])/5)+'\n')
f.write("Average f1-Measure: "+str(sum(Linear_scores['test_f1_macro'])/5)+'\n')
f.write(str(Linear_scores)+'\n')


###PAAC
#5-fold cross validation

kf = KFold(n_splits=5, random_state=0, shuffle=True)
for train_index, test_index in kf.split(Subfamily_paac_x):
    print("TRAIN:", len(train_index), "TEST:", len(test_index))
    X_train, X_test = Subfamily_paac_x[train_index], Subfamily_paac_x[test_index]
    y_train, y_test = Subfamily_y[train_index], Subfamily_y[test_index]

clf = RandomForestClassifier(max_depth=20, random_state=0)
scoring = {'precision_macro','recall_macro','f1_macro'}
Linear_scores = cross_validate(clf, Subfamily_paac_x, Subfamily_y, scoring=scoring, cv=5)
f1 = open("RF_Subfamilies_PAAC.txt","w+")
f1.write('max_depth=20\n')
f1.write("Average Precision: "+str(sum(Linear_scores['test_precision_macro'])/5)+'\n')
f1.write("Average Recall: "+str(sum(Linear_scores['test_recall_macro'])/5)+'\n')
f1.write("Average f1-Measure: "+str(sum(Linear_scores['test_f1_macro'])/5)+'\n')
f1.write(str(Linear_scores)+'\n')

#################one_hot
#5-fold cross validation

kf = KFold(n_splits=5, random_state=0, shuffle=True)
for train_index, test_index in kf.split(Subfamily_onehot_x ):
    print("TRAIN:", len(train_index), "TEST:", len(test_index))
    X_train, X_test = Subfamily_onehot_x [train_index], Subfamily_onehot_x [test_index]
    y_train, y_test = Subfamily_y[train_index], Subfamily_y[test_index]


clf = RandomForestClassifier(max_depth=20, random_state=0)
scoring = {'precision_macro','recall_macro','f1_macro'}
Linear_scores = cross_validate(clf, Subfamily_onehot_x , Subfamily_y, scoring=scoring, cv=5)
f1 = open("RF_Subfamilies_onehot.txt","w+")
f1.write('max_depth=20\n')
f1.write("Average Precision: "+str(sum(Linear_scores['test_precision_macro'])/5)+'\n')
f1.write("Average Recall: "+str(sum(Linear_scores['test_recall_macro'])/5)+'\n')
f1.write("Average f1-Measure: "+str(sum(Linear_scores['test_f1_macro'])/5)+'\n')
f1.write(str(Linear_scores)+'\n')
