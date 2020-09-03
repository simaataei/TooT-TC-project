from sklearn import svm
from Bio import SeqIO
from classification.Extract_lables import extract_subfamilies, extract_families
from Feature_exctraction.AAC_extraction import extract_aac
from Feature_exctraction.PAAC_exraction import extract_paac
from Feature_exctraction.One_hot_encoding import extract_one_hot
from sklearn.model_selection import train_test_split


##################family classification##############
#read the data
family_data = list(SeqIO.parse("../Dataset/Selected30_families_tcdb2.fasta", "fasta"))


#extract lables
family_y = extract_families(family_data)




#extract features
family_aac_x = extract_aac(family_data)
family_paac_x = extract_paac(family_data)
family_onehot_x = extract_one_hot(family_data)

X_train_aac, X_test_aac, y_train_acc, y_test_aac = train_test_split(family_aac_x, family_y, test_size=0.2, random_state=42)

