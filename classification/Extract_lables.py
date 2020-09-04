from Bio import SeqIO
import numpy as np


#family map
def extract_families(data):
    family_map = {} #key:family number, value: id
    i = 0
    j = 0
  #  z = 0
    y = np.zeros((len(data),1))
    #family_y = np.zeros((len(data),1))


    for item in data:
         famil = item.id.split('|')[3].split('.')[0:3]
         family = '.'.join(famil)
      #   family_y[z] = float(family)
         #z += 1
         if family not in family_map:
             family_map[family] = i
             i += 1

    for item in data:
        famil = item.id.split('|')[3].split('.')[0:3]
        family = '.'.join(famil)
        y[j] = family_map[family]
        j += 1

    return y

def extract_subfamilies(data):
    y = np.zeros((len(data), 1))
    sub_family_map = {}  # key:family number, value: id
    i = 0
    j = 0
    y = np.zeros((len(data),1))
    for item in data:
        sub_famil = item.id.split('|')[3].split('.')[0:4]
        subfamily = '.'.join(sub_famil)

        if subfamily not in sub_family_map:
            sub_family_map[subfamily] = i
            i += 1

    for item in data:
        sub_famil = item.id.split('|')[3].split('.')[0:4]
        subfamily = '.'.join(sub_famil)
        y[j] = sub_family_map[subfamily]
        j += 1


    return y
