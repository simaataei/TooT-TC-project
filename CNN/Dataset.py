import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
#mport deepchem as dc




class NumbersDataset(Dataset):
    def __init__(self):
        self.samples = list(range(1, 1001))

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        return self.samples[idx]



dataset = NumbersDataset()
print(len(dataset))
print(dataset[100])
print(dataset[122:361])

a=1