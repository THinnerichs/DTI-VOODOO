import numpy as np

import pickle

import torch
import torch.nn as nn
import torch.utils.data as data

import DTI_data_preparation




def train():


class MolecularDTIDataBuilder:
    def __init__(self, num_proteins=None, drug_mode='trfm'):
        # super(MolecularDTIData, self).__init__()
        print('Loading data.')
        self.drug_list = np.array(DTI_data_preparation.get_drug_list())
        self.protein_list = np.array(DTI_data_preparation.get_human_proteins())[:num_proteins]

        # build drug features
        drug_filename = '../models/drug_representation/' + drug_mode + '.npy'
        if drug_mode=='trfm' or drug_mode == 'rnn':
            self.drug_encodings = torch.from_numpy(np.load(drug_filename))
        else:
            print("No valid mode selected for drug to SMILES encoding.")
            raise ValueError

        # build protein features
        filename = '../models/protein_representation/results/prot_to_encoding_dict'
        protein_to_feature_dict = pickle.load(filename)
        self.protein_encodings = torch.Tensor([protein_to_feature_dict[protein] for protein in self.protein_list])

        y_dti_data = DTI_data_preparation.get_DTIs(drug_list=self.drug_list, protein_list=self.protein_list)
        y_dti_data = y_dti_data.reshape((len(self.protein_list), len(self.drug_list)))
        self.y_dti_data = np.transpose(y_dti_data)

        # calculate dimenions of data
        self.num_proteins = len(self.protein_list)
        self.num_drugs = len(self.drug_list)
        print('Done.\n')

    def get(self, indices):
        return_list = []
        for index in indices:
            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            drug_encoding = self.drug_encodings[drug_index, :]
            protein_encoding = self.protein_encodings[protein_index, :]

            y = int(self.y_dti_data[drug_index, protein_index])
            return_list.append((torch.cat([protein_encoding, drug_encoding], 0), y))
        return return_list

    def __len__(self):
        return self.num_drugs * self.num_proteins

class MolecularDTIDataset(data.Dataset):
    def __init__(self, data):
        super(MolecularDTIDataset, self).__init__()
        self.data = data

    def __getitem__(self, index):
        return self.data[index]

    def __len__(self):
        return len(self.data)



class MolecularPredNet(nn.Module):
    def __init__(self):
        super(MolecularPredNet, self).__init__()

        self.fc1 = nn.Linear(8192 + 1024, 128)
        self.fc2 = nn.Linear(128, 1)

        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.sigmoid(self.fc2(x))

        return x

