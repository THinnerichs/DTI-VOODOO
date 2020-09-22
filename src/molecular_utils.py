import numpy as np

import pickle

import torch
import torch.nn as nn
import torch.utils.data as data

import DTI_data_preparation

import sys
from tqdm import tqdm




class MolecularDTIDataBuilder:
    def __init__(self, num_proteins=None, drug_mode='trfm'):
        # super(MolecularDTIData, self).__init__()
        print('Loading data...')
        self.drug_list = np.array(DTI_data_preparation.get_drug_list())
        print(len(self.drug_list), ' drugs present.')
        self.protein_list = np.array(DTI_data_preparation.get_human_proteins())[:num_proteins]
        print(len(self.protein_list), ' proteins present.')

        # build drug features
        drug_filename = '../models/drug_representation/' + drug_mode + '.npy'
        if drug_mode=='trfm' or drug_mode == 'rnn':
            self.drug_encodings = torch.from_numpy(np.load(drug_filename))
        else:
            print("No valid mode selected for drug to SMILES encoding.")
            raise ValueError

        # build protein features
        filename = '../models/protein_representation/results/prot_to_encoding_dict'
        with open(file=filename+'.pkl', mode='rb') as f:
            protein_to_feature_dict = pickle.load(f)
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
        for index in tqdm(indices):
            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            drug_encoding = self.drug_encodings[drug_index, :]
            protein_encoding = self.protein_encodings[protein_index, :]

            y = int(self.y_dti_data[drug_index, protein_index])
            return_list.append((torch.cat([protein_encoding, drug_encoding], 0), y))
        return return_list

    def __len__(self):
        return self.num_drugs * self.num_proteins

'''
class MolecularAminoacidMolPredDTI_DataBuilder:
    def __init__(self, num_proteins=None, drug_mode='trfm'):
        print('Loading data...')
        self.drug_list = np.array(DTI_data_preparation.get_drug_list())
        print(len(self.drug_list), ' drugs present.')
        self.protein_list = np.array(DTI_data_preparation.get_human_PhenomeNET_proteins())
        print(len(self.protein_list), ' proteins present.')


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
        for index in tqdm(indices):
            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            drug_encoding = self.drug_encodings[drug_index, :]
            protein_encoding = self.protein_encodings[protein_index, :]

            y = int(self.y_dti_data[drug_index, protein_index])
            return_list.append((torch.cat([protein_encoding, drug_encoding], 0), y))
        return return_list
'''


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

        self.fc1 = nn.Linear(8192 + 1024, 256)
        self.fc2 = nn.Linear(256,128)
        self.fc3 = nn.Linear(128,128)
        self.fc4 = nn.Linear(128,128)
        self.fc5 = nn.Linear(128, 1)

        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

        self.dropout = nn.Dropout(0.3)

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.relu(self.fc2(x))
        x = self.dropout(x)
        x = self.relu(self.fc3(x))
        x = self.dropout(x)
        x = self.relu(self.fc4(x))
        x = self.fc5(x)
        # x = self.sigmoid(x)

        return x


def train(model, device, train_loader, optimizer, epoch, weight_dict={0:1., 1:1.}):
    print('Training on {} samples...'.format(len(train_loader.dataset)))
    sys.stdout.flush()
    model.train()
    return_loss = 0
    for batch_idx, (features, labels) in enumerate(train_loader):
        optimizer.zero_grad()
        features = features.to(device)
        labels = labels.float().to(device)
        output = model(features)

        weight_vec = torch.ones([1]) * weight_dict[1]

        loss = nn.BCEWithLogitsLoss(pos_weight=weight_vec.to(output.device))(output, labels.view(-1, 1))
        return_loss += loss
        loss.backward()
        optimizer.step()
        if batch_idx % 10 == 0:
            print('Train epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(epoch,
                                                                           batch_idx * output.size(0),
                                                                           len(train_loader.dataset),
                                                                           100. * batch_idx / len(train_loader),
                                                                           loss.item()))
            sys.stdout.flush()
    return return_loss

def predicting(model, device, loader):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    print('Make prediction for {} samples...'.format(len(loader.dataset)))
    with torch.no_grad():
        for batch_idx, (features, labels) in enumerate(loader):
            features = features.to(device)
            output = model(features).sigmoid()
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            total_labels = torch.cat((total_labels, labels.view(-1, 1).float().cpu()), 0)

    return total_labels.round().numpy().flatten(),np.around(total_preds.numpy()).flatten()

