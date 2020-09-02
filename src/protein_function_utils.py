import numpy as np

import gensim

import torch
import torch.nn as nn
import torch.utils.data as data

from tqdm import tqdm
import pickle
import sys

import DTI_data_preparation


class ProteinFunctionDTIDataBuilder:
    def __init__(self, config, num_proteins=None, drug_mode='trfm'):
        print('Loading protein function data...')
        self.drug_list = np.array(DTI_data_preparation.get_drug_list())
        print(len(self.drug_list), ' drugs present.')
        self.protein_list = np.array(DTI_data_preparation.get_human_prot_func_proteins())[:num_proteins]
        print(len(self.protein_list), ' proteins present.')

        self.include_uberon = config.include_uberon
        self.include_GO = config.include_GO
        self.include_phenotype = config.include_phenotype

        self.num_drugs = len(self.drug_list)
        self.num_proteins = len(self.protein_list)

        # build drug features
        drug_filename = '../models/drug_representation/' + drug_mode + '.npy'
        if drug_mode == 'trfm' or drug_mode == 'rnn':
            self.drug_encodings = torch.from_numpy(np.load(drug_filename))
        else:
            print("No valid mode selected for drug to SMILES encoding.")
            raise ValueError

        print('Building protein embeddings...')
        self.build_protein_embeddings()

        self.degree_features = DTI_data_preparation.get_protein_degree_percentile(self.protein_list, n=100)

        y_dti_data = DTI_data_preparation.get_DTIs(drug_list=self.drug_list, protein_list=self.protein_list)
        y_dti_data = y_dti_data.reshape((len(self.protein_list), len(self.drug_list)))
        self.y_dti_data = np.transpose(y_dti_data)

        print('Done.\n')

    def build_protein_embeddings(self):
        DL2vec_path_prefix = '../data/DL2vec/DL2vec_embeddings/'

        uberon_model_filename = DL2vec_path_prefix + 'uberon_intersection_ppi_embedding'
        GO_model_filename = DL2vec_path_prefix + 'go_intersection_ppi_embedding'
        phenotype_model_filename = DL2vec_path_prefix + 'mp_intersection_ppi_embedding'

        '''
        else:
            uberon_model_filename = DL2vec_path_prefix + 'uberon_intersection_embedding'
            GO_model_filename = DL2vec_path_prefix + 'go_intersection_embedding'
            phenotype_model_filename = DL2vec_path_prefix + 'mp_intersection_embedding'
        '''

        # load models
        uberon_model = gensim.models.Word2Vec.load(uberon_model_filename)
        GO_model = gensim.models.Word2Vec.load(GO_model_filename)
        phenotype_model = gensim.models.Word2Vec.load(phenotype_model_filename)

        # Build wordvector dicts
        uberon_model = uberon_model.wv
        GO_model = GO_model.wv
        phenotype_model = phenotype_model.wv

        uberon_embeddings = []
        GO_embeddings = []
        phenotype_embeddings = []
        for protein in self.protein_list:
            # organism, protein_id = protein.strip().split('.')
            protein_id = protein

            uberon_embeddings.append(uberon_model[protein_id])
            GO_embeddings.append(GO_model[protein_id])
            phenotype_embeddings.append(phenotype_model[protein_id])

        self.uberon_embeddings = torch.Tensor(uberon_embeddings)
        self.GO_embeddings = torch.Tensor(GO_embeddings)
        self.phenotype_embeddings = torch.Tensor(phenotype_embeddings)
        return np.array(uberon_embeddings), np.array(GO_embeddings), np.array(phenotype_embeddings)

    def get(self, indices):
        return_list = []
        for index in tqdm(indices):
            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            drug_encoding = self.drug_encodings[drug_index, :]

            y = int(self.y_dti_data[drug_index, protein_index])
            return_list.append((torch.cat([(self.uberon_embeddings[protein_index] if self.include_uberon else torch.zeros(self.uberon_embeddings[protein_index].size())),
                                           (self.GO_embeddings[protein_index] if self.include_GO else torch.zeros(self.GO_embeddings[protein_index].size())),
                                           (self.phenotype_embeddings[protein_index] if self.include_phenotype else torch.zeros(self.phenotype_embeddings[protein_index].size())),
                                           torch.Tensor(self.degree_features[protein_index,:]),
                                           drug_encoding], 0), y))
        return return_list

class ProteinFunctionDTIDataset(data.Dataset):
    def __init__(self, data):
        super(ProteinFunctionDTIDataset, self).__init__()
        self.data = data

    def __getitem__(self, index):
        return self.data[index]

    def __len__(self):
        return len(self.data)

class ProteinFunctionPredNet(nn.Module):
    def __init__(self):
        super(ProteinFunctionPredNet, self).__init__()

        self.fc1 = nn.Linear(600 + 100 + 1024, 256)
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
        x = self.dropout(x)
        x = self.fc5(x)
        # x = self.sigmoid(x)

        return x





