import numpy as np
from math import sqrt, log2
from scipy import stats

import torch
import torch.nn as nn
import torch.nn.functional as F

from torch_geometric.data import Dataset, Data, InMemoryDataset

from tqdm import tqdm
import sys
import pickle

import DTI_data_preparation
from PPI_utils import get_PPI_graph
from protein_function_utils import ProteinFunctionDTIDataBuilder



class SimpleDTINetworkData():
    """
    Args:
        root (string): Root directory where the dataset should be saved.
        name (string): The name of the dataset (:obj:`"Cora"`,
            :obj:`"CiteSeer"`, :obj:`"PubMed"`).
        transform (callable, optional): A function/transform that takes in an
            :obj:`torch_geometric.data.Data` object and returns a transformed
            version. The data object will be transformed before every access.
            (default: :obj:`None`)
        pre_transform (callable, optional): A function/transform that takes in
            an :obj:`torch_geometric.data.Data` object and returns a
            transformed version. The data object will be transformed before
            being saved to disk. (default: :obj:`None`)
    """

    def __init__(self, num_proteins=None, min_score=700, transform=None, pre_transform=None):
        # super(FullNetworkDataset, self).__init__(root='../data/torch_raw/')

        print("Loading data ...")
        self.drug_list = np.array(DTI_data_preparation.get_drug_list())
        print(len(self.drug_list), "drugs present")
        self.protein_list = np.array(DTI_data_preparation.get_human_proteins())[:num_proteins]
        print(len(self.protein_list), "proteins present\n")

        # PPI data
        print("Loading PPI graph ...")
        PPI_graph = get_PPI_graph(min_score=min_score)
        PPI_graph = PPI_graph.subgraph(self.protein_list)

        print("Building index dict ...")
        self.protein_to_index_dict = {protein: index for index, protein in enumerate(self.protein_list)}
        print("Building edge list ...")
        forward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2])
                              for node1, node2 in list(PPI_graph.edges())]
        backward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2])
                               for node2, node1 in list(PPI_graph.edges())]
        self.edge_list = torch.tensor(np.transpose(np.array(forward_edges_list + backward_edges_list)), dtype=torch.long)
        print("Building feature matrix ...")
        self.feature_matrix = torch.tensor(DTI_data_preparation.get_PPI_node_feature_mat_list(self.protein_list), dtype=torch.float)
        print("feature_matrix.size()", self.feature_matrix.size())
        self.num_PPI_features = self.feature_matrix.shape[1]

        # DDI data
        print("Loading DDI features ...")
        self.DDI_features = DTI_data_preparation.get_DDI_feature_list(self.drug_list)
        print(self.DDI_features.shape)

        # DTI data
        print("Loading DTI links ...")
        y_dti_data = DTI_data_preparation.get_DTIs(drug_list=self.drug_list, protein_list=self.protein_list)
        y_dti_data = y_dti_data.reshape((len(self.protein_list), len(self.drug_list)))
        self.y_dti_data = np.transpose(y_dti_data)
        print(self.y_dti_data.shape)

        # calculate dimensions of network
        self.num_proteins = len(PPI_graph.nodes())
        self.num_drugs = len(self.drug_list)
        print("Finished.\n")

    def get_protein_to_index_dict(self):
        return self.protein_to_index_dict

    def get_protein_list(self):
        return self.protein_list

    def get_drug_list(self):
        return self.drug_list

    def get(self, indices):
        data_list = []

        # for index in tqdm(indices):
        for i in range(len(indices)):
            if (i+1)%int(len(indices)/10) == 0:
                print('Finished {} percent.'.format(str(int(i/len(indices)*100))), end='\r')

            index = indices[i]
            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            # build protein mask
            protein_mask = np.zeros(self.num_proteins)
            protein_mask[protein_index] = 1
            protein_mask = torch.tensor(protein_mask, dtype=torch.bool)

            y = int(self.y_dti_data[drug_index, protein_index])

            DDI_features = torch.tensor(self.DDI_features[:, drug_index], dtype=torch.float).view(1, self.num_drugs)

            full_PPI_graph = Data(x=self.feature_matrix, edge_index=self.edge_list, y=y)
            full_PPI_graph.DDI_features = DDI_features
            full_PPI_graph.protein_mask = protein_mask
            # full_PPI_graph.__num_nodes__ = self.num_proteins

            data_list.append(full_PPI_graph)

        return data_list

    def __len__(self):
        return self.num_proteins * self.num_drugs

class MolPredDTINetworkData():
    def __init__(self, config):
        self.config = config

        if config.num_proteins == -1 or config.num_proteins == 11574:
            config.num_proteins = None


        print("Loading data ...")
        self.drug_list = np.array(DTI_data_preparation.get_drug_list())
        print(len(self.drug_list), "drugs present")
        self.protein_list = np.array(DTI_data_preparation.get_human_proteins())[:config.num_proteins]
        print(len(self.protein_list), "proteins present\n")

        # PPI data
        print("Loading PPI graph ...")
        PPI_graph = DTI_data_preparation.get_PPI_DTI_graph_intersection()
        PPI_graph = PPI_graph.subgraph(self.protein_list)

        # calculate dimensions of network
        self.num_proteins = len(PPI_graph.nodes())
        self.num_drugs = len(self.drug_list)

        print("Building index dict ...")
        self.protein_to_index_dict = {protein: index for index, protein in enumerate(self.protein_list)}
        print("Building edge list ...")
        forward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2])
                              for node1, node2 in list(PPI_graph.edges())]
        backward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2])
                               for node2, node1 in list(PPI_graph.edges())]
        self.edge_list = torch.tensor(np.transpose(np.array(forward_edges_list + backward_edges_list)),
                                      dtype=torch.long)
        print("Building feature matrix ...")
        filename = '../models/molecular_predictor/pred_fold_' + str(config.fold)
        with open(file=filename+'.pkl', mode='rb') as f:
            self.mol_predictions = pickle.load(f).reshape((self.num_drugs, -1))
        self.num_PPI_features = 1

        if config.drug_mode == 'semsim' or self.config.drug_mode == 'all':
            print("Loading semsim results...")
            self.semsim_matrix = DTI_data_preparation.get_side_effect_similarity_feature_list(self.drug_list)

        # DDI data
        print("Loading DDI features ...")
        self.DDI_features = DTI_data_preparation.get_DDI_feature_list(self.drug_list)
        print(self.DDI_features.shape)

        # DTI data
        print("Loading DTI links ...")
        y_dti_data = DTI_data_preparation.get_DTIs(drug_list=self.drug_list, protein_list=self.protein_list)
        y_dti_data = y_dti_data.reshape((len(self.protein_list), len(self.drug_list)))
        self.y_dti_data = np.transpose(y_dti_data)
        print(self.y_dti_data.shape)
        print("Finished.\n")

    def get_protein_to_index_dict(self):
        return self.protein_to_index_dict

    def get_protein_list(self):
        return self.protein_list

    def get_drug_list(self):
        return self.drug_list

    def get(self, indices):
        data_list = []

        # for index in tqdm(indices):
        for i in range(len(indices)):
            if (i + 1) % int(len(indices) / 10) == 0:
                print('Finished {} percent.'.format(str(int(i / len(indices) * 100))), end='\r')

            index = indices[i]
            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            # build protein mask
            protein_mask = np.zeros(self.num_proteins)
            protein_mask[protein_index] = 1
            protein_mask = torch.tensor(protein_mask, dtype=torch.bool)

            y = int(self.y_dti_data[drug_index, protein_index])

            DDI_features = torch.tensor(self.DDI_features[:, drug_index], dtype=torch.float).view(1, self.num_drugs)
            semsim_features = None
            if self.config.drug_mode == 'semsim' or self.config.drug_mode == 'all':
                semsim_features = torch.tensor(self.semsim_matrix[drug_index, :], dtype=torch.float).view(1, self.num_drugs)


            feature_array = torch.tensor(self.mol_predictions[drug_index, :self.num_proteins], dtype=torch.float).round().view(-1, 1)
            full_PPI_graph = Data(x=feature_array, edge_index=self.edge_list, y=y)
            full_PPI_graph.DDI_features = DDI_features
            full_PPI_graph.protein_mask = protein_mask
            if self.config.drug_mode == 'semsim' or self.config.drug_mode == 'all':
                full_PPI_graph.semsim_features = semsim_features
            # full_PPI_graph.__num_nodes__ = self.num_proteins

            data_list.append(full_PPI_graph)

        return data_list

    def __len__(self):
        return self.num_proteins * self.num_drugs

    def get_labels(self, indices):
        data_list = []

        # for index in tqdm(indices):
        for i in range(len(indices)):
            if (i + 1) % int(len(indices) / 10) == 0:
                print('Finished {} percent of labels.'.format(str(int(i / len(indices) * 100))), end='\r')

            index = indices[i]

            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins
            y = int(self.y_dti_data[drug_index, protein_index])
            data_list.append(y)

        return np.array(data_list)

class ProtFuncDTINetworkData:
    def __init__(self, config):
        self.config = config

        print("Loading data ...")
        self.drug_list = np.array(DTI_data_preparation.get_drug_list())
        print(len(self.drug_list), "drugs present")
        self.protein_list = np.array(DTI_data_preparation.get_human_prot_func_proteins())[:config.num_proteins]
        print(len(self.protein_list), "proteins present\n")

        # PPI data
        print("Loading PPI graph ...")
        PPI_graph = DTI_data_preparation.get_PPI_DTI_graph_intersection()
        PPI_graph = PPI_graph.subgraph(self.protein_list)

        # calculate dimensions of network
        self.num_proteins = len(PPI_graph.nodes())
        self.num_drugs = len(self.drug_list)

        print("Building index dict ...")
        self.protein_to_index_dict = {protein: index for index, protein in enumerate(self.protein_list)}
        print("Building edge list ...")
        forward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2])
                              for node1, node2 in list(PPI_graph.edges())]
        backward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2])
                               for node2, node1 in list(PPI_graph.edges())]
        self.edge_list = torch.tensor(np.transpose(np.array(forward_edges_list + backward_edges_list)),
                                      dtype=torch.long)
        self.num_PPI_features = 1

        # DDI data
        print("Loading DDI features ...")
        self.DDI_features = DTI_data_preparation.get_DDI_feature_list(self.drug_list)
        # add self-interactions for better performance
        self.DDI_features[np.identity(self.num_drugs)==1] = 1
        print(self.DDI_features.shape)

        # DTI data
        print("Loading DTI links ...")
        y_dti_data = DTI_data_preparation.get_DTIs(drug_list=self.drug_list, protein_list=self.protein_list)
        self.y_dti_data = y_dti_data.reshape((len(self.drug_list), len(self.protein_list)))
        print(self.y_dti_data.shape)

        print("Building feature matrix ...")
        self.train_prots = config.train_prots
        self.train_mask = np.zeros(self.num_proteins)
        self.train_mask[self.train_prots] = 1
        # self.test_prots = config.test_prots

        print('DDI_features', self.DDI_features[:10, :10])

        self.feature_matrix = np.zeros((self.num_drugs, self.num_proteins))
        epsilon = 0.00001
        for drug_index in tqdm(range(self.num_drugs)):
            drug_interactors = np.arange(len(self.drug_list))[self.DDI_features[drug_index, :] == 1]
            for drug_interactor in drug_interactors:
                # self.feature_matrix[drug_index, :] += self.train_mask * self.y_dti_data[drug_interactor, :]
                self.feature_matrix[drug_index, (self.train_mask * self.y_dti_data[drug_interactor, :]) == 1] = 1

            # print(list(self.feature_matrix[drug_index, :]))
            # print(list(self.y_dti_data[drug_index, :]))

            # self.feature_matrix[drug_index, :] = self.feature_matrix[drug_index, :] / (self.feature_matrix[drug_index, :].max() + epsilon)


        print('sum', ((self.feature_matrix[:, self.train_prots] - self.y_dti_data[:, self.train_prots])**2).sum())

        print(self.feature_matrix[:,self.train_prots][self.y_dti_data[:,self.train_prots]==0].shape)
        print(self.feature_matrix[:,self.train_prots][self.y_dti_data[:,self.train_prots]==0].sum())

        '''
        # Set data to true labels for sanity test
        for drug_index in tqdm(range(self.num_drugs)):
            self.feature_matrix[drug_index, :] = self.y_dti_data[drug_index, :]
        '''

        if not config.pretrain:
            print('Building protfunc data...')
            self.ProtFuncDataBuilder = ProteinFunctionDTIDataBuilder(config, num_proteins=config.num_proteins)

        print("Finished.\n")

        """
        test with open(file='graph_testing', mode='a') as f:
            print('\nTests')
            print(self.DDI_features.sum(), self.DDI_features.sum(axis=0), file=f)
            for drug_index in range(len(self.drug_list)):
                diff = (1 - self.train_mask) * self.feature_matrix[drug_index, :] - (1 - self.train_mask) * self.y_dti_data[drug_index, :]
                print(drug_index, self.feature_matrix[drug_index, :].sum(), self.y_dti_data[drug_index, :].sum(), np.linalg.norm(diff), file=f)
        """

    def get(self, indices):
        data_list = []

        protfunc_data = None
        if not self.config.pretrain:
            print('Loading protfunc data...')
            protfunc_data = self.ProtFuncDataBuilder.get(indices)

        # for index in tqdm(indices):
        for i in range(len(indices)):
            if (i + 1) % int(len(indices) / 10) == 0:
                print('Finished {} percent.'.format(str(int(i / len(indices) * 100))), end='\r')

            index = indices[i]
            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            # build protein mask
            protein_mask = np.zeros(self.num_proteins)
            protein_mask[protein_index] = 1
            protein_mask = torch.tensor(protein_mask, dtype=torch.bool)

            y = int(self.y_dti_data[drug_index, protein_index])

            feature_array = torch.tensor(self.feature_matrix[drug_index, :], dtype=torch.float).round().view(-1, 1)
            full_PPI_graph = Data(x=feature_array, edge_index=self.edge_list, y=y)
            full_PPI_graph.protein_mask = protein_mask
            if not self.config.pretrain:
                full_PPI_graph.protfunc_data = protfunc_data[i][0]

            # full_PPI_graph.__num_nodes__ = self.num_proteins

            data_list.append(full_PPI_graph)

        return data_list

    def __len__(self):
        return self.num_proteins * self.num_drugs

    def get_labels(self, indices):
        data_list = []

        # for index in tqdm(indices):
        for i in range(len(indices)):
            if (i + 1) % int(len(indices) / 10) == 0:
                print('Finished {} percent of labels.'.format(str(int(i / len(indices) * 100))), end='\r')

            index = indices[i]

            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins
            y = int(self.y_dti_data[drug_index, protein_index])
            data_list.append(y)

        return np.array(data_list)

class QuickProtFuncDTINetworkData:
    def __init__(self, config):
        self.config = config

        # write data first
        print("Preparing protein data ...")
        DTI_data_preparation.write_human_protein_list(min_score=config.PPI_min_score, mode=config.mode)
        DTI_data_preparation.write_human_prot_func_protein_list(mode=config.mode)

        print("Loading data ...")
        self.drug_list = np.array(DTI_data_preparation.get_drug_list(config.mode))
        print(len(self.drug_list), "drugs present")
        self.protein_list = np.array(DTI_data_preparation.get_human_prot_func_proteins())[:config.num_proteins]
        print(len(self.protein_list), "proteins present\n")

        # PPI data
        print("Loading PPI graph ...")
        self.PPI_graph = get_PPI_graph(min_score=config.PPI_min_score)
        self.PPI_graph = self.PPI_graph.subgraph(self.protein_list)

        # calculate dimensions of network
        self.num_proteins = len(self.protein_list)
        self.num_drugs = len(self.drug_list)


    def build_data(self, config):

        print('PPI_graph nodes/edges:', len(self.PPI_graph.nodes()), len(self.PPI_graph.edges()))

        print("Building index dict ...")
        self.protein_to_index_dict = {protein: index for index, protein in enumerate(self.protein_list)}
        print("Building edge list ...")
        forward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2])
                              for node1, node2 in list(self.PPI_graph.edges())]
        backward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2])
                               for node2, node1 in list(self.PPI_graph.edges())]
        self.edge_list = torch.tensor(np.transpose(np.array(forward_edges_list + backward_edges_list)),
                                      dtype=torch.long)
        self.num_PPI_features = 1

        print('Building edge feature attributes ...')
        forward_edge_feature_list = [1-self.PPI_graph[node1][node2]['score']/1000 for node1, node2 in list(self.PPI_graph.edges())]
        backward_edge_feature_list = [1-self.PPI_graph[node1][node2]['score']/1000 for node2, node1 in list(self.PPI_graph.edges())]
        self.edge_attr = torch.tensor(forward_edge_feature_list + backward_edge_feature_list, dtype=torch.float)# .view(-1,1)
        # self.edge_attr = torch.ones((self.edge_list.size(1),1), dtype=torch.float)


        # DDI data
        # print("Loading DDI features ...")
        # self.DDI_features = DTI_data_preparation.get_DDI_feature_list(self.drug_list)
        # add self-interactions for better performance
        # self.DDI_features[np.identity(self.num_drugs)==1] = 1
        # print(self.DDI_features.shape)

        # SemSim data
        # print("Loading semantic similarity data ...")
        # self.semsim_feature_matrix = DTI_data_preparation.get_side_effect_similarity_feature_list(self.drug_list)
        # print(self.semsim_feature_matrix[np.identity(self.num_drugs)==1])
        # print('semsim.shape', self.semsim_feature_matrix.shape)

        # DTI data
        print("Loading DTI links ...")
        y_dti_data = DTI_data_preparation.get_DTIs(drug_list=self.drug_list, protein_list=self.protein_list,
                                                   mode=config.mode)
        self.y_dti_data = y_dti_data.reshape((len(self.drug_list), len(self.protein_list)))
        print(self.y_dti_data.shape)

        print("Building feature matrix ...")
        self.train_prots = config.train_prots
        self.train_mask = np.zeros(self.num_proteins)
        self.train_mask[self.train_prots] = 1
        # self.test_prots = config.test_prots

        self.feature_matrix = np.zeros((self.num_drugs, self.num_proteins))
        epsilon = 0.00001

        '''
        for drug_index in tqdm(range(self.num_drugs)):
            drug_interactors = np.arange(len(self.drug_list))[self.DDI_features[drug_index, :] == 1]
            for drug_interactor in drug_interactors:
                # self.feature_matrix[drug_index, :] += self.train_mask * self.y_dti_data[drug_interactor, :]
                self.feature_matrix[drug_index, (self.train_mask * self.y_dti_data[drug_interactor, :]) == 1] = 1

            # print(list(self.feature_matrix[drug_index, :]))
            # print(list(self.y_dti_data[drug_index, :]))

            # self.feature_matrix[drug_index, :] = self.feature_matrix[drug_index, :] / (self.feature_matrix[drug_index, :].max() + epsilon)
        '''

        '''
        print('Building semsim PPI node features')
        for drug_index in tqdm(range(self.num_drugs)):
            for drug_interactor, drug_factor in enumerate(self.semsim_feature_matrix[drug_index, :]):
                self.feature_matrix[drug_index, :] += drug_factor * (self.train_mask * self.y_dti_data[drug_interactor, :])

            self.feature_matrix[drug_index,:] = self.feature_matrix[drug_index,:]/(self.feature_matrix.max() + epsilon)

        print('semsim_PPI.max', self.feature_matrix.max(), self.feature_matrix.mean())
        '''

        '''
        self.jaccard_sim_feature_matrix = DTI_data_preparation.get_jaccard_side_effect_similarity_feature_list(self.drug_list)
        print('Building jaccard similarity PPI node features ...')
        for drug_index in tqdm(range(self.num_drugs)):
            for drug_interactor, drug_factor in enumerate(self.jaccard_sim_feature_matrix[drug_index, :]):
                self.feature_matrix[drug_index, :] += drug_factor * (
                            self.train_mask * self.y_dti_data[drug_interactor, :])

            self.feature_matrix[drug_index, :] = self.feature_matrix[drug_index, :] / (self.feature_matrix.max() + epsilon)

        print('jacc_sim_PPI.max', self.feature_matrix.max(), self.feature_matrix.mean())
        '''

        '''
        histogram = np.zeros((24))
        for sim in self.feature_matrix.flatten():
            histogram[int(sim*8)] += 1
        print('histogram', histogram)
        '''

        # self.feature_matrix = self.feature_matrix/self.feature_matrix.max()

        # un-comment for DL2vec features
        self.drug_features = DTI_data_preparation.get_DL2vec_features(self.drug_list)
        self.protein_features = DTI_data_preparation.get_DL2vec_features(self.protein_list)
        self.num_PPI_features = self.drug_features.shape[1]*2 + 10

        print('feature shape', self.drug_features.shape, self.protein_features.shape)

        self.prot_func_features = DTI_data_preparation.get_protein_function_embeddings(protein_list=self.protein_list)

        self.node_degree_protein_feature = torch.tensor(DTI_data_preparation.get_protein_function_embeddings(protein_list=self.protein_list))


        '''
        # Set data to true labels for sanity test
        for drug_index in tqdm(range(self.num_drugs)):
            self.feature_matrix[drug_index, :] = self.y_dti_data[drug_index, :]
        '''

        if not config.pretrain:
            print('Building protfunc data...')
            self.ProtFuncDataBuilder = ProteinFunctionDTIDataBuilder(config, num_proteins=config.num_proteins)

        print("Finished.\n")

        """
        test with open(file='graph_testing', mode='a') as f:
            print('\nTests')
            print(self.DDI_features.sum(), self.DDI_features.sum(axis=0), file=f)
            for drug_index in range(len(self.drug_list)):
                diff = (1 - self.train_mask) * self.feature_matrix[drug_index, :] - (1 - self.train_mask) * self.y_dti_data[drug_index, :]
                print(drug_index, self.feature_matrix[drug_index, :].sum(), self.y_dti_data[drug_index, :].sum(), np.linalg.norm(diff), file=f)
        """

    def get(self):
        data_list = []

        indices = list(range(self.num_drugs))
        protfunc_data = None
        if not self.config.pretrain:
            print('Loading protfunc data...')
            protfunc_data = self.ProtFuncDataBuilder.get(indices)

        # for index in tqdm(indices):
        for drug_index in indices:
            # build protein mask

            y = torch.tensor(self.y_dti_data[drug_index, :]).view(-1)

            # feature_array = torch.tensor(self.feature_matrix[drug_index, :], dtype=torch.float).view(-1, 1)
            # feature_array = torch.tensor(self.y_dti_data[drug_index, :], dtype=torch.float).view(-1,1)

            # uncomment for DL2vec
            drug_feature = np.vstack([self.drug_features[drug_index, :]]*self.num_proteins)
            drug_feature = torch.tensor(drug_feature)

            # Dl2vec
            # protein_feature = self.protein_features
            protein_feature = self.prot_func_features

            # input node degree
            degree_feature = self.node_degree_protein_feature

            # uncomment for DL2vec
            feature_array = torch.stack([degree_feature, drug_feature, protein_feature], dim=1)
            # feature_array = torch.tensor(feature_array, dtype=torch.float)


            full_PPI_graph = Data(x=feature_array,
                                  edge_index=self.edge_list,
                                  edge_attr=self.edge_attr,
                                  y=y)

            # full_PPI_graph.__num_nodes__ = self.num_proteins

            data_list.append(full_PPI_graph)

        return data_list

    def __len__(self):
        return self.num_drugs

class DTIGraphDataset(Dataset):
    def __init__(self, data_list):
        super(DTIGraphDataset, self).__init__('/data/torch_raw/')
        # self.data, self.slices = self.collate(data_list)
        self.data_list = data_list

    def get(self, idx):
        return self.data_list[idx]

    def __len__(self):
        return len(self.data_list)

    def _download(self):
        pass

    def _process(self):
        pass

'''
def train(model, optimizer, loader, device):
    model.train()
    total_loss = 0
    for data in loader:
        optimizer.zero_grad()
        data = data.to(device)
        out = model(data)
        loss = F.nll_loss(out, data.y.view(-1))
        loss.backward()
        total_loss += loss.item() * loader.batch_size
        optimizer.step()
    return total_loss / len(loader.dataset)
'''

def train(model, device, train_loader, optimizer, epoch, weight_dict={0:1., 1:1.}):
    print('Training on {} samples...'.format(len(train_loader.dataset)))
    sys.stdout.flush()
    model.train()
    return_loss = 0
    for batch_idx, data in enumerate(train_loader):
        optimizer.zero_grad()
        output = model(data)
        # print('max/min:', output.max(), output.sigmoid().max(), output.min(), output.sigmoid().min())
        y = torch.Tensor([graph_data.y for graph_data in data]).float().to(output.device)

        weight_vec = torch.ones([1]) * weight_dict[1]

        loss = nn.BCEWithLogitsLoss(pos_weight=weight_vec.to(output.device))(output, y.view(-1, 1))
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

def quick_train(model, device, train_loader, optimizer, epoch, neg_to_pos_ratio, train_mask):
    print('Training on {} samples...'.format(len(train_loader.dataset)))
    sys.stdout.flush()
    model.train()
    return_loss = 0
    for batch_idx, data in enumerate(train_loader):
        optimizer.zero_grad()

        output = model(data)
        # print('max/min:', output.max(), output.sigmoid().max(), output.min(), output.sigmoid().min())

        y = torch.Tensor(np.array([graph_data.y.numpy() for graph_data in data])).float().to(output.device)

        # print('y.size', y[:, train_mask].size())
        # print('output.size', output[:, train_mask].size())

        pos_weights = torch.Tensor([neg_to_pos_ratio]) # * 0.75

        # print('check', output.min(), output.max(), y.min(), y.max())

        loss = nn.BCEWithLogitsLoss(pos_weight=pos_weights.to(device))(input=output[:, train_mask==1].view(-1, 1), target=y[:, train_mask==1].view(-1, 1),)
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
        for data in loader:
            # data = data.to(device)
            output = model(data).sigmoid()
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            y = torch.Tensor([graph_data.y for graph_data in data])
            total_labels = torch.cat((total_labels, y.view(-1, 1).float().cpu()), 0)
    return total_labels.round().numpy().flatten(),np.around(total_preds.numpy()).flatten()

def quick_predicting(model, device, loader):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    print('Make prediction for {} samples...'.format(len(loader.dataset)))
    with torch.no_grad():
        for data in loader:
            # data = data.to(device)
            output = model(data).sigmoid()
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            y = torch.Tensor(np.array([graph_data.y.numpy() for graph_data in data]))
            total_labels = torch.cat((total_labels.view(-1,
                                                        1), y.view(-1, 1).float().cpu()), 0)

    print('total_preds.max/min', total_labels.max(), total_labels.min())
    return total_labels.round().numpy().flatten(), np.around(total_preds.numpy()).flatten()


def rmse(y, f):
    rmse = sqrt(((y - f)**2).mean(axis=0))
    return rmse
def mse(y, f):
    mse = ((y - f)**2).mean(axis=0)
    return mse
def pearson(y, f):
    rp = np.corrcoef(y, f)[0,1]
    return rp
def spearman(y, f):
    rs = stats.spearmanr(y, f)[0]
    return rs
def ci(y, f):
    ind = np.argsort(y)
    y = y[ind]
    f = f[ind]
    i = len(y)-1
    j = i-1
    z = 0.0
    S = 0.0
    while i > 0:
        while j >= 0:
            if y[i] > y[j]:
                z = z+1
                u = f[i] - f[j]
                if u > 0:
                    S = S + 1
                elif u == 0:
                    S = S + 0.5
            j = j - 1
        i = i - 1
        j = i-1
    ci = S/z
    return ci
def cross_entropy(p,q):
    return -sum([p[i] * log2(q[i]) for i in range(len(p))])
'''
def eval_acc(model, loader, device):
    model.eval()

    correct = 0
    for data in loader:
        data = data.to(device)
        with torch.no_grad():
            pred = model(data).max(1)[1]
        correct += pred.eq(data.y.view(-1)).sum().item()
    return correct / len(loader.dataset)


def eval_loss(model, loader, device):
    model.eval()

    loss = 0
    for data in loader:
        data = data.to(device)
        with torch.no_grad():
            out = model(data)
        loss += F.nll_loss(out, data.y.view(-1), reduction='sum').item()
    return loss / len(loader.dataset)
'''
