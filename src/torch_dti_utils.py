import numpy as np
from math import sqrt, log2
from scipy import stats

import gensim

import torch
import torch.nn as nn
import torch.nn.functional as F

from torch_geometric.data import Dataset, Data

from tqdm import tqdm
import sys
import pickle

import DTI_data_preparation
from PPI_utils import get_PPI_graph
import PhenomeNET_DL2vec_utils
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

        '''
        # DDI data
        print("Loading DDI features ...")
        self.DDI_features = DTI_data_preparation.get_DDI_feature_list(self.drug_list)
        # add self-interactions for better performance
        self.DDI_features[np.identity(self.num_drugs)==1] = 1
        print(self.DDI_features.shape)
        '''

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

        print("Loading data ...")
        self.drug_list = np.array(DTI_data_preparation.get_drug_list(config.mode))
        print(len(self.drug_list), "drugs present")
        # self.protein_list = np.array(DTI_data_preparation.get_human_prot_func_proteins())[:config.num_proteins]
        # get protein lists for each ontology
        uberon_protein_list = PhenomeNET_DL2vec_utils.get_PhenomeNET_protein_list(mode='uberon')
        GO_protein_list = PhenomeNET_DL2vec_utils.get_PhenomeNET_protein_list(mode='GO')
        MP_protein_list = PhenomeNET_DL2vec_utils.get_PhenomeNET_protein_list(mode='MP')

        dti_graph = DTI_data_preparation.get_human_DTI_graph(mode=config.mode)
        self.PPI_graph = get_PPI_graph(min_score=config.PPI_min_score)
        self.protein_list = np.array(list(set(self.PPI_graph.nodes()) & set(dti_graph.nodes()) & (set(uberon_protein_list) | set(GO_protein_list) | set(MP_protein_list))))
        print(len(self.protein_list), "proteins present.\n")
        # self.protein_list = np.array(DTI_data_preparation.get_human_PhenomeNET_proteins())#[:config.num_proteins])

        if config.include_mol_features:
            print("Building molecular features")
            # build drug features
            if config.drug_mode == 'trfm' or config.drug_mode == 'rnn':
                drug_filename = '../models/drug_representation/' + config.drug_mode + '.npy'
                self.drug_mol_encodings = torch.from_numpy(np.load(drug_filename))
            else:
                print("No valid mode selected for drug to SMILES encoding.")
                raise ValueError

            # build molecular protein features
            filename = '../models/protein_representation/results/prot_to_encoding_dict'
            with open(file=filename + '.pkl', mode='rb') as f:
                protein_to_feature_dict = pickle.load(f)

            self.protein_list = np.array(list(set(self.protein_list) & set(protein_to_feature_dict.keys())))
            self.protein_mol_encodings = torch.Tensor([protein_to_feature_dict[protein] for protein in self.protein_list])

            print(len(self.protein_list), "proteins present with mol_pred_intersection.\n")

        # PPI data
        print("Loading PPI graph ...")
        self.PPI_graph = self.PPI_graph.subgraph(self.protein_list)

        # calculate dimensions of network
        self.num_proteins = len(self.protein_list)
        self.num_drugs = len(self.drug_list)

        config.num_drugs = self.num_drugs
        config.num_proteins = self.num_proteins


    def build_data(self, config):

        print('PPI_graph nodes/edges:', len(self.PPI_graph.nodes()), len(self.PPI_graph.edges()))

        print("Building index dict ...")
        self.protein_to_index_dict = {protein: index for index, protein in enumerate(self.protein_list)}
        print("Building edge list ...")
        forward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2]) for node1, node2 in list(self.PPI_graph.edges())]
        backward_edges_list = [(self.protein_to_index_dict[node1], self.protein_to_index_dict[node2]) for node2, node1 in list(self.PPI_graph.edges())]

        self.edge_list = torch.tensor(np.transpose(np.array(forward_edges_list + backward_edges_list)), dtype=torch.long)
        self.num_PPI_features = 1

        print('Building edge feature attributes ...')
        forward_edge_feature_list = [self.PPI_graph[node1][node2]['score']/1000 for node1, node2 in list(self.PPI_graph.edges())]
        backward_edge_feature_list = [self.PPI_graph[node1][node2]['score']/1000 for node2, node1 in list(self.PPI_graph.edges())]
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
        y_dti_data = DTI_data_preparation.get_DTIs(drug_list=self.drug_list, protein_list=self.protein_list, mode=config.mode)
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

        # self.feature_matrix =    self.feature_matrix/self.feature_matrix.max()

        self.drug_features = DTI_data_preparation.get_DL2vec_features(self.drug_list)
        self.protein_features = DTI_data_preparation.get_DL2vec_features(self.protein_list)

        self.num_PPI_features = 200# +100

        # print('feature shape', self.drug_features.shape, self.protein_features.shape)

        # self.prot_func_features = DTI_data_preparation.get_protein_function_embeddings(protein_list=self.protein_list)

        self.node_degree_protein_feature = torch.tensor(DTI_data_preparation.get_protein_degree_percentile(protein_list=self.protein_list, n=100))

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

        DL2vec_path_prefix = '../data/PhenomeNET_data/'

        drug_model_filename = DL2vec_path_prefix + 'drug_embedding_model'
        uberon_model_filename = DL2vec_path_prefix + 'uberon_embedding_model'
        GO_model_filename = DL2vec_path_prefix + 'GO_embedding_model'
        MP_model_filename = DL2vec_path_prefix + 'MP_embedding_model'

        # load models
        drug_model = gensim.models.Word2Vec.load(drug_model_filename)
        uberon_model = gensim.models.Word2Vec.load(uberon_model_filename)
        GO_model = gensim.models.Word2Vec.load(GO_model_filename)
        MP_model = gensim.models.Word2Vec.load(MP_model_filename)

        # Build wordvector dicts
        drug_model = drug_model.wv
        uberon_model = uberon_model.wv
        GO_model = GO_model.wv
        MP_model = MP_model.wv

        drug_embeddings = []
        uberon_embeddings = []
        GO_embeddings = []
        MP_embeddings = []
        for protein in self.protein_list:
            # organism, protein_id = protein.strip().split('.')
            protein_id = protein

            if protein_id in uberon_model.vocab.keys():
                uberon_embeddings.append(uberon_model[protein_id])
            else:
                uberon_embeddings.append(torch.zeros((200)))
            if protein_id in GO_model.vocab.keys():
                GO_embeddings.append(GO_model[protein_id])
            else:
                GO_embeddings.append(torch.zeros((200)))
            if protein_id in MP_model.vocab.keys():
                MP_embeddings.append(MP_model[protein_id])
            else:
                GO_embeddings.append(torch.zeros((200)))

        for drug_id in self.drug_list:
            drug_embeddings.append((drug_model[drug_id]))

        self.drug_embeddings = torch.Tensor(drug_embeddings)
        self.uberon_embeddings = torch.Tensor(uberon_embeddings)
        self.GO_embeddings = torch.Tensor(GO_embeddings)
        self.MP_embeddings = torch.Tensor(MP_embeddings)

        self.protein_embeddings = torch.cat([self.uberon_embeddings, self.GO_embeddings, self.MP_embeddings], dim=1)

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
            # drug_feature = np.vstack([self.drug_features[drug_index, :]]*self.num_proteins)
            # drug_feature = torch.tensor(drug_feature)
            # drug_feature = torch.tensor(self.drug_embeddings[drug_index, :])

            # Dl2vec
            protein_feature = torch.tensor(self.protein_features)
            # protein_feature = self.prot_func_features

            # input node degree
            degree_feature = self.node_degree_protein_feature

            # feature_array = torch.cat([degree_feature, drug_feature, protein_feature], dim=1)
            # feature_array = torch.cat([drug_feature, protein_feature], dim=1)
            feature_array = protein_feature
            # feature_array = torch.tensor(degree_feature, dtype=torch.float)

            molecular_drug_feature = self.drug_mol_encodings[drug_index,:]


            full_PPI_graph = Data(x=self.protein_embeddings,
                                  edge_index=self.edge_list,
                                  edge_attr=self.edge_attr,
                                  y=y)

            full_PPI_graph.drug_feature = self.drug_embeddings[drug_index, :]
            full_PPI_graph.drug_mol_feature = molecular_drug_feature
            full_PPI_graph.protein_mol_feature = self.protein_mol_encodings

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


def BCELoss_ClassWeights(input, target, pos_weight):
    # input (n, d)
    # target (n, d)
    # class_weights (1, d)
    input = torch.clamp(input,min=1e-7,max=1-1e-7)
    target = target.view(-1,1)
    weighted_bce = - pos_weight*target * torch.log(input) - 1*(1 - target) * torch.log(1 - input)
    # weighted_bce = weighted_bce / (pos_weight+1)
    final_reduced_over_batch = weighted_bce.sum(axis=0)
    return final_reduced_over_batch

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
        return_loss += loss.item()
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

def quick_train(config, model, device, train_loader, optimizer, epoch, neg_to_pos_ratio, train_mask):
    print('Training on {} samples...'.format(len(train_loader.dataset)))
    sys.stdout.flush()
    model.train()
    return_loss = 0

    for batch_idx, data in enumerate(train_loader):
        optimizer.zero_grad()

        output = model(data)


        # print('max/min:', output.max(), output.sigmoid().max(), output.min(), output.sigmoid().min())

        y = torch.Tensor(np.array([graph_data.y.numpy() for graph_data in data])).float().to(output.device)

        '''
        help_mask = np.around(np.array(y.to('cpu')) * train_mask).astype(np.int)
        for i in range(help_mask.shape[0]):
            # determine number of positive samples per drug/graph
            num_choices = help_mask.sum(axis=1)[i]
            # choose num_choices indices from num_proteins samples (masked by train_mask) without replacement and set their entries to 1 in help mask
            indices = np.arange(help_mask.shape[1])[help_mask[i,:]==0]
            help_mask[i,np.random.choice(indices, num_choices, replace=False)] = 1
        '''


        # print('y.size', y[:, train_mask].size())
        # print('output.size', output[:, train_mask].size())

        pos_weights = torch.Tensor([neg_to_pos_ratio]) # * 0.75

        # print('check', output.min(), output.max(), y.min(), y.max())

        # my implementation of BCELoss
        output = torch.clamp(output, min=1e-7, max=1 - 1e-7)

        pos_weight = neg_to_pos_ratio
        neg_weight = 1
        loss = BCELoss_ClassWeights(input=output[:, train_mask==1].view(-1,1), target=y[:,train_mask==1].view(-1,1), pos_weight=pos_weight)
        loss = loss/(config.num_drugs*config.num_proteins)

        # loss = nn.BCEWithLogitsLoss(pos_weight=pos_weights.to(device))(input=output[:, train_mask==1].view(-1, 1), target=y[:, train_mask==1].view(-1, 1),)
        # loss = nn.BCELoss(reduction='mean')(input=output[help_mask==1].view(-1, 1), target=y[help_mask==1].view(-1, 1))
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
            output = model(data)#.sigmoid()
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
            output = model(data)#.sigmoid()
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            y = torch.Tensor(np.array([graph_data.y.numpy() for graph_data in data]))
            total_labels = torch.cat((total_labels.view(-1,1), y.view(-1, 1).float().cpu()), 0)

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

