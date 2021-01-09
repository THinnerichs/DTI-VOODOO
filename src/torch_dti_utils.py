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


        if config.yamanishi_test:
            print("Loading Yamanishi data ...")
            self.drug_list, self.protein_list, self.y_dti_data = DTI_data_preparation.get_yamanishi_data(self.drug_list, self.protein_list)

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

        # DTI data
        if not config.yamanishi_test:
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

        self.num_PPI_features = 200# +100


        self.node_degree_protein_feature = torch.tensor(DTI_data_preparation.get_protein_degree_percentile(protein_list=self.protein_list, n=100))



        print("Finished.\n")

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

        # for index in tqdm(indices):
        for drug_index in indices:
            # build protein mask

            y = torch.tensor(self.y_dti_data[drug_index, :]).view(-1)

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

def quick_predicting(model, device, loader, round=True):
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

    if round:
        return total_labels.round().numpy().flatten(), np.around(total_preds.numpy()).flatten()
    else:
        return total_labels.numpy().flatten(), total_preds.numpy().flatten()


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

