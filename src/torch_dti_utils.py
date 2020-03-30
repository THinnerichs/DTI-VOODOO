import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

from torch_geometric.data import Dataset, Data

from tqdm import tqdm

import DTI_data_preparation



class FullNetworkDataset(Dataset):
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

    def __init__(self, num_proteins=None, transform=None, pre_transform=None):
        # super(FullNetworkDataset, self).__init__(root, transform, pre_transform)

        print("Loading data ...")
        self.drug_list = np.array(DTI_data_preparation.get_drug_list())
        print(len(self.drug_list), "drugs present")
        self.protein_list = np.array(DTI_data_preparation.get_human_proteins())[:num_proteins]
        print(len(self.protein_list), "proteins present\n")

        # PPI data
        print("Loading PPI graph ...")
        PPI_graph = DTI_data_preparation.get_PPI_DTI_graph_intersection()
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

    '''
    @property
    def raw_file_names(self):
        print("raw_file_names")
        return []

    @property
    def processed_file_names(self):
        print("processed_file_names")
        return ['../pytorch_data/PPI_network.dataset']
    '''

    def transform(self, data):
        return data

    def set_graph_train_mask(self, indizes):
        self.full_PPI_graph_Data.train_idx = torch.tensor(indizes, dtype=torch.long)

    def set_graph_test_mask(self, indizes):
        self.full_PPI_graph_Data.test_idx = torch.tensor(indizes, dtype=torch.long)

    def get_protein_to_index_dict(self):
        return self.protein_to_index_dict

    def get_protein_list(self):
        return self.protein_list

    def get_drug_list(self):
        return self.drug_list

    def get(self, indices):
        data_list = []

        for index in tqdm(indices):
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
            full_PPI_graph.__num_nodes__ = self.num_proteins

            data_list.append(full_PPI_graph)

        return data_list

    def __len__(self):
        return self.num_proteins * self.num_drugs

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

def train(model, device, train_loader, optimizer, epoch):
    print('Training on {} samples...'.format(len(train_loader.dataset)))
    model.train()
    for batch_idx, data in enumerate(train_loader):
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data)
        print('output.size()', output.size(), data.y.view(-1, 1).size())
        loss = nn.MSELoss()(output, data.y.view(-1, 1).float().to(device))
        loss.backward()
        optimizer.step()
        if batch_idx % 1 == 0:
            print('Train epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(epoch,
                                                                           batch_idx * output.size(0),
                                                                           len(train_loader.dataset),
                                                                           100. * batch_idx / len(train_loader),
                                                                           loss.item()))

def predicting(model, device, loader):
    model.eval()
    total_preds = torch.Tensor()
    total_labels = torch.Tensor()
    print('Make prediction for {} samples...'.format(len(loader.dataset)))
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            output = model(data)
            total_preds = torch.cat((total_preds, output.cpu()), 0)
            total_labels = torch.cat((total_labels, data.y.view(-1, 1).cpu()), 0)
    return total_labels.numpy().flatten(),total_preds.numpy().flatten()


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
