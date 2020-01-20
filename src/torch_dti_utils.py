import numpy as np

import torch
import torch_geometric
from torch_geometric.data import Dataset
from torch_geometric.data import Data

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

    def __init__(self, root, transform=None, pre_transform=None):
        super(FullNetworkDataset, self).__init__(root, transform, pre_transform)

        print("Loading data ...")
        drug_list = np.array(DTI_data_preparation.get_drug_list())
        print("Get protein list ...")
        protein_list = np.array(DTI_data_preparation.get_human_proteins())

        print(len(drug_list))
        print(len(protein_list))

        # PPI data
        print("Loading PPI graph ...")
        PPI_graph = DTI_data_preparation.get_PPI_DTI_graph_intersection()
        self.full_PPI_graph_Data = torch_geometric.utils.from_networkx(PPI_graph)
        print(self.full_PPI_graph_Data)

        # DDI data
        print("Loading DDI features ...")
        DDI_features = DTI_data_preparation.get_DDI_feature_list(drug_list)
        print(DDI_features.shape)

        # DTI data
        print("Loading DTI links ...")
        y_dti_data = DTI_data_preparation.get_DTIs(drug_list=drug_list, protein_list=protein_list)
        y_dti_data = y_dti_data.reshape((len(protein_list), len(drug_list)))
        y_dti_data = np.transpose(y_dti_data)
        print(y_dti_data.shape)

        # calculate dimensions of network
        self.num_proteins = len(PPI_graph.nodes())
        self.num_drugs = len(drug_list)
        print("Finished.\n")

    '''
    @property
    def raw_file_names(self):
        return []

    @property
    def processed_file_names(self):
        return ['../pytorch_data/PPI_network.dataset']

    def download(self):
        pass

    def process(self):
        pass
    '''

    def set_graph_train_mask(self, indizes):
        self.full_PPI_graph_Data.train_idx = torch.tensor(indizes, dtype=torch.long)

    def set_graph_test_mask(self, indizes):
        self.full_PPI_graph_Data.test_idx = torch.tensor(indizes, dtype=torch.long)

    def __getitem__(self, index):
        print("Index: {}".format(index))


    def __len__(self):
        return self.num_proteins * self.num_drugs

