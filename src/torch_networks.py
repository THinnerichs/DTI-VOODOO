import numpy as np

import torch
import torch.nn.functional as F
import torch_geometric
# import torch_geometric.nn as nn


class SimpleConvGCN(torch.nn.Module):
    def __init__(self, num_drugs, num_features, GCN_num_outchannels=32, embedding_layers_sizes = [32, 64], dropout=0.2):
        super(SimpleConvGCN, self).__init__()

        # DDI feature layers
        self.fc1 = torch.nn.Linear(num_drugs + GCN_num_outchannels, 64)
        self.fc2 = torch.nn.Linear(64,16)
        self.fc3 = torch.nn.Linear(16,1)

        # mask feature

        # GCN layers
        self.conv1 = torch_geometric.nn.GCNConv(num_features, num_features, cached=False)
        self.conv2 = torch_geometric.nn.GCNConv(num_features, num_features*2, cached=False)
        self.fc_g1 = torch.nn.Linear(num_features*4, 1028)
        self.fc_g2 = torch.nn.Linear(1028, GCN_num_outchannels)

        self.relu = torch.nn.ReLU()
        self.dropout = torch.nn.Dropout(dropout)

        # self.conv1.weight = torch.nn.Parameter(self.conv1.weight.byte())

        # print(self.conv1.weight)
        # print(type(self.conv1.weight))



    def forward(self, PPI_data_object):
        DDI_feature = PPI_data_object.DDI_features
        protein_mask = PPI_data_object.protein_mask
        PPI_x, PPI_edge_index, PPI_batch = PPI_data_object.x, PPI_data_object.edge_index, PPI_data_object.batch

        # PPI graph network
        PPI_x = self.conv1(PPI_x, PPI_edge_index)
        PPI_x = F.relu(PPI_x)
        # PPI_x = F.dropout(PPI_x, training=self.training)
        PPI_x = self.conv2(PPI_x, PPI_edge_index)
        PPI_x = F.relu(PPI_x)

        PPI_x = torch_geometric.nn.global_max_pool(PPI_x, PPI_batch)
        print("PPI_shape:", PPI_x.shape)

        # flatten
        PPI_x = self.fc_g1(PPI_x)
        PPI_x = self.relu(PPI_x)
        PPI_x = self.dropout(PPI_x)
        PPI_x = self.fc_g2(PPI_x)
        PPI_x = self.dropout(PPI_x)
        print("PPI_shape:", PPI_x.shape)

        print(self.fc1.weight)
        print(self.fc1.weight.shape)

        # DDI feature network

        DDI_x = self.fc1(DDI_feature)
        DDI_x = F.relu(DDI_x)
        DDI_x = F.relu(self.fc2(DDI_x))
        DDI_x = F.relu(self.fc3(DDI_x))
        print("DDI_shape:", DDI_x.shape)
        raise Exception
        return x
