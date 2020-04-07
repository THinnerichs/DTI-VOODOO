import numpy as np

import torch
import torch.nn.functional as F
import torch_geometric
import torch_geometric.nn as nn


class SimpleConvGCN(torch.nn.Module):
    def __init__(self, num_drugs, num_prots, num_features, GCN_num_outchannels=32, embedding_layers_sizes = [32, 64], dropout=0.2):
        super(SimpleConvGCN, self).__init__()

        self.num_drugs = num_drugs
        self.num_prots = num_prots

        # DDI feature layers
        self.fc1 = torch.nn.Linear(num_drugs + GCN_num_outchannels, 64)
        self.fc2 = torch.nn.Linear(64,16)
        self.fc3 = torch.nn.Linear(16,1)

        # mask feature

        # GCN layers
        self.conv1 = torch_geometric.nn.GCNConv(num_features, num_features, cached=False)
        self.conv2 = torch_geometric.nn.GCNConv(num_features, num_features*2, cached=False)
        self.fc_g1 = torch.nn.Linear(num_features*2, 1028)
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

        batch_size = DDI_feature.size(0)

        # print(DDI_feature.shape)
        # print('protein_mask.size()', protein_mask.size())


        # PPI graph network
        PPI_x = self.conv1(PPI_x, PPI_edge_index)
        PPI_x = F.relu(PPI_x)
        # PPI_x = F.dropout(PPI_x, training=self.training)
        PPI_x = self.conv2(PPI_x, PPI_edge_index)
        PPI_x = F.relu(PPI_x)

        PPI_x = PPI_x.view((batch_size, self.num_prots, PPI_x.shape[-1]))

        PPI_x = torch_geometric.nn.global_max_pool(PPI_x, PPI_batch.view((batch_size, -1)))

        protein_mask = protein_mask.view((batch_size, 1, -1)).float()

        # multiply for flattening
        PPI_x = torch.bmm(protein_mask, PPI_x)
        PPI_x = PPI_x.view((batch_size, -1))

        # flatten
        PPI_x = self.fc_g1(PPI_x)
        PPI_x = self.relu(PPI_x)
        PPI_x = self.dropout(PPI_x)
        PPI_x = self.fc_g2(PPI_x)
        PPI_x = self.dropout(PPI_x)

        # DDI feature network
        DDI_x = torch.cat((PPI_x, DDI_feature), 1)
        DDI_x = self.fc1(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc2(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc3(DDI_x)
        # DDI_x = torch.sigmoid(DDI_x)
        return DDI_x

class TopKPoolingSimpleGCN(torch.nn.Module):
    def __init__(self, num_drugs, num_prots, num_features, GCN_num_outchannels=32, embedding_layers_sizes = [32, 64], dropout=0.2):
        super(TopKPoolingSimpleGCN, self).__init__()

        self.num_drugs = num_drugs
        self.num_prots = num_prots

        # DDI feature layers
        self.fc1 = torch.nn.Linear(num_drugs + GCN_num_outchannels, 64)
        self.fc2 = torch.nn.Linear(64,16)
        self.fc3 = torch.nn.Linear(16,1)

        # mask feature

        # GCN layers
        self.conv1 = torch_geometric.nn.GCNConv(num_features, num_features, cached=False)
        self.conv2 = torch_geometric.nn.GCNConv(num_features, num_features*2, cached=False)
        self.pooling1 = torch_geometric.nn.TopKPooling(num_features, min_score=None)
        self.pooling2 = torch_geometric.nn.TopKPooling(num_features*2, min_score=None)
        self.fc_g1 = torch.nn.Linear(num_features*2, 1028)
        self.fc_g2 = torch.nn.Linear(1028, GCN_num_outchannels)

        self.relu = torch.nn.ReLU()
        self.dropout = torch.nn.Dropout(dropout)

    def forward(self, PPI_data_object):
        DDI_feature = PPI_data_object.DDI_features
        protein_mask = PPI_data_object.protein_mask
        PPI_x, PPI_edge_index, PPI_batch = PPI_data_object.x, PPI_data_object.edge_index, PPI_data_object.batch

        batch_size = DDI_feature.size(0)

        # PPI graph network
        PPI_out = self.conv1(PPI_x, PPI_edge_index)
        PPI_out = F.relu(PPI_out)
        # PPI_out = F.dropout(PPI_out, training=self.training)
        out, edge_index, _, batch, _, _ = self.pooling1(PPI_out, PPI_edge_index, None, PPI_batch, attn=PPI_x)
        PPI_out = self.conv2(PPI_out, PPI_edge_index)
        out, edge_index, _, batch, _, _ = self.pooling2(PPI_out, PPI_edge_index, None, PPI_batch, attn=PPI_x)

        PPI_out = F.relu(PPI_out)
        PPI_out = PPI_out.view((batch_size, self.num_prots, PPI_out.shape[-1]))

        protein_mask = protein_mask.view((batch_size, 1, -1)).float()

        # multiply for flattening
        PPI_out = torch.bmm(protein_mask, PPI_out)
        PPI_out = PPI_out.view((batch_size, -1))

        # flatten
        PPI_out = self.fc_g1(PPI_out)
        PPI_out = self.relu(PPI_out)
        PPI_out = self.dropout(PPI_out)
        PPI_out = self.fc_g2(PPI_out)
        PPI_out = self.dropout(PPI_out)

        # DDI feature network
        DDI_x = torch.cat((PPI_out, DDI_feature), 1)
        DDI_x = self.fc1(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc2(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc3(DDI_x)
        # DDI_x = torch.sigmoid(DDI_x)
        return DDI_x

class ResTopKGCN(torch.nn.Module):
    def __init__(self, num_drugs, num_prots, num_features, GCN_num_outchannels=32, embedding_layers_sizes = [32, 64], dropout=0.2):
        super(ResTopKGCN, self).__init__()
        self.num_drugs = num_drugs
        self.num_prots = num_prots

        self.conv1 = nn.GraphConv(num_features, 128)
        # self.pool1 = nn.TopKPooling(128, ratio=1.0)
        self.conv2 = nn.GraphConv(128, 128)
        # self.pool2 = nn.TopKPooling(128, min_score=None)
        self.conv3 = nn.GraphConv(128, 128)
        # self.pool3 = nn.TopKPooling(128, min_score=None)

        self.lin1 = torch.nn.Linear(128 + self.num_drugs, 128)
        self.lin2 = torch.nn.Linear(128, 64)
        self.lin3 = torch.nn.Linear(64, 1)


    def forward(self, data):
        DDI_feature = data.DDI_features
        protein_mask = data.protein_mask
        x, edge_index, batch = data.x, data.edge_index, data.batch
        batch_size = DDI_feature.size(0)

        gmp = torch_geometric.nn.global_max_pool
        gap = torch_geometric.nn.global_add_pool

        protein_mask = protein_mask.view((batch_size, 1, -1)).float()

        print('size check')
        x = F.relu(self.conv1(x, edge_index))
        # x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
        # x = x.view((batch_size, -1))
        x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = F.relu(self.conv2(x, edge_index))
        # x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
        # x2 = x.view((batch_size, -1))
        x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = F.relu(self.conv3(x, edge_index))
        # x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
        # x3 = x.view((batch_size, -1))
        x3 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = x1 + x2 + x3


        x = x.view((batch_size, self.num_prots, 128))
        x = torch.bmm(protein_mask, x)
        x = x.view((batch_size, -1))

        x = torch.cat((x, DDI_feature), 1)

        x = F.relu(self.lin1(x))
        x = F.dropout(x, p=0.5, training=self.training)
        x = F.relu(self.lin2(x))
        x = self.lin3(x)

        return x

class ChebConvNet(torch.nn.Module):
    def __init__(self, num_drugs, num_prots, num_features, GCN_num_outchannels=256, embedding_layers_sizes=[32, 64],
                 dropout=0.2):
        super(ChebConvNet, self).__init__()

        self.num_drugs = num_drugs
        self.num_prots = num_prots

        # DDI feature layers
        self.fc1 = torch.nn.Linear(num_drugs + GCN_num_outchannels, 128)
        self.fc2 = torch.nn.Linear(128, 32)
        self.fc3 = torch.nn.Linear(32, 8)
        self.fc4 = torch.nn.Linear(8, 1)

        # mask feature

        # GCN layers
        self.conv1 = torch_geometric.nn.ChebConv(num_features, num_features, cached=False)
        self.conv2 = torch_geometric.nn.ChebConv(num_features, num_features * 2, cached=False)
        self.fc_g1 = torch.nn.Linear(num_features * 2, 1028)
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

        batch_size = DDI_feature.size(0)

        # print(DDI_feature.shape)
        # print('protein_mask.size()', protein_mask.size())

        # PPI graph network
        PPI_x = self.conv1(PPI_x, PPI_edge_index)
        PPI_x = F.relu(PPI_x)
        # PPI_x = F.dropout(PPI_x, training=self.training)
        PPI_x = self.conv2(PPI_x, PPI_edge_index)
        PPI_x = F.relu(PPI_x)

        PPI_x = PPI_x.view((batch_size, self.num_prots, PPI_x.shape[-1]))

        PPI_x = torch_geometric.nn.global_max_pool(PPI_x, PPI_batch.view((batch_size, -1)))

        protein_mask = protein_mask.view((batch_size, 1, -1)).float()

        # multiply for flattening
        PPI_x = torch.bmm(protein_mask, PPI_x)
        PPI_x = PPI_x.view((batch_size, -1))

        # flatten
        PPI_x = self.fc_g1(PPI_x)
        PPI_x = self.relu(PPI_x)
        PPI_x = self.dropout(PPI_x)
        PPI_x = self.fc_g2(PPI_x)
        PPI_x = self.dropout(PPI_x)

        # DDI feature network
        DDI_x = torch.cat((PPI_x, DDI_feature), 1)
        DDI_x = self.fc1(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc2(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc3(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc4(DDI_x)
        # DDI_x = torch.sigmoid(DDI_x)
        return DDI_x


class SAGEConvNet(torch.nn.Module):
    def __init__(self):
        pass

