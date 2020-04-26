import sys

import torch
import torch.nn.functional as F
import torch_geometric
import torch_geometric.nn as nn


class TemplateSimpleNet(torch.nn.Module):
    def __init__(self, config, num_drugs, num_prots, num_features, conv_method, GCN_num_outchannels=128, dropout=0.2):
        super(TemplateSimpleNet, self).__init__()
        self.batch_size = config.batch_size

        self.num_drugs = num_drugs
        self.num_prots = num_prots

        # DDI feature layers
        self.fc1 = torch.nn.Linear(num_drugs + GCN_num_outchannels, 256)
        self.fc2 = torch.nn.Linear(256,64)
        self.fc3 = torch.nn.Linear(64,16)
        self.fc4 = torch.nn.Linear(16,1)

        # mask feature

        # GCN layers
        if 'GCNConv' in conv_method:
            self.conv1 = nn.GCNConv(num_features, num_features*4, cached=False)
            self.conv2 = nn.GCNConv(num_features*4, num_features*16, cached=False)
            self.conv3 = nn.GCNConv(num_features*16, num_features*128, cached=False)
        elif 'ChebConv' in conv_method:
            self.conv1 = nn.ChebConv(num_features, num_features*4, 3)
            self.conv2 = nn.ChebConv(num_features*4, num_features*16, 3)
            self.conv3 = nn.ChebConv(num_features*16, num_features*128, 3)
        elif 'SAGEConv' in conv_method:
            self.conv1 = nn.SAGEConv(num_features, num_features*8)
            self.conv2 = nn.SAGEConv(num_features*4, num_features*64)
            self.conv3 = nn.SAGEConv(num_features*16, num_features*128)
        elif 'GraphConv' in conv_method:
            self.conv1 = nn.GraphConv(num_features, num_features*8)
            self.conv2 = nn.GraphConv(num_features*4, num_features*64)
            self.conv3 = nn.GraphConv(num_features*16, num_features*128)
        elif 'GATConv' in conv_method:
            self.conv1 = nn.GATConv(num_features, num_features*4, heads=5)
            self.conv2 = nn.GATConv(num_features*4, num_features*16, heads=5)
            self.conv3 = nn.GATConv(num_features*16, num_features*128, heads=5)
        elif 'TAGConv' in conv_method:
            self.conv1 = nn.TAGConv(num_features, num_features*8)
            self.conv2 = nn.TAGConv(num_features*4, num_features*64)
            self.conv3 = nn.TAGConv(num_features*16, num_features*128)
        elif 'ARMAConv' in conv_method:
            self.conv1 = nn.ARMAConv(num_features, num_features*8)
            self.conv2 = nn.ARMAConv(num_features*4, num_features*64)
            self.conv3 = nn.ARMAConv(num_features*16, num_features*128)
        elif 'SGConv' in conv_method:
            self.conv1 = nn.SGConv(num_features, num_features*8)
            self.conv2 = nn.SGConv(num_features*4, num_features*64)
            self.conv3 = nn.SGConv(num_features*16, num_features*128)
        elif 'FeaStConv' in conv_method:
            self.conv1 = nn.FeaStConv(num_features, num_features*4, heads=5)
            self.conv2 = nn.FeaStConv(num_features*4, num_features*16, heads=5)
            self.conv3 = nn.FeaStConv(num_features*16, num_features*128, heads=5)
        else:
            print("No valid model selected.")
            sys.stdout.flush()
            raise ValueError


        self.fc_g1 = torch.nn.Linear(num_features*128, 512)
        self.fc_g2 = torch.nn.Linear(512, GCN_num_outchannels)

        self.relu = torch.nn.ReLU()
        self.dropout = torch.nn.Dropout(dropout)

    def forward(self, PPI_data_object):
        # DDI_feature = PPI_data_object.DDI_features
        protein_mask = PPI_data_object.protein_mask
        PPI_x, PPI_edge_index, PPI_batch = PPI_data_object.x, PPI_data_object.edge_index, PPI_data_object.batch

        batch_size = self.batch_size


        # print(DDI_feature.shape)
        # print('protein_mask.size()', protein_mask.size())

        # PPI graph network
        PPI_x = self.conv1(PPI_x, PPI_edge_index)
        PPI_x = F.relu(PPI_x)
        # PPI_x = F.dropout(PPI_x, training=self.training)
        PPI_x = self.conv2(PPI_x, PPI_edge_index)
        PPI_x = F.relu(PPI_x)
        PPI_x = F.relu(self.conv3(PPI_x, PPI_edge_index))

        PPI_x = PPI_x.view((batch_size, self.num_prots, PPI_x.shape[-1]))

        # PPI_x = torch_geometric.nn.global_max_pool(PPI_x, PPI_batch.view((batch_size, -1)))

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
        # DDI_x = torch.cat((PPI_x, DDI_feature), 1)
        DDI_x = self.fc1(PPI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc2(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc3(DDI_x)
        DDI_x = F.relu(DDI_x)
        DDI_x = self.fc4(DDI_x)
        # DDI_x = torch.sigmoid(DDI_x)
        return DDI_x

class ResTemplateNet(torch.nn.Module):
    def __init__(self, config, num_drugs, num_prots, num_features, conv_method, out_channels=64, dropout=0.2):
        super(ResTemplateNet, self).__init__()
        self.batch_size = config.batch_size

        self.num_drugs = num_drugs
        self.num_prots = num_prots

        self.out_channels = out_channels

        if 'ResGCNConv' in conv_method:
            self.conv1 = nn.GCNConv(num_features, out_channels, cached=False)
            self.conv2 = nn.GCNConv(out_channels, out_channels, cached=False)
            self.conv3 = nn.GCNConv(out_channels, out_channels, cached=False)
        elif 'ResChebConv' in conv_method:
            self.conv1 = nn.ChebConv(num_features, out_channels, 3)
            self.conv2 = nn.ChebConv(out_channels, out_channels, 3)
            self.conv3 = nn.ChebConv(out_channels, out_channels, 3)
        elif 'ResSAGEConv' in conv_method:
            self.conv1 = nn.SAGEConv(num_features, out_channels)
            self.conv2 = nn.SAGEConv(out_channels, out_channels)
            self.conv3 = nn.SAGEConv(out_channels, out_channels)
        elif 'ResGraphConv' in conv_method:
            self.conv1 = nn.GraphConv(num_features, out_channels)
            self.conv2 = nn.GraphConv(out_channels, out_channels)
            self.conv3 = nn.GraphConv(out_channels, out_channels)
        elif 'ResGATConv' in conv_method:
            self.conv1 = nn.GATConv(num_features, out_channels, heads=5)
            self.conv2 = nn.GATConv(out_channels, out_channels, heads=5)
            self.conv3 = nn.GATConv(out_channels, out_channels, heads=5)
        elif 'ResTAGConv' in conv_method:
            self.conv1 = nn.TAGConv(num_features, out_channels)
            self.conv2 = nn.TAGConv(out_channels, out_channels)
            self.conv3 = nn.TAGConv(out_channels, out_channels)
        elif 'ResARMAConv' in conv_method:
            self.conv1 = nn.ARMAConv(num_features, out_channels)
            self.conv2 = nn.ARMAConv(out_channels, out_channels)
            self.conv3 = nn.ARMAConv(out_channels, out_channels)
        elif 'ResSGConv' in conv_method:
            self.conv1 = nn.SGConv(num_features, out_channels)
            self.conv2 = nn.SGConv(out_channels, out_channels)
            self.conv3 = nn.SGConv(out_channels, out_channels)
        elif 'ResFeaStConv' in conv_method:
            self.conv1 = nn.FeaStConv(num_features, out_channels, heads=5)
            self.conv2 = nn.FeaStConv(out_channels, out_channels, heads=5)
            self.conv3 = nn.FeaStConv(out_channels, out_channels, heads=5)
        else:
            print("No valid model selected.")
            sys.stdout.flush()
            raise ValueError


        self.lin1 = torch.nn.Linear(out_channels + self.num_drugs, 256)
        self.lin2 = torch.nn.Linear(256, 64)
        self.lin3 = torch.nn.Linear(64, 16)
        self.lin4 = torch.nn.Linear(16, 1)


    def forward(self, data):
        # DDI_feature = data.DDI_features
        protein_mask = data.protein_mask
        x, edge_index, batch = data.x, data.edge_index, data.batch
        batch_size = self.batch_size

        gmp = torch_geometric.nn.global_max_pool
        gap = torch_geometric.nn.global_add_pool

        protein_mask = protein_mask.view((batch_size, 1, -1)).float()

        x = F.relu(self.conv1(x, edge_index))
        # x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
        x1 = x.view((batch_size, -1))
        # x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = F.relu(self.conv2(x, edge_index))
        # x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)
        x2 = x.view((batch_size, -1))
        # x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = F.relu(self.conv3(x, edge_index))
        # x, edge_index, _, batch, _, _ = self.pool1(x, edge_  index, None, batch)
        x3 = x.view((batch_size, -1))
        # x3 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)

        x = x1 + x2 + x3

        x = x.view((batch_size, self.num_prots, self.out_channels))
        x = torch.bmm(protein_mask, x)
        x = x.view((batch_size, -1))

        # x = torch.cat((x, DDI_feature), 1)

        x = F.relu(self.lin1(x))
        # x = F.dropout(x, p=0.5, training=self.training)
        x = F.relu(self.lin2(x))
        x = F.relu(self.lin3(x))
        x = self.lin4(x)

        return x

