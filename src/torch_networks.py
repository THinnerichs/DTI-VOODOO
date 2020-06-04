import sys

import torch
import torch.nn.functional as F
import torch_geometric
import torch_geometric.nn as nn

from protein_function_utils import ProteinFunctionPredNet


class TemplateSimpleNet(torch.nn.Module):
    def __init__(self, config, num_drugs, num_prots, num_features, conv_method, GCN_num_outchannels=64, dropout=0.2):
        super(TemplateSimpleNet, self).__init__()

        self.num_drugs = num_drugs
        self.num_prots = num_prots


        # mask feature

        # GCN laye4s
        if 'GCNConv' in conv_method:
            self.conv1 = nn.GCNConv(num_features, num_features*1, cached=False)
            # self.conv2 = nn.GCNConv(num_features*8, num_features*32, cached=False)
            # self.conv3 = nn.GCNConv(num_features*16, num_features*128, cached=False)
        elif 'ChebConv' in conv_method:
            self.conv1 = nn.ChebConv(num_features, num_features*8, 3)
            self.conv2 = nn.ChebConv(num_features*8, num_features*32, 3)
            # self.conv3 = nn.ChebConv(num_features*16, num_features*128, 3)
        elif 'SAGEConv' in conv_method:
            self.conv1 = nn.SAGEConv(num_features, num_features*8)
            self.conv2 = nn.SAGEConv(num_features*8, num_features*32)
            # self.conv3 = nn.SAGEConv(num_features*16, num_features*128)
        elif 'GraphConv' in conv_method:
            self.conv1 = nn.GraphConv(num_features, num_features*4)
            self.conv2 = nn.GraphConv(num_features*4, num_features*16)
            self.conv3 = nn.GraphConv(num_features*16, num_features*128)
        elif 'GATConv' in conv_method:
            self.conv1 = nn.GATConv(num_features, num_features*8, heads=5)
            self.conv2 = nn.GATConv(num_features*8, num_features*32, heads=5)
            # self.conv3 = nn.GATConv(num_features*16, num_features*128, heads=5)
        elif 'TAGConv' in conv_method:
            self.conv1 = nn.TAGConv(num_features, num_features*4)
            self.conv2 = nn.TAGConv(num_features*4, num_features*16)
            self.conv3 = nn.TAGConv(num_features*16, num_features*128)
        elif 'ARMAConv' in conv_method:
            self.conv1 = nn.ARMAConv(num_features, num_features*4)
            self.conv2 = nn.ARMAConv(num_features*4, num_features*16)
            self.conv3 = nn.ARMAConv(num_features*16, num_features*128)
        elif 'SGConv' in conv_method:
            self.conv1 = nn.SGConv(num_features, num_features*4)
            self.conv2 = nn.SGConv(num_features*4, num_features*16)
            self.conv3 = nn.SGConv(num_features*16, num_features*128)
        elif 'FeaStConv' in conv_method:
            self.conv1 = nn.FeaStConv(num_features, num_features*4, heads=5)
            self.conv2 = nn.FeaStConv(num_features*4, num_features*16, heads=5)
            self.conv3 = nn.FeaStConv(num_features*16, num_features*128, heads=5)
        elif 'Node2Vec' in conv_method:
            self.conv1 = nn.Node2Vec
        else:
            print("No valid model selected.")
            sys.stdout.flush()
            raise ValueError


        self.fc_g1 = torch.nn.Linear(num_features*1, 64)
        self.fc_g2 = torch.nn.Linear(64, GCN_num_outchannels)



        # DDI feature layers
        self.fc1 = torch.nn.Linear(num_drugs + GCN_num_outchannels, 64)
        self.fc2 = torch.nn.Linear(64, 64)
        self.fc3 = torch.nn.Linear(64, 64)
        self.fc4 = torch.nn.Linear(64, 64)
        self.fc5 = torch.nn.Linear(64, 1)

        # batch norm
        self.bn_g1 = torch.nn.BatchNorm1d(64)
        self.bn_g2 = torch.nn.BatchNorm1d(64)
        self.bn_1 = torch.nn.BatchNorm1d(64)
        self.bn_2 = torch.nn.BatchNorm1d(64)
        self.bn_3 = torch.nn.BatchNorm1d(64)
        self.bn_4 = torch.nn.BatchNorm1d(64)

        self.relu = torch.nn.ReLU()
        self.dropout = torch.nn.Dropout(dropout)

    def forward(self, PPI_data_object):
        # DDI_feature = PPI_data_object.DDI_features
        protein_mask = PPI_data_object.protein_mask
        PPI_x, PPI_edge_index, PPI_batch = PPI_data_object.x, PPI_data_object.edge_index, PPI_data_object.batch

        protein_mask = protein_mask.view((-1, 1, self.num_prots)).float()

        batch_size = protein_mask.size(0)


        # print(DDI_feature.shape)
        # print('protein_mask.size()', protein_mask.size())

        # PPI graph network

        PPI_x = self.conv1(PPI_x, PPI_edge_index)
        # PPI_x = F.relu(PPI_x)
        # PPI_x = self.conv2(PPI_x, PPI_edge_index)
        # PPI_x = F.relu(PPI_x)
        # PPI_x = F.relu(self.conv3(PPI_x, PPI_edge_index))

        PPI_x = PPI_x.view((batch_size, self.num_prots, PPI_x.shape[-1]))

        # PPI_x = torch_geometric.nn.global_max_pool(PPI_x, PPI_batch.view((batch_size, -1)))


        # multiply for flattening
        PPI_x = torch.bmm(protein_mask, PPI_x)
        PPI_x = PPI_x.view((batch_size, -1))

        return PPI_x

        # print('PPI_x.sum()', PPI_x.sum())

        # flatten

        x = self.relu(self.bn_g1(self.fc_g1(PPI_x)))
        # print('PPI_x.sum()', x.sum())
        # x_1 = self.relu(self.fc_g2(x) + PPI_x)

        # x = self.relu(self.fc1(x_1))
        # x_2 = self.relu(self.fc2(x) + x_1)
        # x = self.relu(self.fc3(x_2))
        # x_3 = self.relu(self.fc4(x) + x_2)

        x = self.relu(self.bn_1(self.fc1(x)))
        # print('x.sum()', x.sum())
        x = self.relu(self.bn_2(self.fc2(x)))
        # print('x.sum()', x.sum())
        x = self.relu(self.bn_3(self.fc3(x)))
        # print('x.sum()', x.sum())
        x = self.relu(self.bn_4(self.fc4(x)))
        x = self.fc5(x)

        # DDI feature network
        # DDI_x = torch.cat((PPI_x, DDI_feature), 1)
        # DDI_x = torch.sigmoid(DDI_x)
        return x

class ResTemplateNet(torch.nn.Module):
    def __init__(self, config, num_drugs, num_prots, num_features, conv_method, out_channels=128, dropout=0.2):
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

        gmp = torch_geometric.nn.global_max_pool
        gap = torch_geometric.nn.global_add_pool

        protein_mask = protein_mask.view((-1, 1, self.num_prots)).float()
        batch_size = protein_mask.size(0)

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

class CombinedProtFuncInteractionNetwork(torch.nn.Module):
    def __init__(self, config, epoch):
        super(CombinedProtFuncInteractionNetwork, self).__init__()

        self.config = config

        # initialize models
        if 'Res' in config.arch:
            self.interaction_model = ResTemplateNet(config,
                                                    num_drugs=0,
                                                    num_prots=config.num_proteins,
                                                    num_features=1,
                                                    conv_method=config.arch,
                                                    out_channels=128)
        else:
            self.interaction_model = TemplateSimpleNet(config,
                                                       num_drugs=0,
                                                       num_prots=config.num_proteins,
                                                       num_features=1,
                                                       conv_method=config.arch)

        self.protfunc_model = ProteinFunctionPredNet()

        print('Loading interaction model...')
        state_dict_path = '../models/PPI_network_' + (config.model_id +'_' if config.model_id else '') + config.arch + '_' + str(epoch) + '_epochs_model_fold_' + str(config.fold) + '.model'
        print(state_dict_path)
        self.interaction_model.load_state_dict(torch.load(state_dict_path))
        self.interaction_model.fc5 = torch.nn.Linear(128, 128)

        print('Loading protfunc model...')
        state_dict_path = '../models/protein_function_predictor/prot_func_pred_'+  'model_fold_'+str(config.fold)+'.model'
        print(state_dict_path)
        self.protfunc_model.load_state_dict(torch.load(state_dict_path))
        self.protfunc_model.fc5 = torch.nn.Linear(128, 128)
        print('Done.')

        self.fc1 = torch.nn.Linear(256, 1)

    def train(self):
        super(CombinedProtFuncInteractionNetwork, self).train()
        self.protfunc_model.train()
        self.interaction_model.train()

    def eval(self):
        super(CombinedProtFuncInteractionNetwork, self).eval()
        self.protfunc_model.eval()
        self.interaction_model.eval()

    def parameters(self):
        return list(self.interaction_model.parameters()) + list(self.protfunc_model.parameters()) + list(super(CombinedProtFuncInteractionNetwork, self).parameters())

    def forward(self, data):
        # protein_mask = data.protein_mask
        # x, edge_index, batch = data.x, data.edge_index, data.batch

        # protein_mask = protein_mask.view((-1, 1, self.num_prots)).float()
        # batch_size = protein_mask.size(0)

        inter_x = self.interaction_model(data)
        protfunc_x = self.protfunc_model(data.protfunc_data)

        x = torch.cat((inter_x, protfunc_x), 0)
        x = self.fc1(x)
        return x



class QuickTemplateSimpleNet(torch.nn.Module):
    def __init__(self, config, num_drugs, num_prots, num_features, conv_method, dropout=0.2):
        super(QuickTemplateSimpleNet, self).__init__()

        self.num_drugs = num_drugs
        self.num_prots = num_prots


        # mask feature

        # GCN laye4s
        if 'GCNConv' in conv_method:
            self.conv1 = nn.GCNConv(num_features, 4, cached=False)
            self.conv2 = nn.GCNConv(4, 16, cached=False)
            self.conv3 = nn.GCNConv(16, 32, cached=False)
            self.conv4 = nn.GCNConv(32, 64, cached=False)
            self.conv5 = nn.GCNConv(64, 128, cached=False)
            self.conv6 = nn.GCNConv(128, 128, cached=False)
            self.conv7 = nn.GCNConv(128, 128, cached=False)
            self.conv8 = nn.GCNConv(128, 1, cached=False)
        elif 'ChebConv' in conv_method:
            self.conv1 = nn.ChebConv(num_features, num_features*4, 3)
            self.conv2 = nn.ChebConv(num_features*4, num_features*16, 3)
            self.conv3 = nn.ChebConv(num_features*16, num_features*1, 3)
        elif 'SAGEConv' in conv_method:
            self.conv1 = nn.SAGEConv(num_features, num_features*4)
            self.conv2 = nn.SAGEConv(num_features*4, num_features*16)
            self.conv3 = nn.SAGEConv(num_features*16, num_features*1)
        elif 'GraphConv' in conv_method:
            self.conv1 = nn.GraphConv(num_features, num_features*4)
            self.conv2 = nn.GraphConv(num_features*4, num_features*16)
            self.conv3 = nn.GraphConv(num_features*16, num_features*1)
        elif 'GATConv' in conv_method:
            self.conv1 = nn.GATConv(num_features, num_features*4, heads=config.heads)
            self.conv2 = nn.GATConv(num_features*4*config.heads, num_features*16, heads=config.heads)
            self.conv3 = nn.GATConv(num_features*16*config.heads, num_features*1, heads=1)
        elif 'TAGConv' in conv_method:
            self.conv1 = nn.TAGConv(num_features, num_features*4)
            self.conv2 = nn.TAGConv(num_features*4, num_features*16)
            self.conv3 = nn.TAGConv(num_features*16, num_features*1)
        elif 'ARMAConv' in conv_method:
            self.conv1 = nn.ARMAConv(num_features, num_features*4)
            self.conv2 = nn.ARMAConv(num_features*4, num_features*16)
            self.conv3 = nn.ARMAConv(num_features*16, num_features*1)
        elif 'SGConv' in conv_method:
            self.conv1 = nn.SGConv(num_features, num_features*4)
            self.conv2 = nn.SGConv(num_features*4, num_features*16)
            self.conv3 = nn.SGConv(num_features*16, num_features*1)
        elif 'FeaStConv' in conv_method:
            self.conv1 = nn.FeaStConv(num_features, num_features*4, heads=5)
            self.conv2 = nn.FeaStConv(num_features*4, num_features*16, heads=5)
            self.conv3 = nn.FeaStConv(num_features*16, num_features*1, heads=5)
        elif 'SplineConv' in conv_method:
            self.conv1 = nn.SplineConv(num_features, 16, dim=1, kernel_size=5)
            self.conv2 = nn.SplineConv(16, 32, dim=1, kernel_size=5)
            self.conv3 = nn.SplineConv(32, 64, dim=1, kernel_size=7)
            self.conv4 = nn.SplineConv(64, 128, dim=1, kernel_size=7)
            self.conv5 = nn.SplineConv(128, 128, dim=1, kernel_size=11)
            self.conv6 = nn.SplineConv(128, 1, dim=1, kernel_size=11)
        elif 'Node2Vec' in conv_method:
            self.conv1 = nn.Node2Vec()
        else:
            print("No valid model selected.")
            sys.stdout.flush()
            raise ValueError

        self.bn_1 = nn.BatchNorm(4*config.heads)
        self.bn_2 = nn.BatchNorm(16*config.heads)

        self.relu = torch.nn.ReLU()
        self.dropout = torch.nn.Dropout(dropout)


    def forward(self, PPI_data_object):
        # DDI_feature = PPI_data_object.DDI_features
        PPI_x, PPI_edge_index, PPI_batch, edge_attr, edge_weight = PPI_data_object.x, PPI_data_object.edge_index, PPI_data_object.batch, PPI_data_object.edge_attr, PPI_data_object.edge_weight

        # print(DDI_feature.shape)
        # print('protein_mask.size()', protein_mask.size())

        # PPI graph network

        # batch_size = int(PPI_x.size(0)/self.num_prots)

        '''
        PPI_x = F.elu(self.conv1(PPI_x, PPI_edge_index, edge_weight=edge_weight))
        PPI_x = F.elu(self.conv2(PPI_x, PPI_edge_index, edge_weight=edge_weight))
        PPI_x = F.elu(self.conv3(PPI_x, PPI_edge_index, edge_weight=edge_weight))
        PPI_x = F.elu(self.conv4(PPI_x, PPI_edge_index, edge_weight=edge_weight))
        PPI_x = F.elu(self.conv5(PPI_x, PPI_edge_index, edge_weight=edge_weight))
        PPI_x = F.elu(self.conv6(PPI_x, PPI_edge_index, edge_weight=edge_weight))
        PPI_x = F.elu(self.conv7(PPI_x, PPI_edge_index, edge_weight=edge_weight))
        PPI_x = self.conv8(PPI_x, PPI_edge_index, edge_weight=edge_weight)

        PPI_x = PPI_x.view((-1, self.num_prots))

        return PPI_x
        '''

        edge_index = PPI_edge_index
        x = PPI_x
        x = F.elu(self.conv1(x, edge_index, edge_attr))
        x = F.elu(self.conv2(x, edge_index, edge_attr))
        x = F.elu(self.conv3(x, edge_index, edge_attr))
        x = self.conv4(x, edge_index, edge_attr)
        x = F.elu(self.conv5(x, edge_index, edge_attr))
        x = self.conv6(x, edge_index, edge_attr)

        # x = F.dropout(x, training=self.training)
        x = x.view((-1, self.num_prots))



        return x







