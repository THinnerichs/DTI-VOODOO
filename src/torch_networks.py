import sys

import torch
import torch.nn.functional as F
import torch_geometric
import torch_geometric.nn as nn
from torch_geometric.transforms import add_self_loops

from protein_function_utils import ProteinFunctionPredNet
from HPO_predictor import HPOPredNet

class QuickTemplateNodeFeatureNet(torch.nn.Module):
    def __init__(self, config, num_drugs, num_prots, num_features, conv_method, dropout=0.2):
        super(QuickTemplateNodeFeatureNet, self).__init__()

        self.config = config

        self.num_drugs = num_drugs
        self.num_prots = num_prots
        self.num_features = num_features

        # mask feature

        # GCN layers
        if 'GCNConv' in conv_method:
            self.conv1 = nn.GCNConv(200, 200, cached=True, improved=True)
            self.conv2 = nn.GCNConv(200, 200, cached=True, improved=True)

            weight1 = torch.zeros((200,200))
            weight2 = torch.zeros((200,200))
            for i in range(200):
                weight1[i,i] = 1
                weight2[i,i] = 1

            self.conv1.weight = torch.nn.Parameter(weight1)
            self.conv2.weight = torch.nn.Parameter(weight2)

            # self.conv3 = nn.GCNConv(1, 1, cached=False, add_self_loops=True)
        elif 'GCNSideConv' in conv_method:
            self.conv1 = nn.GCNConv(200, 200, cached=True, add_self_loops=False)
            self.conv2 = nn.GCNConv(200, 200, cached=True, add_self_loops=False)
            self.conv3 = nn.GCNConv(200, 200, cached=True, add_self_loops=False)
        elif 'GENConv' in conv_method:
            conv1 = nn.GENConv(200,200, aggr='softmax', t=1.0, learn_t=True, num_layers=2, norm='layer')
            norm1 = torch.nn.LayerNorm(200, elementwise_affine=True)

            conv2 = nn.GENConv(200, 200, aggr='softmax', t=1.0, learn_t=True, num_layers=2, norm='layer')
            norm2 = torch.nn.LayerNorm(200, elementwise_affine=True)

            conv3 = nn.GENConv(200, 200, aggr='softmax', t=1.0, learn_t=True, num_layers=2, norm='layer')
            norm3 = torch.nn.LayerNorm(200, elementwise_affine=True)
            conv4 = nn.GENConv(200, 200, aggr='softmax', t=1.0, learn_t=True, num_layers=2, norm='layer')
            norm4 = torch.nn.LayerNorm(200, elementwise_affine=True)
            conv5 = nn.GENConv(200, 200, aggr='softmax', t=1.0, learn_t=True, num_layers=2, norm='layer')
            norm5 = torch.nn.LayerNorm(200, elementwise_affine=True)

            act = torch.nn.LeakyReLU(0.2, inplace=True)

            self.conv1 = nn.DeepGCNLayer(conv1, norm1, act, block='res', dropout=0.5)
            self.conv2 = nn.DeepGCNLayer(conv2, norm2, act, block='res', dropout=0.5)
            self.conv3 = nn.DeepGCNLayer(conv3, norm3, act, block='res', dropout=0.5)
            self.conv4 = nn.DeepGCNLayer(conv4, norm4, act, block='res', dropout=0.1)
            self.conv5 = nn.DeepGCNLayer(conv5, norm5, act, block='res', dropout=0.1)

        elif 'GATConv' in conv_method:
            self.conv1 = nn.GATConv(200, 200, heads=4, dropout=0.2, add_self_loops=False)
            self.conv2 = nn.GATConv(200*4, 200, heads=1, dropout=0.2, add_self_loops=False)
            # self.conv3 = nn.GATConv(8*2, 1, heads=1)
        elif 'APPNP' in conv_method:
            self.conv1 = nn.APPNP(K=50, alpha=0.15)
        else:
            print("No valid model selected.")
            sys.stdout.flush()
            raise ValueError

        state_dict_path = '../models/HPO_models/hpo_pred_fold_' + str(config.fold) + '_model'
        self.HPO_model = HPOPredNet(include_indications=self.config.include_indications)
        state_dict = torch.load(state_dict_path)
        from collections import OrderedDict
        new_state_dict = OrderedDict()
        for k, v in state_dict.items():
            name = k[7:]
            new_state_dict[name] = v
        self.HPO_model.load_state_dict(new_state_dict)

        for param in self.HPO_model.parameters():
            param.requires_grad = False

        self.mol_protein_model = torch.nn.Sequential(
            torch.nn.Linear(8192, 256),
            torch.nn.Dropout(0.5),
            # nn.BatchNorm1d(256),
            torch.nn.LeakyReLU(0.2, inplace=True),
            torch.nn.Linear(256, 200),
            # nn.Dropout(0.5),
            # nn.BatchNorm1d(50),
            # nn.LeakyReLU(0.2, inplace=True),
            # nn.Linear(256, 1),
            # nn.Sigmoid()
        )
        self.mol_drug_model = torch.nn.Sequential(
            torch.nn.Linear(1024, 256),
            torch.nn.Dropout(0.5),
            # nn.BatchNorm1d(256),
            torch.nn.LeakyReLU(0.2, inplace=True),
            torch.nn.Linear(256, 200),
            # nn.BatchNorm1d(50),
            # nn.Dropout(0.5),
            # nn.LeakyReLU(0.2, inplace=True),
            # nn.Linear(256, 1),
            # nn.Sigmoid()
        )
        self.protein_linear1 = torch.nn.Linear(400, 200)

        self.drug_linear1 = torch.nn.Linear(400, 200)
        self.drug_linear2 = torch.nn.Linear(200, 200)

        self.overall_linear1 = torch.nn.Linear(400, 200)
        self.overall_linear2 = torch.nn.Linear(200, 1)
        # self.overall_linear3 = torch.nn.Linear(16, 1)

        self.relu = torch.nn.ReLU()
        self.activation = torch.nn.LeakyReLU(0.2)
        self.sigmoid = torch.nn.Sigmoid()
        self.dropout = torch.nn.Dropout(dropout)

        self.sim = torch.nn.CosineSimilarity(dim=1)
        self.eps = 1e-7

    def forward(self, PPI_data_object):
        # DDI_feature = PPI_data_object.DDI_features
        PPI_x, PPI_edge_index, PPI_batch, edge_attr = PPI_data_object.x, PPI_data_object.edge_index, PPI_data_object.batch, PPI_data_object.edge_attr
        if self.config.include_indications:
            drug_feature = PPI_data_object.drug_feature.view(-1, self.num_features*2)
        else:
            drug_feature = PPI_data_object.drug_feature.view(-1, self.num_features)

        batch_size = drug_feature.size(0)

        drug_mol_feature = PPI_data_object.drug_mol_feature.view(batch_size, -1)

        PPI_x = self.HPO_model.model2(PPI_x)
        PPI_mol_x = self.mol_protein_model(PPI_data_object.protein_mol_feature)
        PPI_x = self.activation(torch.cat([PPI_x, PPI_mol_x], dim=1))
        PPI_x = self.protein_linear1(PPI_x)
        PPI_x = PPI_x.view(-1,200)

        # PPI_x = self.dropout(PPI_x)
        # PPI_x = F.elu(self.linear3(PPI_x))

        drug_feature = self.HPO_model.model(drug_feature).view(batch_size, 1, -1)
        # drug_mol_feature = drug_feature.repeat(1,self.num_prots,1).view(batch_size*self.num_prots,-1)
        drug_mol_feature = self.mol_drug_model(drug_mol_feature).view(batch_size, 1, -1)
        drug_feature = self.activation(torch.cat([drug_mol_feature, drug_feature], dim=2))
        drug_feature = self.activation(self.drug_linear1(drug_feature))
        drug_feature = self.drug_linear2(drug_feature)
        drug_feature = drug_feature.repeat(1,self.num_prots,1).view(batch_size*self.num_prots,-1)


        PPI_x = self.conv1(PPI_x, PPI_edge_index)
        PPI_x = self.conv2(PPI_x, PPI_edge_index)
        # PPI_x = self.conv3(PPI_x, PPI_edge_index)
        # PPI_x = self.conv4(PPI_x, PPI_edge_index)
        # PPI_x = self.conv5(PPI_x, PPI_edge_index)

        PPI_x = self.sim(drug_feature, PPI_x).unsqueeze(-1)
        cat_feature = PPI_x.view((-1, self.num_prots))

        return self.sigmoid(cat_feature)
