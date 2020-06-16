import numpy as np
import math
from sklearn.model_selection import KFold
from sklearn import metrics

from tqdm import tqdm
import argparse

import torch
import torch.nn as nn
import torch.utils.data as data

import DTI_data_preparation
from molecular_utils import train, predicting
import dti_utils


class HPODTIDataBuilder:
    def __init__(self, config):
        self.config = config

        # write data first
        print("Preparing protein data ...")
        # DTI_data_preparation.write_human_protein_list(min_score=config.PPI_min_score, mode=config.mode)
        # DTI_data_preparation.write_human_prot_func_protein_list(mode=config.mode)

        print("Loading data ...")
        self.drug_list = np.array(DTI_data_preparation.get_drug_list(config.mode))
        print(len(self.drug_list), "drugs present")
        self.protein_list = np.array(DTI_data_preparation.get_human_prot_func_proteins())[:config.num_proteins]
        print(len(self.protein_list), "proteins present\n")

        # PPI data
        print("Loading PPI graph ...")
        self.PPI_graph = DTI_data_preparation.get_PPI_DTI_graph_intersection()
        self.PPI_graph = self.PPI_graph.subgraph(self.protein_list)

        # calculate dimensions of network
        self.num_proteins = len(self.PPI_graph.nodes())
        self.num_drugs = len(self.drug_list)


    def build_data(self, config):


        # print('Building edge feature attributes ...')
        # forward_edge_feature_list = [1-self.PPI_graph[node1][node2]['score']/1000 for node1, node2 in list(self.PPI_graph.edges())]
        # backward_edge_feature_list = [1-self.PPI_graph[node1][node2]['score']/1000 for node2, node1 in list(self.PPI_graph.edges())]
        # self.edge_attr = torch.tensor(forward_edge_feature_list + backward_edge_feature_list, dtype=torch.float)# .view(-1,1)
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

        self.drug_features = DTI_data_preparation.get_DL2vec_features(self.drug_list)
        self.protein_features = DTI_data_preparation.get_DL2vec_features(self.protein_list)
        self.num_PPI_features = self.drug_features.shape[1]*2

        print('feature shape', self.drug_features.shape, self.protein_features.shape)

        print("Finished.\n")

    def get(self, indices):
        data_list = []

        # indices = list(range(self.num_drugs))
        for index in tqdm(indices):

            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            # build protein mask

            y = int(self.y_dti_data[drug_index, protein_index])

            # feature_array = torch.tensor(self.feature_matrix[drug_index, :], dtype=torch.float).view(-1, 1)
            # feature_array = torch.tensor(self.y_dti_data[drug_index, :], dtype=torch.float).view(-1,1)
            drug_feature = torch.tensor(self.drug_features[drug_index, :])
            protein_feature = torch.tensor(self.protein_features[protein_index, :])

            data_list.append((torch.cat((drug_feature, protein_feature), 0), y))

        return data_list

    def __len__(self):
        return self.num_drugs

class DTIGraphDataset(data.Dataset):
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


class HPOPredNet(nn.Module):
    def __init__(self):
        super(HPOPredNet, self).__init__()

        self.fc1 = nn.Linear(200, 128)
        self.fc2 = nn.Linear(128, 128)
        self.fc3 = nn.Linear(128, 128)
        self.fc4 = nn.Linear(128, 128)
        self.fc5 = nn.Linear(128, 1)

        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

        self.dropout = nn.Dropout(0.3)

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.relu(self.fc2(x))
        x = self.dropout(x)
        x = self.relu(self.fc3(x))
        x = self.dropout(x)
        x = self.relu(self.fc4(x))
        x = self.dropout(x)
        x = self.fc5(x)
        # x = self.sigmoid(x)

        return x


def siamese_drug_protein_network(config):
    model_st = 'siamese_drug_protein_netowrk'

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    if config.num_proteins <= 0:
        config.num_proteins = None

    dti_data = HPODTIDataBuilder(config)

    # generate indices for proteins
    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    X = np.zeros((dti_data.num_proteins, 1))

    # build for help matrix for indices
    help_matrix = np.arange(dti_data.num_drugs * dti_data.num_proteins)
    help_matrix = help_matrix.reshape((dti_data.num_drugs, dti_data.num_proteins))

    results = []
    fold = 0
    for train_protein_indices, test_protein_indices in kf.split(X):
        fold += 1
        print("Fold:", fold)
        if config.fold != -1 and config.fold != fold:
            continue

        dti_data.build_data(config)

        # build train data over whole dataset with help matrix
        train_indices = help_matrix[:, train_protein_indices].flatten()
        test_indices = help_matrix[:, test_protein_indices].flatten()
        print(train_indices.shape, test_indices.shape)

        train_dataset = dti_data.get(train_indices)
        test_dataset = dti_data.get(test_indices)

        train_dataset = DTIGraphDataset(train_dataset)
        test_dataset = DTIGraphDataset(test_dataset)

        # Calculate weights
        positives = dti_data.y_dti_data.flatten()[train_indices].sum()
        len_to_sum_ratio = (len(train_indices) - positives) / positives
        weight_dict = {0: 1.,
                       1: len_to_sum_ratio}

        train_loader = data.DataLoader(train_dataset, batch_size=config.batch_size, shuffle=True)
        test_loader = data.DataLoader(test_dataset, batch_size=config.batch_size)

        model = HPOPredNet()
        model = nn.DataParallel(model).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=config.lr)

        # storing best results
        best_loss = math.inf
        best_test_loss = math.inf
        best_epoch = -1
        best_test_ci = 0

        # model_file_name = '../models/' + model_st + '_' + config.node_features + '_' + str(
            # config.num_proteins) + '_fold_' + str(fold) + '_model.model'

        sys.stdout.flush()

        ret = None
        for epoch in range(1, config.num_epochs + 1):
            loss = train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch, weight_dict=weight_dict)
            print('Train Loss:', loss)

            if epoch%config.num_epochs == 10:
                print('Predicting for validation data...')
                # file='../results/protfunc_pred_results_' + str(config.num_epochs)+'_epochs'
                with open(file=file, mode='a') as f:
                    train_labels, train_predictions = predicting(model, device, train_loader)
                    print('Train: Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(train_labels, train_predictions),
                          dti_utils.dti_auroc(train_labels, train_predictions),
                          dti_utils.dti_f1_score(train_labels, train_predictions),
                          metrics.matthews_corrcoef(train_labels, train_predictions))#@TODO, file=f)

                    test_labels, test_predictions = predicting(model, device, test_loader)
                    print('Test: Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(test_labels, test_predictions),
                          dti_utils.dti_auroc(test_labels, test_predictions),
                          dti_utils.dti_f1_score(test_labels, test_predictions),
                          metrics.matthews_corrcoef(test_labels, test_predictions))#@TODO, file=f)

                    metrics_func_list = [metrics.accuracy_score, dti_utils.dti_auroc, dti_utils.dti_f1_score,
                                         metrics.matthews_corrcoef]
                    ret = [list_fun(test_labels, test_predictions) for list_fun in metrics_func_list]

                    test_loss = metrics.log_loss(test_labels, test_predictions, eps=0.000001)
                    if test_loss < best_loss:
                        best_loss = test_loss
                        best_epoch = epoch

                        print('rmse improved at epoch ', best_epoch, '; best_test_loss, best_test_ci:', best_test_loss,
                              best_test_ci, model_st)
                    else:
                        print(test_loss, 'No improvement since epoch ', best_epoch, ';', model_st)
            sys.stdout.flush()

        return

        print('Build overall data...')
        sys.stdout.flush()
        overall_dataset = dti_data.get(np.arange(dti_data.num_drugs * dti_data.num_proteins))
        overall_dataloader = data.DataLoader(overall_dataset, batch_size=config.batch_size)

        print('Predicting...')
        sys.stdout.flush()
        labels, predictions = predicting(model, device, overall_dataloader)
        filename = '../models/protein_function_predictor/pred_fold_'+str(fold)
        with open(file=filename+'.pkl', mode='wb') as f:
            pickle.dump(predictions, f, pickle.HIGHEST_PROTOCOL)


        model_filename = '../models/protein_function_predictor/prot_func_pred_'+ (config.model_id +'_' if config.model_id else '') + 'model_fold_'+str(fold)+'.model'
        torch.save(model.state_dict(), model_filename)
        print("Done.")
        sys.stdout.flush()

        results.append(ret)

    return

    results = np.array(results)
    results = [(results[:, i].mean(), results[:, i].std()) for i in range(results.shape[1])]

    results_file_name = '../results/protein_function_model' + '_'+ (config.model_id +'_' if config.model_id else '') + str(config.num_proteins) + '_model_results'

    print('Overall Results:')
    print('Model\tacc\tauroc\tf1\tmatt')
    print(model_st + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)))

    with open(results_file_name, 'a') as f:
        print('Model\tacc\tauroc\tf1\tmatt', file=f, end='\n')
        print(model_st + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)), file=f, end='\n')

    print("Done.")
    sys.stdout.flush()

if __name__ == '__main__':

    # Add parser arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_proteins", type=int, default=-1)
    # parser.add_argument("--arch", type=str, default='GCNConv')
    # parser.add_argument("--node_features", type=str, default='simple')

    parser.add_argument("--num_epochs", type=int, default=20)
    parser.add_argument("--batch_size", type=int, default=1024)
    # parser.add_argument("--num_folds", type=int, default=5)
    parser.add_argument("--lr", type=float, default=0.001)

    parser.add_argument("--model_id", type=str, default='')
    parser.add_argument("--model", type=str, default='protein')
    parser.add_argument("--mode", type=str, default='')
    parser.add_argument("--fold", type=int, default=1)

    config = parser.parse_args()

    if config.model=='protein':
        siamese_drug_protein_network(config)
    elif config.model=='drug':
        # drug_split_protein_function_predictor(config)
        pass
    else:
        print('No valid model selected...')
        raise ValueError
