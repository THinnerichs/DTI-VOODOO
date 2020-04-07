import sys
from tqdm import tqdm


import numpy as np
import networkx as nx
import math


from sklearn.model_selection import KFold
from sklearn import metrics

from torch_dti_utils import *
from torch_networks import *

import torch
import torch_geometric.nn as nn
import torch_geometric.data as data

import argparse

import dti_utils



def transductive_missing_target_predictor(config,
                                          plot=False,
                                          embedding_layer_sizes=[32, 64]):

    # activate device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    num_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0

    # get full protein num
    config.num_proteins = None if config.num_proteins==-1 else config.num_proteins

    print("Loading data ...")
    network_data = None
    if config.node_features == 'simple':
        network_data = SimpleDTINetworkData(num_proteins=config.num_proteins)
    else:
        print("No valid node feature option selected.")
        sys.stdout.flush()
        raise ValueError
    # dataset is present in dimension (num_drugs * num_proteins)

    print("Finished.")

    # generate indices for proteins
    kf = KFold(n_splits=config.num_folds, random_state=42, shuffle=True)
    X = np.zeros((network_data.num_proteins,1))

    # build for help matrix for indices
    help_matrix = np.arange(network_data.num_drugs * network_data.num_proteins)
    help_matrix = help_matrix.reshape((network_data.num_drugs, network_data.num_proteins))

    print('Model:', config.arch)

    results = []
    fold = 0
    for train_protein_indices, test_protein_indices in kf.split(X):
        fold += 1
        print("Fold:", fold)

        # build train data over whole dataset with help matrix
        train_indices = help_matrix[:, train_protein_indices].flatten()
        test_indices = help_matrix[:, test_protein_indices].flatten()
        print(train_indices.shape, test_indices.shape)

        train_dataset = network_data.get(train_indices)
        test_dataset = network_data.get(test_indices)

        train_dataset = DTIGraphDataset(train_dataset)
        test_dataset = DTIGraphDataset(test_dataset)

        # Calculate weights
        len_to_sum_ratio = len(train_indices)/network_data.y_dti_data.flatten()[train_indices].sum()
        weight_dict = {0: 1.,
                       1: len_to_sum_ratio}


        # train_size = int(0.8 * len(train_dataset))
        # valid_size = len(train_dataset) - train_size
        # train_dataset, valid_dataset = torch.utils.data.random_split(train_dataset, [train_size, valid_size])


        # build DataLoaders
        train_loader = data.DataListLoader(train_dataset, config.batch_size, shuffle=True)
        # valid_loader = data.DataLoader(valid_dataset, config.batch_size, shuffle=False)
        valid_loader = None
        test_loader = data.DataListLoader(test_dataset, config.batch_size, shuffle=False)

        model = None
        if config.arch=='SimpleGCN':
            model = SimpleConvGCN(num_drugs=network_data.num_drugs,
                                  num_prots=network_data.num_proteins,
                                  num_features=network_data.num_PPI_features)
        elif config.arch=='TopKSimpleGCN':
            model = TopKPoolingSimpleGCN(num_drugs=network_data.num_drugs,
                                         num_prots=network_data.num_proteins,
                                         num_features=network_data.num_PPI_features)
        else:
            print("No valid architecture selected.")
            sys.stdout.flush()
            raise ValueError
        model = nn.DataParallel(model).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

        # storing best results
        best_loss = math.inf
        best_test_loss = math.inf
        best_epoch = -1
        best_test_ci = 0

        model_st = 'transductive_simple_node_feature'
        model_file_name = '../models/'+model_st+'_'+config.node_features + '_'+ str(config.num_proteins)+'_fold_'+str(fold)+'_model.model'

        sys.stdout.flush()

        ret = None
        for epoch in range(1, config.num_epochs + 1):
            train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch, weight_dict=weight_dict)

            print('Predicting for validation data...')
            labels, predictions = predicting(model, device, test_loader)
            predictions = np.around(predictions)
            print('labels', labels, 'predictions', predictions)
            print(labels.min(), labels.max(), predictions.min(), predictions.max())
            print('Validation:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                  metrics.accuracy_score(labels, predictions),
                  dti_utils.dti_auroc(labels, predictions),
                  dti_utils.dti_f1_score(labels, predictions),
                  metrics.matthews_corrcoef(labels, predictions))

            val = mse(labels, predictions)
            if val < best_loss:
                best_loss = val
                best_epoch = epoch + 1
                torch.save(model.state_dict(), model_file_name)
                print('predicting for test data')
                G, P = predicting(model, device, test_loader)
                ret = [rmse(G, P), mse(G, P), pearson(G, P), spearman(G, P), ci(G, P)]
                P = np.around(P)
                metrics_func_list = [metrics.accuracy_score, dti_utils.dti_auroc, dti_utils.dti_f1_score, metrics.matthews_corrcoef]
                metrics_list = [list_fun(G, P) for list_fun in metrics_func_list]
                ret += metrics_list

                # write results to results file
                best_test_loss = ret[1]
                best_test_ci = ret[-1]
                print('Test:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics_list)
                print('rmse improved at epoch ', best_epoch, '; best_test_loss, best_test_ci:', best_test_loss,
                      best_test_ci, model_st)
            else:
                print(ret[1], 'No improvement since epoch ', best_epoch, '; best_test_mse,best_test_ci:', best_test_loss,
                      best_test_ci, model_st)
            sys.stdout.flush()
        results.append(ret)

    results_file_name = '../results/' + config.arch + '_' + config.node_features + '_' + str(config.num_proteins) + '_model_results'

    results = np.array(results)
    results = [(results[:, i].mean(), results[:, i].std()) for i in range(results.shape[1])]

    print('Overall Results:')
    print('Model\trmse\tmse\tpearson\tspearman\tacc\tauroc\tf1\tmatt')
    print(config.arch+'\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)))

    with open(results_file_name, 'a') as f:
        print('Model\trmse\tmse\tpearson\tspearman\tacc\tauroc\tf1\tmatt', file=f)
        print(config.arch+'\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)), file=f)

    print("Done.")
    sys.stdout.flush()

def inductive_missing_target_predictor(config,
                                       ):
    pass



if __name__ == '__main__':
    # Add parser arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_proteins", type=int, default=-1)
    parser.add_argument("--arch", type=str, default='SimpleGCN')
    parser.add_argument("--node_features", type=str, default='simple')


    parser.add_argument("--num_epochs", type=int, default=3)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--num_folds", type=int, default=5)

    config = parser.parse_args()


    # Run classifier
    transductive_missing_target_predictor(config)

