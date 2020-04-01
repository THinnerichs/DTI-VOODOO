import numpy as np
import networkx as nx
from tqdm import tqdm
import math

from sklearn.model_selection import KFold
from sklearn import metrics

from torch_dti_utils import *
from torch_networks import *

import torch
import torch.nn as nn
import torch_geometric.data as data

import argparse


def enlightened_missing_target_predictor(config,
                                         # results_filename='../results/torched_results_log',
                                         # num_epochs=3,
                                         # batch_size=512,
                                         plot=False,
                                         embedding_layer_sizes=[32, 64],
                                         embedding_method='gcn'):




    # activate device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    num_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0
    print('device', device)
    print('num_gpus', num_gpus)
    print(torch.cuda.is_available())

    raise Exception

    print("Loading data ...")
    dataset = FullNetworkDataset(num_proteins=config.num_proteins)
    # dataset is present in dimension (num_drugs * num_proteins)
    print("Finished.")

    # get full protein num
    config.num_proteins = None if config.num_proteins==-1 else config.num_proteins

    # generate indices for proteins
    kf = KFold(n_splits=config.num_folds, random_state=42, shuffle=True)
    X = np.zeros((dataset.num_proteins,1))

    # build for help matrix for indices
    help_matrix = np.arange(dataset.num_drugs * dataset.num_proteins)
    help_matrix = help_matrix.reshape((dataset.num_drugs, dataset.num_proteins))

    val_losses, accs, durations = [], [], []
    fold = 0
    for train_protein_indices, test_protein_indices in kf.split(X):
        fold += 1
        print("Fold:", fold)


        # build train data over whole dataset with help matrix
        train_indices = help_matrix[:, train_protein_indices].flatten()
        test_indices = help_matrix[:, test_protein_indices].flatten()
        print(train_indices.shape, test_indices.shape)

        train_dataset = dataset.get(train_indices)
        test_dataset = dataset.get(test_indices)

        train_size = int(0.8 * len(train_dataset))
        valid_size = len(train_dataset) - train_size
        train_dataset, valid_dataset = torch.utils.data.random_split(train_dataset, [train_size, valid_size])

        # build DataLoaders
        train_loader = data.DataLoader(train_dataset, config.batch_size, shuffle=True)
        valid_loader = data.DataLoader(valid_dataset, config.batch_size, shuffle=False)
        test_loader = data.DataLoader(test_dataset, config.batch_size, shuffle=False)

        # if torch.cuda.is_available():
            # torch.cuda.synchronize(device=device)

        model = SimpleConvGCN(num_drugs=dataset.num_drugs,
                              num_prots=dataset.num_proteins,
                              num_features=dataset.num_PPI_features)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
        model = nn.DataParallel(model, device_ids=list(range(num_gpus)))

        # storing best results
        best_loss = math.inf
        best_test_loss = math.inf
        best_epoch = -1
        best_test_ci = 0

        model_st = 'transductive_simple_node_feature'
        model_file_name = '../models/'+model_st+'_'+str(config.num_proteins)+'_model.model'
        results_file_name = '../results/'+model_st+'_'+str(config.num_proteins)+'_model.model'

        for epoch in range(1, config.num_epochs + 1):
            train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch)

            # print('{:02d}/{:03d}: Val Loss: {:.4f}, Test Accuracy: {:.3f}'.format(fold, epoch, val_loss, test_acc))
            print('Predicting for validation data...')
            labels, predictions = predicting(model, device, valid_loader)
            predictions = np.around(predictions)
            print('Validation:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                  metrics.accuracy_score(labels, predictions),
                  metrics.roc_auc_score(labels, predictions),
                  metrics.f1_score(labels, predictions),
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
                metrics_func_list = [metrics.accuracy_score, metrics.roc_auc_score, metrics.f1_score, metrics.matthews_corrcoef]
                metrics_list = [list_fun(G, P) for list_fun in metrics_func_list]
                ret += metrics_list

                # write results to results file
                with open(results_file_name, 'w') as f:
                    f.write(','.join(map(str, [fold]+ret)))
                best_test_loss = ret[1]
                best_test_ci = ret[-1]
                print('Test:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics_list)
                print('rmse improved at epoch ', best_epoch, '; best_test_loss, best_test_ci:', best_test_loss,
                      best_test_ci, model_st)
            else:
                print(ret[1], 'No improvement since epoch ', best_epoch, '; best_test_mse,best_test_ci:', best_test_loss,
                      best_test_ci, model_st)

        print('Done.')

    loss, acc = np.array(val_losses), np.array(accs)
    print("Acc mean:", acc.mean())

    print("Done.")


if __name__ == '__main__':
    # Add parser arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_proteins", type=int, default=-1)

    parser.add_argument("--num_epochs", type=int, default=3)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--num_folds", type=int, default=5)

    config = parser.parse_args()


    # Run classifier
    enlightened_missing_target_predictor(config)

