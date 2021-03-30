import sys
from tqdm import tqdm
import time

import numpy as np
import networkx as nx
import math


from sklearn.model_selection import KFold
from sklearn import metrics

import torch
import torch_geometric.nn as nn
import torch_geometric.data as data

import argparse

# import other files
from torch_dti_utils import *
from torch_networks import *

from protein_function_utils import ProteinFunctionPredNet, ProteinFunctionDTIDataBuilder
import dti_utils

import DTI_data_preparation
import PPI_utils
import DDI_utils


def quickened_missing_target_predictor(config,
                                       plot=False):

    # activate device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    num_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0
    config.num_gpus = num_gpus

    np.random.seed(42)

    # build data
    config.num_proteins = None if config.num_proteins==-1 else config.num_proteins
    network_data = QuickProtFuncDTINetworkData(config=config)

    # get full protein num
    num_drugs = network_data.num_drugs
    num_proteins = network_data.num_proteins

    # dataset is present in dimension (num_drugs * num_proteins)

    # generate indices for proteins
    kf = KFold(n_splits=config.num_folds, random_state=42, shuffle=True)
    X = np.zeros((num_proteins,1))

    # build for help matrix for indices
    help_matrix = np.arange(num_drugs * num_proteins)
    help_matrix = help_matrix.reshape((num_drugs, num_proteins))

    print('Model:', config.arch)


    results = []
    fold = 0
    for train_protein_indices, test_protein_indices in kf.split(X):
        fold += 1
        if config.fold != -1 and fold != config.fold:
            continue
        print("Fold:", fold)

        config.train_prots = train_protein_indices
        network_data.build_data(config)

        # build train data over whole dataset with help matrix

        print('Fetching data...')
        train_dataset = network_data.get()
        print('\nDone.\n')

        train_dataset = DTIGraphDataset(train_dataset)

        # Calculate weights
        positives = network_data.y_dti_data[:, train_protein_indices].sum()
        len_to_sum_ratio = (network_data.num_drugs * len(train_protein_indices)-positives)/positives #Get negatives/positives ratio
        print('Neg/pos ratio:', len_to_sum_ratio)

        train_mask = network_data.train_mask
        test_mask = 1-network_data.train_mask

        # train_size = int(0.8 * len(train_dataset))
        # valid_size = len(train_dataset) - train_size
        # train_dataset, valid_dataset = torch.utils.data.random_split(train_dataset, [train_size, valid_size])

        # build DataLoaders
        print('Building data loader ...')
        train_loader = data.DataListLoader(train_dataset, config.batch_size, shuffle=True)
        # valid_loader = data.DataLoader(valid_dataset, config.batch_size, shuffle=False)

        print('Initializing model ...')
        model = QuickTemplateNodeFeatureNet(config,
                                            num_drugs=0,
                                            num_prots=network_data.num_proteins,
                                            num_features=network_data.num_PPI_features,
                                            conv_method=config.arch)
        model = nn.DataParallel(model).to(device)
        print("model total parameters", sum(p.numel() for p in model.parameters()))

        optimizer = torch.optim.Adam(model.parameters(), lr=config.lr)#, momentum=0.9)

        # storing best results
        best_loss = math.inf
        best_test_loss = math.inf
        best_epoch = -1
        best_test_ci = 0

        model_st = 'transductive_simple_node_feature'
        # model_file_name = '../models/'+model_st+'_'+config.node_features + '_'+ str(config.num_proteins)+'_fold_'+str(fold)+'_model.model'

        sys.stdout.flush()

        ret = None
        best_AUROC = 0
        for epoch in range(1, config.num_epochs + 1):
            loss = quick_train(config=config,
                               model=model,
                               device=device,
                               train_loader=train_loader,
                               optimizer=optimizer,
                               epoch=epoch,
                               neg_to_pos_ratio=len_to_sum_ratio,
                               train_mask=train_mask)
            print('Train loss:', loss)
            sys.stdout.flush()

            if epoch%3 == 0:
                print('Predicting for validation data...')
                file='../results/quick_pred_' +config.arch+'_fold_'+str(config.fold) + '_results'
                with open(file=file, mode='a') as f:
                    labels, predictions = quick_predicting(model, device, train_loader, round=False)

                    # get train and test predictions
                    train_labels = labels.reshape((num_drugs, num_proteins))[:, train_mask==1].flatten()
                    train_predictions = predictions.reshape((num_drugs, num_proteins))[:, train_mask==1].flatten()

                    test_labels = labels.reshape((num_drugs, num_proteins))[:, train_mask==0].flatten()
                    test_predictions = predictions.reshape((num_drugs, num_proteins))[:, train_mask==0].flatten()

                    print('pred_eval', train_labels.max(), train_predictions.max(), train_labels.min(), train_predictions.min(), train_predictions.shape)
                    print('pred_eval', test_labels.max(), test_predictions.max(), test_labels.min(), test_predictions.min(), test_predictions.shape)

                    print('Train:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(train_labels, train_predictions.round()),
                          dti_utils.dti_auroc(train_labels, train_predictions),
                          dti_utils.micro_AUC_per_prot(train_labels, train_predictions, config.num_drugs),
                          dti_utils.dti_f1_score(train_labels, train_predictions.round()),
                          metrics.matthews_corrcoef(train_labels, train_predictions.round()), file=f)

                    print('Test:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(test_labels, test_predictions.round()),
                          dti_utils.dti_auroc(test_labels, test_predictions),
                          dti_utils.micro_AUC_per_prot(test_labels, test_predictions, config.num_drugs),
                          dti_utils.dti_f1_score(test_labels, test_predictions.round()),
                          metrics.matthews_corrcoef(test_labels, test_predictions.round()), file=f)

                    test_AUROC = dti_utils.dti_auroc(test_labels, test_predictions)
                    if test_AUROC > best_AUROC:
                        model_filename = '../models/graph_models/PPI_network_model_with_mol_features_fold_' + str(fold) + '.model'
                        torch.save(model.state_dict(), model_filename)
            if epoch == 50:

                drug_list_repeated = network_data.drug_list.repeat(config.num_proteins)
                protein_list_repeated = network_data.protein_list.reshape(1,-1).repeat(config.num_drugs, axis=0).reshape(-1)

                pred_loader = data.DataListLoader(train_dataset, config.batch_size, shuffle=False)
                labels, predictions = quick_predicting(model, device, pred_loader, round=False)

                zipped_list = list(zip(drug_list_repeated, protein_list_repeated, labels, predictions))

                pred_filename = '../models/graph_models/PPI_network_model_with_mol_features_fold_' + str(fold) + '_predictions.pkl'
                with open(file=pred_filename, mode='wb') as f:
                    pickle.dump(zipped_list, f, pickle.HIGHEST_PROTOCOL)

            sys.stdout.flush()
        results.append(ret)

        if torch.cuda.is_available():
            torch.cuda.synchronize()

        # model_filename = '../models/PPI_network_' + config.arch + '_model_fold_' + str(fold) + '.model'
        # torch.save(model.state_dict(), model_filename)

    return

if __name__ == '__main__':
    # Add parser arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_proteins", type=int, default=-1)
    parser.add_argument("--arch", type=str, default='GCNConv')
    parser.add_argument("--node_features", type=str, default='MolPred')


    parser.add_argument("--num_epochs", type=int, default=3)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--num_folds", type=int, default=5)
    parser.add_argument("--lr", type=float, default=0.001)

    parser.add_argument("--fold", type=int, default=-1)

    parser.add_argument("--split_mode", type=str, default='standard')

    parser.add_argument("--mode", type=str, default='')
    parser.add_argument("--PPI_min_score", type=int, default=700)

    parser.add_argument("--drug_mode", type=str, default='trfm')
    parser.add_argument("--include_mol_features", action='store_true')

    parser.add_argument("--yamanishi_test", action='store_true')
    parser.add_argument("--biosnap_test", action='store_true')
    parser.add_argument("--include_indications", action='store_true')


    config = parser.parse_args()

    # Run classifier
    quickened_missing_target_predictor(config)


