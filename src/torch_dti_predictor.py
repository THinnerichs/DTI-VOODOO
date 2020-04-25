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

import DTI_data_preparation
import PPI_utils
import DDI_utils



def transductive_missing_target_predictor(config,
                                          plot=False):

    # activate device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    num_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0

    np.random.seed(42)

    # get full protein num
    config.num_proteins = None if config.num_proteins==-1 else config.num_proteins
    num_proteins = config.num_proteins if config.num_proteins else 11574
    num_drugs = 641

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

        network_data = MolPredDTINetworkData(config=config)

        # build train data over whole dataset with help matrix
        train_indices = help_matrix[:, train_protein_indices].flatten()
        test_indices = help_matrix[:, test_protein_indices].flatten()
        print(train_indices.shape, test_indices.shape)

        train_labels = network_data.get_labels(train_indices)
        zipped_label_ind_array = list(zip(train_indices, train_labels))
        positive_label_indices = np.array([index for index, label in zipped_label_ind_array if label == 1])
        negative_label_indices = np.array([index for index, label in zipped_label_ind_array if label == 0])
        train_indices = np.random.choice(negative_label_indices, int(config.neg_sample_ratio * len(negative_label_indices)))
        train_indices = np.concatenate((train_indices, positive_label_indices), axis=0)

        print(train_indices.shape, test_indices.shape)

        print('Fetching train data...\n')
        train_dataset = network_data.get(train_indices)
        print('Fetching test data...\n')
        test_dataset = network_data.get(test_indices)

        train_dataset = DTIGraphDataset(train_dataset)
        test_dataset = DTIGraphDataset(test_dataset)

        # Calculate weights
        positives = network_data.y_dti_data.flatten()[train_indices].sum()
        len_to_sum_ratio = (len(train_indices)-positives)/positives #Get negatives/positives ratio
        weight_dict = {0: 1.,
                       1: len_to_sum_ratio}

        # train_size = int(0.8 * len(train_dataset))
        # valid_size = len(train_dataset) - train_size
        # train_dataset, valid_dataset = torch.utils.data.random_split(train_dataset, [train_size, valid_size])

        # build DataLoaders
        train_loader = data.DataListLoader(train_dataset, config.batch_size, shuffle=True)
        # valid_loader = data.DataLoader(valid_dataset, config.batch_size, shuffle=False)
        test_loader = data.DataListLoader(test_dataset, config.batch_size, shuffle=False)

        model = None
        if 'Res' not in config.arch:
            model = TemplateSimpleNet(num_drugs=network_data.num_drugs,
                                      num_prots=network_data.num_proteins,
                                      num_features=network_data.num_PPI_features,
                                      conv_method=config.arch)
        elif 'Res' in config.arch:
            model = ResTemplateNet(num_drugs=network_data.num_drugs,
                                   num_prots=network_data.num_proteins,
                                   num_features=network_data.num_PPI_features,
                                   conv_method=config.arch,
                                   out_channels=64)

        model = nn.DataParallel(model).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

        # storing best results
        best_loss = math.inf
        best_test_loss = math.inf
        best_epoch = -1
        best_test_ci = 0

        model_st = 'transductive_simple_node_feature'
        # model_file_name = '../models/'+model_st+'_'+config.node_features + '_'+ str(config.num_proteins)+'_fold_'+str(fold)+'_model.model'

        sys.stdout.flush()

        ret = None
        for epoch in range(1, config.num_epochs + 1):
            loss = train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch, weight_dict=weight_dict)
            print('Train loss:', loss)
            sys.stdout.flush()

            if epoch%10 == 0:
                print('Predicting for validation data...')
                file='../results/full_interactions_results_' +config.arch+'_'+ str(num_proteins) + '_prots_'+str(epoch)+'_epochs'
                with open(file=file, mode='a') as f:
                    train_labels, train_predictions = predicting(model, device, train_loader)
                    print('Train:', config.neg_sample_ratio,'Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(train_labels, train_predictions),
                          dti_utils.dti_auroc(train_labels, train_predictions),
                          dti_utils.dti_f1_score(train_labels, train_predictions),
                          metrics.matthews_corrcoef(train_labels, train_predictions), file=f)

                    test_labels, test_predictions = predicting(model, device, test_loader)
                    print('Test:', config.neg_sample_ratio,'Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(test_labels, test_predictions),
                          dti_utils.dti_auroc(test_labels, test_predictions),
                          dti_utils.dti_f1_score(test_labels, test_predictions),
                          metrics.matthews_corrcoef(test_labels, test_predictions), file=f)

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
        results.append(ret)

        if torch.cuda.is_available():
            torch.cuda.synchronize()

    return
    results_file_name = '../results/' + config.arch + '_' + config.node_features + '_' + str(config.num_proteins) + '_model_results'

    results = np.array(results)
    results = [(results[:, i].mean(), results[:, i].std()) for i in range(results.shape[1])]

    print('Overall Results:')
    print('Model\tacc\tauroc\tf1\tmatt')
    print(config.arch+'\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)))

    with open(results_file_name, 'a') as f:
        print('Model\tacc\tauroc\tf1\tmatt', file=f)
        print(config.arch+'\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)), file=f)

    print("Done.")
    sys.stdout.flush()


def transductive_missing_drug_predictor(config,
                                        plot=False):
    # activate device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    num_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0

    # get full protein num
    config.num_proteins = None if config.num_proteins == -1 else config.num_proteins
    num_proteins = config.num_proteins if config.num_proteins else 11574
    num_drugs = 641

    # dataset is present in dimension (num_drugs * num_proteins)

    # generate indices for proteins
    kf = KFold(n_splits=config.num_folds, random_state=42, shuffle=True)
    X = np.zeros((num_drugs, 1))

    # build for help matrix for indices
    help_matrix = np.arange(num_drugs * num_proteins)
    help_matrix = help_matrix.reshape((num_drugs, num_proteins))

    print('Model:', config.arch)

    results = []
    fold = 0
    for train_drug_indices, test_drug_indices in kf.split(X):
        fold += 1
        if config.fold != -1 and fold != config.fold:
            continue
        print("Fold:", fold)

        network_data = MolPredDTINetworkData(config=config)

        # build train data over whole dataset with help matrix
        train_indices = help_matrix[train_drug_indices, :].flatten()
        test_indices = help_matrix[test_drug_indices, :].flatten()
        print(train_indices.shape, test_indices.shape)


        train_labels = network_data.get_labels(train_indices)
        zipped_label_ind_array = list(zip(train_indices, train_labels))
        positive_label_indices = np.array([index for index, label in zipped_label_ind_array if label == 1])
        negative_label_indices = np.array([index for index, label in zipped_label_ind_array if label == 0])
        train_indices = np.random.choice(negative_label_indices, int(config.neg_sample_ratio * len(negative_label_indices)))
        train_indices = np.concatenate((train_indices, positive_label_indices), axis=0)

        print(train_indices.shape, test_indices.shape)


        train_dataset = network_data.get(train_indices)
        test_dataset = network_data.get(test_indices)

        train_dataset = DTIGraphDataset(train_dataset)
        test_dataset = DTIGraphDataset(test_dataset)

        # Calculate weights
        positives = network_data.y_dti_data.flatten()[train_indices].sum()
        len_to_sum_ratio = (len(train_indices) - positives) / positives
        weight_dict = {0: 1.,
                       1: len_to_sum_ratio}

        # train_size = int(0.8 * len(train_dataset))
        # valid_size = len(train_dataset) - train_size
        # train_dataset, valid_dataset = torch.utils.data.random_split(train_dataset, [train_size, valid_size])

        # build DataLoaders
        train_loader = data.DataListLoader(train_dataset, config.batch_size, shuffle=True)
        # valid_loader = data.DataLoader(valid_dataset, config.batch_size, shuffle=False)
        test_loader = data.DataListLoader(test_dataset, config.batch_size, shuffle=False)

        model = None
        if 'Res' not in config.arch:
            model = TemplateSimpleNet(num_drugs=network_data.num_drugs,
                                      num_prots=network_data.num_proteins,
                                      num_features=network_data.num_PPI_features,
                                      conv_method=config.arch)
        elif 'Res' in config.arch:
            model = ResTemplateNet(num_drugs=network_data.num_drugs,
                                   num_prots=network_data.num_proteins,
                                   num_features=network_data.num_PPI_features,
                                   conv_method=config.arch)

        model = nn.DataParallel(model).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

        # storing best results
        best_loss = math.inf
        best_test_loss = math.inf
        best_epoch = -1
        best_test_ci = 0

        model_st = 'transductive_simple_node_feature'
        # model_file_name = '../models/'+model_st+'_'+config.node_features + '_'+ str(config.num_proteins)+'_fold_'+str(fold)+'_model.model'

        sys.stdout.flush()

        ret = None
        for epoch in range(1, config.num_epochs + 1):
            loss = train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch,
                         weight_dict=weight_dict)
            print('Train loss:', loss)

            if epoch % 10 == 0:
                print('Predicting for validation data...')
                file = '../results/full_interactions_drug_split_results_' + config.arch + '_' + str(num_proteins) + '_prots_' + str(epoch) + '_epochs'
                with open(file=file, mode='a') as f:
                    train_labels, train_predictions = predicting(model, device, train_loader)
                    print('Train:', config.neg_sample_ratio, 'Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(train_labels, train_predictions),
                          dti_utils.dti_auroc(train_labels, train_predictions),
                          dti_utils.dti_f1_score(train_labels, train_predictions),
                          metrics.matthews_corrcoef(train_labels, train_predictions), file=f)

                    test_labels, test_predictions = predicting(model, device, test_loader)
                    print('Test:', config.neg_sample_ratio, 'Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(test_labels, test_predictions),
                          dti_utils.dti_auroc(test_labels, test_predictions),
                          dti_utils.dti_f1_score(test_labels, test_predictions),
                          metrics.matthews_corrcoef(test_labels, test_predictions), file=f)

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
        results.append(ret)

        if torch.cuda.is_available():
            torch.cuda.synchronize()


def test_predictor_on_drughub_protein_data(config):
    # activate device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    num_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0

    drug_list = DTI_data_preparation.get_drug_list()
    protein_list = DTI_data_preparation.get_human_proteins()

    drughub_drugs = PPI_utils.get_drughub_drug_list()
    drughub_proteins = set(PPI_utils.get_drughub_protein_list()) & set(protein_list)

    test_protein_indices = np.array(sorted([protein_list.index(protein) for protein in drughub_proteins]))
    train_protein_indices = np.array([i for i in range(len(protein_list)) if i not in test_protein_indices])

    test_drug_indices = np.array(sorted([list(drug_list).index(drug) for drug in drughub_drugs if drug in drug_list]))
    train_drug_indices = np.array([i for i in range(len(drug_list)) if i not in test_drug_indices])

    drughub_dti_graph = PPI_utils.get_drughub_dti_graph()

    y_drughub_dti_data = np.array([int(drughub_dti_graph.has_edge(drug_list[i], protein_list[j]))
                                   for i in test_drug_indices for j in test_protein_indices])
    print('y_data.shape', y_drughub_dti_data.shape)


    # get full protein num
    config.num_proteins = None if config.num_proteins == -1 else config.num_proteins
    num_proteins = config.num_proteins if config.num_proteins else 11574
    num_drugs = 641

    # dataset is present in dimension (num_drugs * num_proteins)

    # generate indices for proteins

    # build for help matrix for indices
    help_matrix = np.arange(num_drugs * num_proteins)
    help_matrix = help_matrix.reshape((num_drugs, num_proteins))

    print('Model:', config.arch)

    results = []

    config.fold = 1
    network_data = MolPredDTINetworkData(config=config)

    # build train data over whole dataset with help matrix
    train_indices = help_matrix[:, train_protein_indices].flatten()
    test_indices = help_matrix[:, test_protein_indices][test_drug_indices, :].flatten()
    print(train_indices.shape, test_indices.shape)

    train_labels = network_data.get_labels(train_indices)
    zipped_label_ind_array = list(zip(train_indices, train_labels))
    positive_label_indices = np.array([index for index, label in zipped_label_ind_array if label == 1])
    negative_label_indices = np.array([index for index, label in zipped_label_ind_array if label == 0])
    train_indices = np.random.choice(negative_label_indices, int(config.neg_sample_ratio * len(negative_label_indices)))
    train_indices = np.concatenate((train_indices, positive_label_indices), axis=0)

    print(train_indices.shape, test_indices.shape)

    train_dataset = network_data.get(train_indices)
    test_dataset = network_data.get(test_indices)

    for i, data_obj  in enumerate(test_dataset):
        data_obj.y = y_drughub_dti_data[i]

    train_dataset = DTIGraphDataset(train_dataset)
    test_dataset = DTIGraphDataset(test_dataset)

    # Calculate weights
    positives = network_data.y_dti_data.flatten()[train_indices].sum()
    len_to_sum_ratio = (len(train_indices) - positives) / positives
    weight_dict = {0: 1.,
                   1: len_to_sum_ratio}

    # train_size = int(0.8 * len(train_dataset))
    # valid_size = len(train_dataset) - train_size
    # train_dataset, valid_dataset = torch.utils.data.random_split(train_dataset, [train_size, valid_size])

    # build DataLoaders
    train_loader = data.DataListLoader(train_dataset, config.batch_size, shuffle=True)
    # valid_loader = data.DataLoader(valid_dataset, config.batch_size, shuffle=False)
    test_loader = data.DataListLoader(test_dataset, config.batch_size, shuffle=False)

    model = None
    if 'Res' not in config.arch:
        model = TemplateSimpleNet(config,
                                  num_drugs=network_data.num_drugs,
                                  num_prots=network_data.num_proteins,
                                  num_features=network_data.num_PPI_features,
                                  conv_method=config.arch)
    elif 'Res' in config.arch:
        model = ResTemplateNet(num_drugs=network_data.num_drugs,
                               num_prots=network_data.num_proteins,
                               num_features=network_data.num_PPI_features,
                               conv_method=config.arch)
    model = nn.DataParallel(model).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

    # storing best results
    best_loss = math.inf
    best_test_loss = math.inf
    best_epoch = -1
    best_test_ci = 0

    # model_file_name = '../models/'+model_st+'_'+config.node_features + '_'+ str(config.num_proteins)+'_fold_'+str(fold)+'_model.model'

    sys.stdout.flush()

    ret = None
    for epoch in range(1, config.num_epochs + 1):
        loss = train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch,
                     weight_dict=weight_dict)
        print('Train loss:', loss)

        if epoch % 10 == 0:
            print('Predicting for validation data...')
            file = '../results/drughub_protein_split_' + config.arch + '_prots_' + str(epoch) + '_epochs'
            with open(file=file, mode='a') as f:
                train_labels, train_predictions = predicting(model, device, train_loader)
                print('Train:', config.neg_sample_ratio,'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(train_labels, train_predictions),
                      dti_utils.dti_auroc(train_labels, train_predictions),
                      dti_utils.dti_f1_score(train_labels, train_predictions),
                      metrics.matthews_corrcoef(train_labels, train_predictions), file=f)

                test_labels, test_predictions = predicting(model, device, test_loader)
                print('Test:', config.neg_sample_ratio,'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(test_labels, test_predictions),
                      dti_utils.dti_auroc(test_labels, test_predictions),
                      dti_utils.dti_f1_score(test_labels, test_predictions),
                      metrics.matthews_corrcoef(test_labels, test_predictions), file=f)

                metrics_func_list = [metrics.accuracy_score, dti_utils.dti_auroc, dti_utils.dti_f1_score,
                                     metrics.matthews_corrcoef]
                ret = [list_fun(test_labels, test_predictions) for list_fun in metrics_func_list]

                test_loss = metrics.log_loss(test_labels, test_predictions, eps=0.000001)

        sys.stdout.flush()
    results.append(ret)

    if torch.cuda.is_available():
        torch.cuda.synchronize()


def test_predictor_on_drughub_drug_data(config):
    # activate device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    num_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0

    drug_list = DTI_data_preparation.get_drug_list()
    protein_list = DTI_data_preparation.get_human_proteins()

    drughub_drugs = PPI_utils.get_drughub_drug_list()
    drughub_proteins = set(PPI_utils.get_drughub_protein_list()) & set(protein_list)

    test_protein_indices = np.array(sorted([protein_list.index(protein) for protein in drughub_proteins]))
    train_protein_indices = np.array([i for i in range(len(protein_list)) if i not in test_protein_indices])

    test_drug_indices = np.array(sorted([list(drug_list).index(drug) for drug in drughub_drugs if drug in drug_list]))
    train_drug_indices = np.array([i for i in range(len(drug_list)) if i not in test_drug_indices])

    drughub_dti_graph = PPI_utils.get_drughub_dti_graph()

    y_drughub_dti_data = np.array([int(drughub_dti_graph.has_edge(drug_list[i], protein_list[j]))
                                   for i in test_drug_indices for j in test_protein_indices])
    print('y_data.shape', y_drughub_dti_data.shape)


    # get full protein num
    config.num_proteins = None if config.num_proteins == -1 else config.num_proteins
    num_proteins = config.num_proteins if config.num_proteins else 11574
    num_drugs = 641

    # dataset is present in dimension (num_drugs * num_proteins)

    # generate indices for proteins

    # build for help matrix for indices
    help_matrix = np.arange(num_drugs * num_proteins)
    help_matrix = help_matrix.reshape((num_drugs, num_proteins))

    print('Model:', config.arch)

    results = []

    config.fold = 1
    network_data = MolPredDTINetworkData(config=config)

    # build train data over whole dataset with help matrix
    train_indices = help_matrix[train_drug_indices, :].flatten()
    test_indices = help_matrix[test_drug_indices, :][:, test_protein_indices].flatten()

    print(train_indices.shape, test_indices.shape)

    train_dataset = network_data.get(train_indices)
    test_dataset = network_data.get(test_indices)

    for i, data_obj in enumerate(test_dataset):
        data_obj.y = y_drughub_dti_data[i]

    train_dataset = DTIGraphDataset(train_dataset)
    test_dataset = DTIGraphDataset(test_dataset)

    # Calculate weights
    positives = network_data.y_dti_data.flatten()[train_indices].sum()
    len_to_sum_ratio = (len(train_indices) - positives) / positives
    weight_dict = {0: 1.,
                   1: len_to_sum_ratio}

    # train_size = int(0.8 * len(train_dataset))
    # valid_size = len(train_dataset) - train_size
    # train_dataset, valid_dataset = torch.utils.data.random_split(train_dataset, [train_size, valid_size])

    # build DataLoaders
    train_loader = data.DataListLoader(train_dataset, config.batch_size, shuffle=True)
    # valid_loader = data.DataLoader(valid_dataset, config.batch_size, shuffle=False)
    test_loader = data.DataListLoader(test_dataset, config.batch_size, shuffle=False)

    model = None
    if 'Res' not in config.arch:
        model = TemplateSimpleNet(num_drugs=network_data.num_drugs,
                                  num_prots=network_data.num_proteins,
                                  num_features=network_data.num_PPI_features,
                                  conv_method=config.arch)
    elif 'Res' in config.arch:
        model = ResTemplateNet(num_drugs=network_data.num_drugs,
                               num_prots=network_data.num_proteins,
                               num_features=network_data.num_PPI_features,
                               conv_method=config.arch)
    model = nn.DataParallel(model).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

    # storing best results
    best_loss = math.inf
    best_test_loss = math.inf
    best_epoch = -1
    best_test_ci = 0

    # model_file_name = '../models/'+model_st+'_'+config.node_features + '_'+ str(config.num_proteins)+'_fold_'+str(fold)+'_model.model'

    sys.stdout.flush()

    ret = None
    for epoch in range(1, config.num_epochs + 1):
        loss = train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch,
                     weight_dict=weight_dict)
        print('Train loss:', loss)

        if epoch % 10 == 0:
            print('Predicting for validation data...')
            file = '../results/drughub_drug_split_' + config.arch + '_prots_' + str(epoch) + '_epochs'
            with open(file=file, mode='a') as f:
                train_labels, train_predictions = predicting(model, device, train_loader)
                print('Train:', config.neg_sample_ratio,'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(train_labels, train_predictions),
                      dti_utils.dti_auroc(train_labels, train_predictions),
                      dti_utils.dti_f1_score(train_labels, train_predictions),
                      metrics.matthews_corrcoef(train_labels, train_predictions), file=f)

                test_labels, test_predictions = predicting(model, device, test_loader)
                print('Test:', config.neg_sample_ratio,'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(test_labels, test_predictions),
                      dti_utils.dti_auroc(test_labels, test_predictions),
                      dti_utils.dti_f1_score(test_labels, test_predictions),
                      metrics.matthews_corrcoef(test_labels, test_predictions), file=f)

                metrics_func_list = [metrics.accuracy_score, dti_utils.dti_auroc, dti_utils.dti_f1_score,
                                     metrics.matthews_corrcoef]
                ret = [list_fun(test_labels, test_predictions) for list_fun in metrics_func_list]

                test_loss = metrics.log_loss(test_labels, test_predictions, eps=0.000001)

        sys.stdout.flush()
    results.append(ret)

    if torch.cuda.is_available():
        torch.cuda.synchronize()

def inductive_missing_target_predictor(config):
    pass



if __name__ == '__main__':
    # Add parser arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_proteins", type=int, default=-1)
    parser.add_argument("--arch", type=str, default='GCNConv')
    parser.add_argument("--node_features", type=str, default='MolPred')
    parser.add_argument("--neg_sample_ratio", type=float, default=0.1)


    parser.add_argument("--num_epochs", type=int, default=3)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--num_folds", type=int, default=5)
    parser.add_argument("--lr", type=float, default=0.001)

    parser.add_argument("--fold", type=int, default=-1)

    parser.add_argument("--mode", type=str, default='standard')
    parser.add_argument("--drug_mode", type=str, default='standard')

    config = parser.parse_args()


    # Run classifier
    if config.mode == 'standard':
        transductive_missing_target_predictor(config)
    elif config.mode == 'drug':
        transductive_missing_drug_predictor(config)
    elif config.mode == 'protein_drughub':
        test_predictor_on_drughub_protein_data(config)
    elif config.mode == 'drug_drughub':
        test_predictor_on_drughub_drug_data(config)

