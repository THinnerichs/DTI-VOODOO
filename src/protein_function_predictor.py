import numpy as np
import math
from sklearn.model_selection import KFold
from sklearn import metrics

import torch

import pickle
import sys
import gc
import argparse

from molecular_utils import train, predicting
import dti_utils
from protein_function_utils import *


def protein_function_predictor(config):
    model_st = 'molecular_predictor'

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    if config.num_proteins <= 0:
        config.num_proteins = None

    dti_data = ProteinFunctionDTIDataBuilder(num_proteins=config.num_proteins)

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
        gc.collect()

        # build train data over whole dataset with help matrix
        train_indices = help_matrix[:, train_protein_indices].flatten()
        test_indices = help_matrix[:, test_protein_indices].flatten()
        print(train_indices.shape, test_indices.shape)

        train_dataset = dti_data.get(train_indices)
        test_dataset = dti_data.get(test_indices)

        train_dataset = ProteinFunctionDTIDataset(train_dataset)
        test_dataset = ProteinFunctionDTIDataset(test_dataset)

        # Calculate weights
        positives = dti_data.y_dti_data.flatten()[train_indices].sum()
        len_to_sum_ratio = (len(train_indices) - positives) / positives
        weight_dict = {0: 1.,
                       1: len_to_sum_ratio}

        train_loader = data.DataLoader(train_dataset, batch_size=config.batch_size, shuffle=True)
        test_loader = data.DataLoader(test_dataset, batch_size=config.batch_size)

        model = ProteinFunctionPredNet()
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
            loss = train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch,
                  weight_dict=weight_dict)
            print('Train Loss:', loss)

            if epoch%config.num_epochs == 0:
                print('Predicting for validation data...')
                file='../results/protfunc_pred_results_' + str(config.num_epochs)+'_epochs'
                with open(file=file, mode='a') as f:
                    train_labels, train_predictions = predicting(model, device, train_loader)
                    print('Train:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(train_labels, train_predictions),
                          dti_utils.dti_auroc(train_labels, train_predictions),
                          dti_utils.dti_f1_score(train_labels, train_predictions),
                          metrics.matthews_corrcoef(train_labels, train_predictions), file=f)

                    test_labels, test_predictions = predicting(model, device, test_loader)
                    print('Test:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
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

        print('Build overall data...')
        sys.stdout.flush()
        gc.collect()
        overall_dataset = dti_data.get(np.arange(dti_data.num_drugs * dti_data.num_proteins))
        overall_dataloader = data.DataLoader(overall_dataset, batch_size=config.batch_size)

        print('Predicting...')
        sys.stdout.flush()
        labels, predictions = predicting(model, device, overall_dataloader)
        filename = '../models/protein_function_predictor/pred_fold_'+str(fold)
        with open(file=filename+'.pkl', mode='wb') as f:
            pickle.dump(predictions, f, pickle.HIGHEST_PROTOCOL)

        gc.collect()

        model_filename = '../models/protein_function_predictor/mol_pred_'+ (config.model_id +'_' if config.model_id else '') + 'model_fold_'+str(fold)+'.model'
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
    parser.add_argument("--fold", type=int, default=-1)

    config = parser.parse_args()

    protein_function_predictor(config)


