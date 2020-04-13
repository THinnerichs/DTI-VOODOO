import numpy as np
import math

from sklearn.model_selection import KFold
from sklearn import metrics

import torch
import torch.utils as utils
import torch.nn as nn
import torch.utils.data as data

from smiles_transformer.pretrain_trfm import TrfmSeq2seq
from smiles_transformer.pretrain_rnn import RNNSeq2Seq
from smiles_transformer.build_vocab import WordVocab
from smiles_transformer.utils import split

import DTI_data_preparation
import DDI_utils
import dti_utils
from molecular_utils import *

from torch_dti_utils import rmse, mse, pearson, spearman, ci

import subprocess
import argparse
import sys



def write_encoded_drugs(drug_list,
                        mode='trfm'):

    drug_SMILES_dict = DTI_data_preparation.get_truncated_drug_to_SMILES_dict()

    pad_index = 0
    unk_index = 1
    eos_index = 2
    sos_index = 3
    mask_index = 4

    # Adapted from https://github.com/DSPsleeporg/smiles-transformer
    pretrained_model_dir = '../models/drug_representation/'
    vocab = WordVocab.load_vocab(pretrained_model_dir + 'vocab.pkl')

    def get_inputs(sm):
        seq_len = 614 # formerly 220
        sm = sm.split()
        if len(sm) > 612: #formerly 218
            print('SMILES is too long ({:d})'.format(len(sm)))
            sm = sm[:109] + sm[-109:]
        ids = [vocab.stoi.get(token, unk_index) for token in sm]
        ids = [sos_index] + ids + [eos_index]
        seg = [1] * len(ids)
        padding = [pad_index] * (seq_len - len(ids))
        ids.extend(padding), seg.extend(padding)
        return ids, seg

    def get_array(smiles):
        x_id, x_seg = [], []
        for sm in smiles:
            a, b = get_inputs(sm)
            x_id.append(a)
            x_seg.append(b)
        return torch.tensor(x_id), torch.tensor(x_seg)

    out_filename = '../models/drug_representation/'+mode+'.npy'
    if mode=='trfm':
        trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
        trfm.load_state_dict(torch.load(pretrained_model_dir+'trfm.pkl'))
        trfm.eval()

        x_split = [split(sm) for sm in [drug_SMILES_dict[drug] for drug in drug_list]]
        xid, xseg = get_array(x_split)
        X = trfm.encode(torch.t(xid))

        print('Trfm size:', X.shape)

        np.save(out_filename, X)
    elif mode=='rnn':
        rnn = RNNSeq2Seq(len(vocab), 256, len(vocab), 3)
        rnn.load_state_dict(torch.load(pretrained_model_dir+'seq2seq.pkl'))
        rnn.eval()

        x_split = [split(sm) for sm in [drug_SMILES_dict[drug] for drug in drug_list]]
        xid, _ = get_array(x_split)
        X = rnn.encode(torch.t(xid))

        print('RNN size:', X.shape)

        np.save(out_filename, X)
    else:
        print("No valid mode selected for drug to SMILES encoding.")
        raise ValueError

def write_encoded_proteins():

    protein_dir = "../models/protein_representation/"
    in_file = protein_dir+'data/PPI_graph_protein_seqs.fasta'
    out_file = protein_dir+'results/output'

    subprocess.call(protein_dir+'predict.sh {} {}'.format(in_file, out_file))

def molecular_predictor(config):

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    if config.num_proteins <= 0:
        config.num_proteins = None

    dti_data = MolecularDTIDataBuilder(num_proteins=config.num_proteins)

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

        # build train data over whole dataset with help matrix
        train_indices = help_matrix[:, train_protein_indices].flatten()
        test_indices = help_matrix[:, test_protein_indices].flatten()
        print(train_indices.shape, test_indices.shape)

        train_dataset = dti_data.get(train_indices)
        test_dataset = dti_data.get(test_indices)

        train_dataset = MolecularDTIDataset(train_dataset)
        test_dataset = MolecularDTIDataset(test_dataset)

        # Calculate weights
        len_to_sum_ratio = len(train_indices)/dti_data.y_dti_data.flatten()[train_indices].sum()
        weight_dict = {0: 1.,
                       1: len_to_sum_ratio}

        train_loader = data.DataLoader(train_dataset, batch_size=config.batch_size, shuffle=True)
        test_loader = data.DataLoader(test_dataset, batch_size=config.batch_size)

        model = MolecularPredNet()
        model = nn.DataParallel(model).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

        # storing best results
        best_loss = math.inf
        best_test_loss = math.inf
        best_epoch = -1
        best_test_ci = 0

        model_st = 'molecular_predictor'
        # model_file_name = '../models/' + model_st + '_' + config.node_features + '_' + str(
            # config.num_proteins) + '_fold_' + str(fold) + '_model.model'

        sys.stdout.flush()

        ret = None
        for epoch in range(1, config.num_epochs + 1):
            loss = train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch,
                  weight_dict=weight_dict)
            print('Train Loss:', loss)

            if epoch%10 == 0:
                print('Predicting for validation data...')
                labels, predictions = predicting(model, device, train_loader)
                predictions = np.around(predictions)
                print('Validation:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(labels, predictions),
                      dti_utils.dti_auroc(labels, predictions),
                      dti_utils.dti_f1_score(labels, predictions),
                      metrics.matthews_corrcoef(labels, predictions))

                G, P = predicting(model, device, test_loader)
                predictions = np.around(predictions)
                val = metrics.log_loss(labels, predictions, eps=0.000001)
                print('predictions: max/min', predictions.max(), predictions.min())
                if val < best_loss:
                    best_loss = val
                    best_epoch = epoch + 1
                    # torch.save(model.state_dict(), model_file_name)
                    print('predicting for test data')
                    # G, P = predicting(model, device, test_loader)
                    # ret = [rmse(G, P), mse(G, P), pearson(G, P), spearman(G, P), ci(G, P)]
                    ret = []
                    P = np.around(P)
                    metrics_func_list = [metrics.accuracy_score, dti_utils.dti_auroc, dti_utils.dti_f1_score,
                                         metrics.matthews_corrcoef]
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
                    print(ret[1], 'No improvement since epoch ', best_epoch, '; best_test_mse,best_test_ci:',
                          best_test_loss,
                          best_test_ci, model_st)
            sys.stdout.flush()
        
        results.append(ret)

    results = np.array(results)
    results = [(results[:, i].mean(), results[:, i].std()) for i in range(results.shape[1])]

    results_file_name = '../results/molecular_model' + '_' + str(config.num_proteins) + '_model_results'

    print('Overall Results:')
    print('Model\trmse\tmse\tpearson\tspearman\tacc\tauroc\tf1\tmatt')
    print(config.arch + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)))

    with open(results_file_name, 'a') as f:
        print('Model\trmse\tmse\tpearson\tspearman\tacc\tauroc\tf1\tmatt', file=f)
        print(config.arch + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)), file=f)

    print("Done.")
    sys.stdout.flush()


if __name__=='__main__':
    drug_list = DTI_data_preparation.get_drug_list()

    # write_encoded_drugs(drug_list, mode='trfm')
    # write_encoded_drugs(drug_list, mode='rnn')

    # write_encoded_proteins()

    # Add parser arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--num_proteins", type=int, default=-1)
    # parser.add_argument("--arch", type=str, default='GCNConv')
    # parser.add_argument("--node_features", type=str, default='simple')

    parser.add_argument("--num_epochs", type=int, default=20)
    parser.add_argument("--batch_size", type=int, default=1024)
    # parser.add_argument("--num_folds", type=int, default=5)
    parser.add_argument("--lr", type=float, default=0.001)

    config = parser.parse_args()

    # Run classifier
    molecular_predictor(config)
