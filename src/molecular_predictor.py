import numpy as np
import math

from sklearn.model_selection import KFold
from sklearn import metrics

from xgboost import XGBClassifier


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
from tqdm import tqdm
import gc



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
    model_st = 'molecular_predictor'

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
        gc.collect()

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

            if False and epoch%10 == 0:
                print('Predicting for validation data...')
                train_labels, train_predictions = predicting(model, device, train_loader)
                print('Train:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(train_labels, train_predictions),
                      dti_utils.dti_auroc(train_labels, train_predictions),
                      dti_utils.dti_f1_score(train_labels, train_predictions),
                      metrics.matthews_corrcoef(train_labels, train_predictions))

                test_labels, test_predictions = predicting(model, device, test_loader)
                print('Test:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(test_labels, test_predictions),
                      dti_utils.dti_auroc(test_labels, test_predictions),
                      dti_utils.dti_f1_score(test_labels, test_predictions),
                      metrics.matthews_corrcoef(test_labels, test_predictions))

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
        filename = '../models/molecular_predictor/pred_fold_'+str(fold)
        with open(file=filename+'.pkl', mode='wb') as f:
            pickle.dump(predictions, f, pickle.HIGHEST_PROTOCOL)

        gc.collect()

        model_filename = '../models/molecular_predictor/mol_pred_'+ (config.model_id +'_' if config.model_id else '') + 'model_fold_'+str(fold)+'.model'
        torch.save(model.state_dict(), model_filename)
        print("Done.")
        sys.stdout.flush()

        results.append(ret)

    results = np.array(results)
    results = [(results[:, i].mean(), results[:, i].std()) for i in range(results.shape[1])]

    results_file_name = '../results/molecular_model' + '_'+ (config.model_id +'_' if config.model_id else '') + str(config.num_proteins) + '_model_results'

    print('Overall Results:')
    print('Model\tacc\tauroc\tf1\tmatt')
    print(model_st + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)))

    with open(results_file_name, 'a') as f:
        print('Model\tacc\tauroc\tf1\tmatt', file=f, end='\n')
        print(model_st + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)), file=f, end='\n')

    print("Done.")
    sys.stdout.flush()

def drug_split_molecular_predictor(config):
    model_st = 'molecular_drug_split_predictor'

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    dti_data = MolecularDTIDataBuilder(num_proteins=None)

    # generate indices for proteins
    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    X = np.zeros((dti_data.num_drugs, 1))

    # build for help matrix for indices
    help_matrix = np.arange(dti_data.num_drugs * dti_data.num_proteins)
    help_matrix = help_matrix.reshape((dti_data.num_drugs, dti_data.num_proteins))

    results = []
    fold = 0
    for train_drug_indices, test_drug_indices in kf.split(X):
        fold += 1
        print("Fold:", fold)

        # build train data over whole dataset with help matrix
        train_indices = help_matrix[train_drug_indices, :].flatten()
        test_indices = help_matrix[test_drug_indices, :].flatten()
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

            if epoch%10 == 0:
                print('Predicting for validation data...')
                train_labels, train_predictions = predicting(model, device, train_loader)
                print('Train:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(train_labels, train_predictions),
                      dti_utils.dti_auroc(train_labels, train_predictions),
                      dti_utils.dti_f1_score(train_labels, train_predictions),
                      metrics.matthews_corrcoef(train_labels, train_predictions))

                test_labels, test_predictions = predicting(model, device, test_loader)
                print('Test:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
                      metrics.accuracy_score(test_labels, test_predictions),
                      dti_utils.dti_auroc(test_labels, test_predictions),
                      dti_utils.dti_f1_score(test_labels, test_predictions),
                      metrics.matthews_corrcoef(test_labels, test_predictions))

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

        overall_dataset = dti_data.get(np.arange(dti_data.num_drugs * dti_data.num_proteins))
        overall_loader = data.DataLoader(overall_dataset, batch_size=config.batch_size)

        labels, predictions = predicting(model, device, overall_loader)
        filename = '../models/molecular_drug_split_predictor/pred_fold_'+str(fold)
        with open(file=filename+'.pkl', mode='wb') as f:
            pickle.dump(predictions, f, pickle.HIGHEST_PROTOCOL)

        model_filename = '../models/molecular_drug_split_predictor/mol_pred_'+ (config.model_id +'_' if config.model_id else '') + 'model_fold_'+str(fold)+'.model'
        torch.save(model.state_dict(), model_filename)

        results.append(ret)

    results = np.array(results)
    results = [(results[:, i].mean(), results[:, i].std()) for i in range(results.shape[1])]

    results_file_name = '../results/molecular_model' + '_'+ (config.model_id +'_' if config.model_id else '') + str(config.num_proteins) + '_model_results'

    print('Overall Results:')
    print('Model\tacc\tauroc\tf1\tmatt')
    print(model_st + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)))

    with open(results_file_name, 'a') as f:
        print('Model\tacc\tauroc\tf1\tmatt', file=f, end='\n')
        print(model_st + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)), file=f, end='\n')

    print("Done.")
    sys.stdout.flush()

def XGBoost_molecular_predictor(config):
    model_st = 'molecular_xgboost_predictor'

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

        print("Building train data...")
        X_train = []
        Y_train = []
        for x,y in tqdm(train_dataset):
            X_train.append(x.numpy())
            Y_train.append(y)
        X_train = np.array(X_train)
        Y_train = np.array(Y_train)
        print("Building test data...")
        X_test = []
        Y_test = []
        for x,y in tqdm(test_dataset):
            X_test.append(x.numpy())
            Y_test.append(y)
        X_test = np.array(X_test)
        Y_test = np.array(Y_test)

        # Calculate weights
        positives = dti_data.y_dti_data.flatten()[train_indices].sum()
        len_to_sum_ratio = (len(train_indices)-positives)/positives # negatives/positives
        weight_dict = {0: 1.,
                       1: len_to_sum_ratio}

        model = XGBClassifier(scale_pos_weight=len_to_sum_ratio,
                              max_depth=4,
                              n_jobs=40,
                              tree_method='gpu_hist')
        model.fit(X_train, Y_train)
        y_pred = np.around(model.predict(X_test))

        sys.stdout.flush()

        test_labels = Y_test
        test_predictions = y_pred

        print('Test:', 'Acc, ROC_AUC, f1, matthews_corrcoef',
              metrics.accuracy_score(test_labels, test_predictions),
              dti_utils.dti_auroc(test_labels, test_predictions),
              dti_utils.dti_f1_score(test_labels, test_predictions),
              metrics.matthews_corrcoef(test_labels, test_predictions))

        metrics_func_list = [metrics.accuracy_score, dti_utils.dti_auroc, dti_utils.dti_f1_score,
                             metrics.matthews_corrcoef]
        ret = [list_fun(test_labels, test_predictions) for list_fun in metrics_func_list]

        test_loss = metrics.log_loss(test_labels, test_predictions, eps=0.000001)
        print(test_loss)
        sys.stdout.flush()

        overall_dataset = dti_data.get(np.arange(dti_data.num_drugs * dti_data.num_proteins))
        overall_data = np.array([x.numpy() for x,y in overall_dataset])

        predictions = np.around(model.predict(overall_data))
        filename = '../models/molecular_xgboost_predictor/pred_fold_'+str(fold)
        with open(file=filename+'.pkl', mode='wb') as f:
            pickle.dump(predictions, f, pickle.HIGHEST_PROTOCOL)

        # model_filename = '../models/molecular_xgboost_predictor/mol_pred_'+ (config.model_id +'_' if config.model_id else '') + 'model_fold_'+str(fold)+'.model'
        # torch.save(model.state_dict(), model_filename)

        results.append(ret)

    results = np.array(results)
    results = [(results[:, i].mean(), results[:, i].std()) for i in range(results.shape[1])]

    results_file_name = '../results/molecular_xgboost_model' + '_'+ (config.model_id +'_' if config.model_id else '') + str(config.num_proteins) + '_model_results'

    print('Overall Results:')
    print('Model\tacc\tauroc\tf1\tmatt')
    print(model_st + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)))

    with open(results_file_name, 'a') as f:
        print('Model\tacc\tauroc\tf1\tmatt', file=f, end='\n')
        print(model_st + '\t' + str(config.num_proteins) + '\t' + '\t'.join(map(str, results)), file=f, end='\n')

    print("Done.")


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

    parser.add_argument("--model_id", type=str, default='')
    parser.add_argument("--model", type=str, default='protein')

    config = parser.parse_args()

    # Run classifier
    if config.model == 'protein':
        molecular_predictor(config)
    elif config.model == 'drug':
        drug_split_molecular_predictor(config)
    elif config.model == 'xgboost':
        XGBoost_molecular_predictor(config)
    else:
        print("No valid model selected.")
        raise ValueError

