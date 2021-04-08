import numpy as np
import math
from sklearn.model_selection import KFold
from sklearn import metrics

import gensim

from tqdm import tqdm
import argparse
import sys
import pickle

import torch
import torch.nn as nn
import torch.utils.data as data

import DTI_data_preparation
import PPI_utils
import PhenomeNET_DL2vec_utils
from molecular_utils import train, predicting
import dti_utils


class HPODTIDataBuilder:
    def __init__(self, config):
        self.config = config

        print("Loading data ...")
        if config.yamanishi_test or config.biosnap_test:
            self.drug_list = np.array(PhenomeNET_DL2vec_utils.get_PhenomeNET_drug_list())
            if config.include_indications:
                indication_drug_list = set(PhenomeNET_DL2vec_utils.get_UMLS_drug_list())
                self.drug_list = np.array(list(set(self.drug_list) & indication_drug_list))
        else:
            self.drug_list = np.array(DTI_data_preparation.get_drug_list(config.mode))
        print(len(self.drug_list), "drugs present")

        # get protein lists for each ontology
        uberon_protein_list = PhenomeNET_DL2vec_utils.get_PhenomeNET_protein_list(mode='uberon')
        GO_protein_list = PhenomeNET_DL2vec_utils.get_PhenomeNET_protein_list(mode='GO')
        MP_protein_list = PhenomeNET_DL2vec_utils.get_PhenomeNET_protein_list(mode='MP')

        dti_graph = DTI_data_preparation.get_human_DTI_graph(mode=config.mode)
        PPI_graph = PPI_utils.get_PPI_graph(min_score=700)

        if config.yamanishi_test or config.biosnap_test:
            self.protein_list = np.array(list(set(PPI_graph.nodes()) & (set(uberon_protein_list) | set(GO_protein_list) | set(MP_protein_list))))
            if config.yamanishi_test:
                print("Loading Yamanishi data ...")
                self.drug_list, self.protein_list, self.y_dti_data = DTI_data_preparation.get_yamanishi_data(self.drug_list, self.protein_list)
            elif config.biosnap_test:
                self.drug_list, self.protein_list, self.y_dti_data = DTI_data_preparation.get_BioSnap_data(self.drug_list, self.protein_list)

            print(self.drug_list.shape, self.y_dti_data.shape, self.protein_list.shape)
        else:
            if config.include_indications:
                print("Drug indications can only be selected with yamanishi test for now.")
                raise ValueError
        print(len(self.protein_list), "proteins present\n")

        # PPI data
        print("Loading PPI graph ...")
        self.PPI_graph = DTI_data_preparation.get_PPI_DTI_graph_intersection()
        self.PPI_graph = self.PPI_graph.subgraph(self.protein_list)

        # calculate dimensions of network
        self.num_proteins = len(self.protein_list)
        self.num_drugs = len(self.drug_list)

        config.num_proteins = self.num_proteins
        config.num_drugs = self.num_drugs


    def build_data(self, config):

        # DTI data
        if not config.yamanishi_test and not config.biosnap_test:
            print("Loading DTI links ...")
            y_dti_data = DTI_data_preparation.get_DTIs(drug_list=self.drug_list, protein_list=self.protein_list,
                                                       mode=config.mode)
            self.y_dti_data = y_dti_data.reshape((len(self.drug_list), len(self.protein_list)))
        print(self.y_dti_data.shape, self.y_dti_data.sum())

        self.feature_matrix = np.zeros((self.num_drugs, self.num_proteins))
        epsilon = 0.00001

        DL2vec_path_prefix = '../data/PhenomeNET_data/'
        drug_indication_prefix = '../data/drug_indications/'

        drug_model_filename = DL2vec_path_prefix + 'drug_embedding_model'
        drug_indication_filename = drug_indication_prefix + 'embedding_model'
        uberon_model_filename = DL2vec_path_prefix + 'uberon_embedding_model'
        GO_model_filename = DL2vec_path_prefix + 'GO_embedding_model'
        MP_model_filename = DL2vec_path_prefix + 'MP_embedding_model'

        # load models
        drug_model = gensim.models.Word2Vec.load(drug_model_filename)
        drug_indication_model = gensim.models.Word2Vec.load(drug_indication_filename)
        uberon_model = gensim.models.Word2Vec.load(uberon_model_filename)
        GO_model = gensim.models.Word2Vec.load(GO_model_filename)
        MP_model = gensim.models.Word2Vec.load(MP_model_filename)

        # Build wordvector dicts
        drug_model = drug_model.wv
        drug_indication_model = drug_indication_model.wv
        uberon_model = uberon_model.wv
        GO_model = GO_model.wv
        MP_model = MP_model.wv

        drug_embeddings = []
        drug_indication_embeddings = []
        uberon_embeddings = []
        GO_embeddings = []
        MP_embeddings = []
        for protein in self.protein_list:
            protein_id = protein

            if protein_id in uberon_model.vocab.keys():
                uberon_embeddings.append(uberon_model[protein_id])
            else:
                uberon_embeddings.append(torch.zeros((200)))
            if protein_id in GO_model.vocab.keys():
                GO_embeddings.append(GO_model[protein_id])
            else:
                GO_embeddings.append(torch.zeros((200)))
            if protein_id in MP_model.vocab.keys():
                MP_embeddings.append(MP_model[protein_id])
            else:
                GO_embeddings.append(torch.zeros((200)))

        for drug_id in self.drug_list:
            drug_embeddings.append((drug_model[drug_id]))

            if config.include_indications:
                drug_indication_embeddings.append((drug_indication_model[drug_id]))

        self.drug_embeddings = torch.Tensor(drug_embeddings)
        self.uberon_embeddings = torch.Tensor(uberon_embeddings)
        self.GO_embeddings = torch.Tensor(GO_embeddings)
        self.MP_embeddings = torch.Tensor(MP_embeddings)

        if config.include_indications:
            self.drug_indication_embeddings = torch.Tensor(drug_indication_embeddings)

        self.protein_embeddings = torch.cat([self.uberon_embeddings, self.GO_embeddings, self.MP_embeddings], dim=1)

        print("Finished.\n")

    def get(self, indices):
        data_list = []

        # indices = list(range(self.num_drugs))
        for index in tqdm(indices):

            drug_index = index // self.num_proteins
            protein_index = index % self.num_proteins

            # build protein mask

            y = int(self.y_dti_data[drug_index, protein_index])

            if self.config.include_indications:
                drug_feature = torch.cat([self.drug_embeddings[drug_index, :], self.drug_indication_embeddings[drug_index, :]])
                # drug_feature = self.drug_indication_embeddings[drug_index, :]
            else:
                drug_feature = self.drug_embeddings[drug_index, :]
            protein_feature = self.protein_embeddings[protein_index, :]

            # additional
            # degree_feature = torch.tensor(self.degree_features[protein_index, :])

            data_list.append((torch.cat((drug_feature, protein_feature), 0).float(), y))

        return data_list

    def __len__(self):
        return self.num_drugs * self.num_proteins

class DTIGraphDataset(data.Dataset):
    def __init__(self, data_list):
        super(DTIGraphDataset, self).__init__()
        self.data_list = data_list

    def __getitem__(self, idx):
        return self.data_list[idx]

    def __len__(self):
        return len(self.data_list)

    def _download(self):
        pass

    def _process(self):
        pass


class HPOPredNet(nn.Module):
    def __init__(self, include_indications=False):
        super(HPOPredNet, self).__init__()

        self.dropout = nn.Dropout(0.5)

        self.include_indications = include_indications

        # siamese network approach
        if self.include_indications:
            self.model = nn.Sequential(
                nn.Linear(400, 128),
                nn.Dropout(0.5),
                nn.BatchNorm1d(128, affine=True),
                nn.LeakyReLU(0.2, inplace=True),
                nn.Linear(128, 200),
                # nn.Dropout(0.5),
                # nn.BatchNorm1d(200),
                # nn.LeakyReLU(0.2, inplace=True),
                # nn.Linear(200, 200),
                # nn.Sigmoid()
            )
        else:
            self.model = nn.Sequential(
                nn.Linear(200, 256),
                # nn.Dropout(0.2),
                nn.BatchNorm1d(256),
                nn.LeakyReLU(0.2, inplace=True),
                nn.Linear(256, 200),
                # nn.Dropout(0.5),
                # nn.BatchNorm1d(50),
                # nn.LeakyReLU(0.2, inplace=True),
                # nn.Linear(256, 1),
                # nn.Sigmoid()
            )
        self.model2 = nn.Sequential(
            nn.Linear(600, 256),
            # nn.Dropout(0.5),
            nn.BatchNorm1d(256, affine=True),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(256, 200),
            # nn.BatchNorm1d(200),
            # nn.Dropout(0.5),
            # nn.LeakyReLU(0.2, inplace=True),
            # nn.Linear(200, 200),
            # nn.Sigmoid()
        )

        self.sim = nn.CosineSimilarity(dim=1)

    def forward(self, x):


        if self.include_indications:
            p1 = self.model(x[:,:400]).view(-1, 200)
            d1 = self.model2(x[:,400:]).view(-1, 200)
        else:
            p1 = self.model(x[:,:200]).view(-1, 200)
            d1 = self.model2(x[:,200:]).view(-1, 200)

        s1 = self.sim(p1, d1)

        out = s1.reshape(-1, 1)

        # out = self.output_sig(s1)

        return out.sigmoid()


def siamese_drug_protein_network(config):
    model_st = 'siamese_drug_protein_network'

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    if config.num_proteins <= 0:
        config.num_proteins = None

    dti_data = HPODTIDataBuilder(config)

    # generate indices for proteins
    kf = KFold(n_splits=config.num_folds, random_state=42, shuffle=True)
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

        print('len(train_dataset)', len(train_dataset))

        # Calculate weights
        positives = dti_data.y_dti_data.flatten()[train_indices].sum()
        neg_to_pos_ratio = (len(train_indices) - positives) / positives

        train_loader = data.DataLoader(train_dataset, batch_size=config.batch_size, shuffle=True)
        test_loader = data.DataLoader(test_dataset, batch_size=config.batch_size)

        model = HPOPredNet(include_indications=config.include_indications)
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
        best_AUROC = 0
        for epoch in range(1, config.num_epochs + 1):
            loss = train(config=config, model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch, neg_to_pos_ratio=neg_to_pos_ratio)
            print('Train Loss:', loss)

            if epoch%2 == 0:
                print('Predicting for validation data...')
                file='../results/HPO_pred_results_' + str(config.num_epochs)+'_epochs'
                with open(file=file, mode='a') as f:
                    train_labels, train_predictions = predicting(model, device, train_loader)
                    print('Train: Acc, ROC_AUC, MicroAUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(train_labels, train_predictions.round()),
                          dti_utils.dti_auroc(train_labels, train_predictions),
                          dti_utils.micro_AUC_per_prot(train_labels, train_predictions, config.num_drugs),
                          dti_utils.dti_f1_score(train_labels, train_predictions.round()),
                          dti_utils.dti_mcc(train_labels, train_predictions.round()), file=f)

                    test_labels, test_predictions = predicting(model, device, test_loader)
                    print('Test: Acc, ROC_AUC, MicroAUC, f1, matthews_corrcoef',
                          metrics.accuracy_score(test_labels, test_predictions.round()),
                          dti_utils.dti_auroc(test_labels, test_predictions),
                          dti_utils.micro_AUC_per_prot(test_labels, test_predictions, config.num_drugs),
                          dti_utils.dti_f1_score(test_labels, test_predictions.round()),
                          dti_utils.dti_mcc(test_labels, test_predictions.round()), file=f)


                    test_AUROC = dti_utils.micro_AUC_per_prot(test_labels, test_predictions, config.num_drugs)

                    if test_AUROC > best_AUROC:
                        state_dict_path = '../models/HPO_models/hpo_pred_fold_'+str(fold)+'_model'
                        torch.save(model.state_dict(), state_dict_path)
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
    parser.add_argument("--num_folds", type=int, default=5)
    parser.add_argument("--lr", type=float, default=0.001)

    parser.add_argument("--model_id", type=str, default='')
    parser.add_argument("--model", type=str, default='protein')
    parser.add_argument("--mode", type=str, default='')
    parser.add_argument("--fold", type=int, default=3)

    parser.add_argument("--yamanishi_test", action='store_true')
    parser.add_argument("--biosnap_test", action='store_true')
    parser.add_argument("--include_indications", action='store_true')

    config = parser.parse_args()

    if config.model=='protein':
        siamese_drug_protein_network(config)
    elif config.model=='drug':
        # drug_split_protein_function_predictor(config)
        pass
    else:
        print('No valid model selected...')
        raise ValueError

