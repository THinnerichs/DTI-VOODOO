import numpy as np
import networkx as nx
import pandas as pd

from tqdm import tqdm

import DTI_data_preparation

from sklearn.model_selection import KFold
from sklearn import metrics

import matplotlib.pyplot as plt

import tensorflow as tf

import keras.backend as K
from keras import layers, regularizers, optimizers, losses, models
from keras.utils import plot_model


from keras_dgl.layers.graph_cnn_layer import GraphCNN
from keras_dgl.layers.graph_attention_cnn_layer import GraphAttentionCNN
from keras_dgl.layers.multi_graph_cnn_layer import MultiGraphCNN
from keras_dgl.layers.multi_graph_attention_cnn_layer import MultiGraphAttentionCNN
from keras_dgl.utils import *

import stellargraph as sg
from stellargraph import globalvar
from stellargraph.mapper import GraphSAGENodeGenerator, \
    Attri2VecNodeGenerator, Attri2VecLinkGenerator, \
    HinSAGENodeGenerator, \
    FullBatchNodeGenerator
from stellargraph.layer import GraphSAGE, Attri2Vec, HinSAGE, GCN, link_classification, GAT
from stellargraph.layer.ppnp import PPNP
from stellargraph.layer.appnp import APPNP
from stellargraph.data import UnsupervisedSampler

import dti_utils


def better_missing_target_predictor(results_filename = '../results/results_log',
                                    nb_epochs = 3,
                                    plot = False,
                                    embedding_layer_sizes = [32, 64],
                                    embedding_method='gcn',
                                    ):
    results_table_filename = "../results/results_table"

    print("Loading data ...")
    drug_list = np.array(DTI_data_preparation.get_drug_list())
    print("Get protein list ...")
    protein_list = np.array(DTI_data_preparation.get_human_proteins())
    print("Finished.\n")

    print("Scaling data ...")
    # side_effect_features = np.tile(DTI_data_preparation.get_side_effect_similarity_feature_list(drug_list),
                                   # (len(protein_list), 1))
    # side_effect_features = side_effect_features.reshape((len(protein_list), len(drug_list), 1430))
    # DDI_features = np.tile(DTI_data_preparation.get_DDI_feature_list(drug_list), (len(protein_list),1))
    DDI_features = DTI_data_preparation.get_DDI_feature_list(drug_list)
    # DDI_features = DDI_features.reshape((len(drug_list), len(protein_list), len(drug_list)))
    print("Finished.\n")

    print("Get DTIs")
    PPI_dti_features = DTI_data_preparation.get_PPI_dti_feature_list(drug_list, protein_list)

    y_dti_data = DTI_data_preparation.get_DTIs(drug_list=drug_list, protein_list=protein_list)
    y_dti_data = y_dti_data.reshape((len(protein_list), len(drug_list)))
    y_dti_data = np.transpose(y_dti_data)
    print("Finished.\n")
    print("Finished loading data.\n")

    # building stellar graph
    print("Building graph data ...")
    PPI_graph = DTI_data_preparation.get_PPI_DTI_graph_intersection()
    PPI_adj_mat = nx.adjacency_matrix(PPI_graph, weight=None).todense()

    # Build Graph Convolution filters
    SYM_NORM = True
    A_norm = preprocess_adj_numpy(PPI_adj_mat, SYM_NORM)
    num_filters = 2
    graph_conv_filters = np.concatenate([A_norm, np.matmul(A_norm, A_norm)], axis=0)
    graph_conv_filters = K.constant(graph_conv_filters)
    print("Adj_mat shape", PPI_adj_mat.shape)

    node_feature_mat = DTI_data_preparation.get_PPI_node_feature_mat_list(protein_list=protein_list)

    print("node_feature shape", node_feature_mat.shape)
    print("Finished.\n")

    # Fold parameters
    skf = KFold(n_splits=5, random_state=42)

    cv_scores = {'acc': [],
                 'auroc': [],
                 'f1-score': []}

    conf_matrix = None
    model = None
    round = 1
    print("Starting folds ...")
    for train, test in skf.split(protein_list):
        print("Round", round)
        round += 1

        train_filter_mask = np.zeros(len(protein_list))
        train_filter_mask[train] = 1

        y_graph_dti_train_data = np.zeros(y_dti_data.shape)
        y_graph_dti_train_data[:,train] = y_dti_data[:,train]

        y_graph_dti_test_data = np.zeros(y_dti_data.shape)
        y_graph_dti_test_data[:,test] = y_dti_data[:,test]



        # Build actual dti model
        GCN_layer_sizes = embedding_layer_sizes[:]
        PPI_input = layers.Input(shape=(node_feature_mat.shape[1],))
        graph_layer = GraphCNN(GCN_layer_sizes.pop(0), num_filters, graph_conv_filters, kernel_regularizer=regularizers.l2(5e-1),
                               activation='elu')(PPI_input)
        graph_layer = layers.Dropout(0.2)(graph_layer)

        for size in GCN_layer_sizes:
            graph_layer = GraphCNN(size, num_filters, graph_conv_filters, kernel_regularizer=regularizers.l2(5e-2),
                                   activation='elu')(graph_layer)

            graph_layer = layers.Dropout(0.2)(graph_layer)

        DDI_input = layers.Input(shape=(len(drug_list),))

        # side_effect_input = layers.Input(shape=(side_effect_features.shape[2],))

        merge_1 = layers.Concatenate(axis=1)([graph_layer,
                                              DDI_input,
                                              # side_effect_input
                                              ])

        dense_1 = layers.Dense(100, activation='relu')(merge_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        dense_1 = layers.Dense(20, activation='relu')(dropout_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        output = layers.Dense(1, activation='sigmoid')(dropout_1)

        model = models.Model(inputs=[PPI_input,
                                     DDI_input,
                                     # side_effect_input
                                     ],
                             outputs=output)

        model.compile(loss=losses.binary_crossentropy,
                      optimizer=optimizers.Adam(lr=0.005),
                      metrics=[dti_utils.dti_f1_score,
                               # dti_utils.dti_auroc,
                               'accuracy'])

        model.summary()

        imb_ratio = (np.prod(y_graph_dti_train_data.shape) - y_graph_dti_train_data.sum()) / y_graph_dti_train_data.sum()

        print('imb_ratio:', imb_ratio)

        class_weight = {0: 1.,
                        1: imb_ratio}

        print("Starting training ...")
        for epoch in range(nb_epochs):
            results = np.array([])

            for j in tqdm(range(len(drug_list))):
                epoch_DDI_feature = np.repeat([DDI_features[j,:]], len(protein_list), axis=0)
                # epoch_side_effect_feature = np.repeat()

                epoch_y_train = y_graph_dti_train_data[j, :]

                model.fit([node_feature_mat,
                           epoch_DDI_feature
                           ],
                          epoch_y_train,
                          sample_weight=train_filter_mask,
                          batch_size=len(protein_list),
                          epochs=1,
                          shuffle=False,
                          class_weight=class_weight,
                          verbose=0)

                epoch_y_pred = model.predict([node_feature_mat,
                                              epoch_DDI_feature
                                              ],
                                             batch_size=len(protein_list))

                epoch_y_pred = (epoch_y_pred.reshape((epoch_y_pred.shape[0])) >= 0.5).astype(int)
                results = np.vstack([results, epoch_y_pred]) if results.size else epoch_y_pred

            epoch_y_train_true = y_dti_data[:, train].flatten()
            epoch_y_train_pred = results[:, train].flatten()

            epoch_y_test_true = y_dti_data[:, test].flatten()
            epoch_y_test_pred = results[:, test].flatten()

            train_f1 = metrics.f1_score(y_true=epoch_y_train_true, y_pred=epoch_y_train_pred)
            test_f1 = metrics.f1_score(y_true=epoch_y_test_true, y_pred=epoch_y_test_pred)

            train_acc = metrics.accuracy_score(y_true=epoch_y_train_true, y_pred=epoch_y_train_pred)
            test_acc = metrics.accuracy_score(y_true=epoch_y_test_true, y_pred=epoch_y_test_pred)

            train_auc = metrics.roc_auc_score(y_true=epoch_y_train_true, y_score=epoch_y_train_pred)
            test_auc = metrics.roc_auc_score(y_true=epoch_y_test_true, y_score=epoch_y_test_pred)

            train_mcc = metrics.matthews_corrcoef(y_true=epoch_y_train_true, y_pred=epoch_y_train_pred)
            test_mcc = metrics.matthews_corrcoef(y_true=epoch_y_test_true, y_pred=epoch_y_test_pred)

            print("Epoch: {:04d}".format(epoch),
                  "train_acc= {:.4f}".format(train_acc), "test_acc= {:.4f}".format(test_acc),
                  "train_f1= {:.4f}".format(train_f1), "test_f1= {:.4f}".format(test_f1),
                  "train_auc= {:.4f}".format(train_auc), "test_auc= {:.4f}".format(test_auc),
                  "train_mcc= {:.4f}".format(train_mcc), "test_mcc= {:.4f}".format(test_mcc),
                  )

            if epoch == nb_epochs-1:

                conf_matrix = metrics.confusion_matrix(y_true=epoch_y_test_true,
                                                       y_pred=(epoch_y_test_pred.reshape((len(epoch_y_test_pred))) >= 0.5).astype(int))

                print("Confusion matrix", conf_matrix)

                tn = conf_matrix[0, 0]
                tp = conf_matrix[1, 1]
                fp = conf_matrix[0, 1]
                fn = conf_matrix[1, 0]

                precision = tp / (tp + fp)
                recall = tp / (tp + fn)
                accuracy = (tp + tn) / (tp + tn + fp + fn)

                auroc = metrics.roc_auc_score(epoch_y_test_true, epoch_y_test_pred)
                f1_score = metrics.f1_score(epoch_y_test_true, epoch_y_test_pred)

                print("accuracy", accuracy * 100)
                print("prec", precision * 100)
                print("recall", recall * 100)
                print("auroc", auroc * 100)
                print("f1-score", f1_score * 100)

                cv_scores['acc'].append(accuracy * 100)
                cv_scores['auroc'].append(auroc * 100)
                cv_scores['f1-score'].append(f1_score * 100)

                print("Mean accuracy:", np.mean(cv_scores['acc']))
                print("Mean auroc:", np.mean(cv_scores['auroc']))
                print("Mean f1-score:", np.mean(cv_scores['f1-score']))

        with open(file=results_table_filename, mode='a') as f:
            print("1" + "\t" + "1" + "\t" + "0" + "\t" + "0" + "\t" +
                  str(len(protein_list)) + "\t" + str(len(drug_list)) + "\t" +
                  str(nb_epochs) + "\t" +
                  str(embedding_layer_sizes) + "\t" +
                  'better_' + embedding_method + '\t' +
                  str(np.mean(cv_scores['acc'])) + "\t" + str(np.mean(cv_scores['auroc'])) + "\t" + str(
                np.mean(cv_scores['f1-score'])),
                  file=f)
        with open(file=results_filename, mode='a') as filehandler:
            print("BETTER DTI " + embedding_method + " PREDICTION", file=filehandler)

            print("Including:", file=filehandler)
            print("- PPIs", file=filehandler)
            # print("- DDIs", file=filehandler)
            print("- side similarity scores", file=filehandler)
            # print("- HMM node features", file=filehandler)
            print("Number of targets:\t", len(protein_list), file=filehandler)
            print("Number of drugs:\t", len(drug_list), file=filehandler)

            print("Mean accuracy: {},\t Std: {}".format(np.mean(cv_scores['acc']), np.std(cv_scores['acc'])), file=filehandler)
            print("Mean auroc: {},\t Std: {}".format(np.mean(cv_scores['auroc']), np.std(cv_scores['auroc'])), file=filehandler)
            print("Mean f1: {},\t Std: {}".format(np.mean(cv_scores['f1-score']), np.std(cv_scores['f1-score'])),
                  file=filehandler)

            print("Epochs: {} ".format(nb_epochs), file=filehandler)
            model.summary(print_fn=lambda x: filehandler.write(x + '\n'))

            # print confusion matrix
            print("Confusion matrix:",
                  conf_matrix,
                  file=filehandler)

            print("------------------------------------------------\n")

            # serialize model to JSON
            if plot:
                tf.keras.utils.plot_model(model,
                                          show_shapes=True,
                                          to_file='../models/better_gcn_model.png')


if __name__ == '__main__':

    better_missing_target_predictor(nb_epochs=20,
                                    plot=False,
                                    embedding_layer_sizes=[16]
                                    )


