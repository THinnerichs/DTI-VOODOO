import numpy as np
import networkx as nx
import pandas as pd

import DTI_data_preparation

from sklearn.model_selection import KFold
from sklearn import metrics

import matplotlib.pyplot as plt

import tensorflow as tf

from tensorflow.keras import models, layers, optimizers, losses, regularizers
# import tensorflow.keras.layers as layers
import tensorflow.keras.backend as K
# import tensorflow.keras.regularizers as regularizers
# import tensorflow.keras.optimizers as optimizers
# import tensorflow.keras.losses as losses
# from tensorflow.keras.utils import plot_model


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



def missing_target_predictor(results_filename='../results/results_log',
                             nb_epochs=3,
                             batch_size=500,
                             plot=False,
                             num_samples=[100,20],
                             embedding_layer_sizes=[32,64],
                             embedding_output_size=32,
                             embedding_method='graphsage',
                             supervised=True):
    results_table_filename = "../results/results_table"

    print("Loading data ...")
    drug_list = np.array(DTI_data_preparation.get_drug_list())
    print("Get protein list ...")
    protein_list = np.array(DTI_data_preparation.get_human_proteins())[:5000]
    print("Finished.\n")

    print("Scaling data ...")
    side_effect_features = np.tile(DTI_data_preparation.get_side_effect_similarity_feature_list(drug_list), (len(protein_list),1))
    side_effect_features = side_effect_features.reshape((len(protein_list), len(drug_list), 1430))
    # DDI_features = np.tile(DTI_data_preparation.get_DDI_feature_list(drug_list), (len(protein_list),1))
    # DDI_features = DDI_features.reshape((len(protein_list), len(drug_list), len(drug_list)))
    print("Finished.\n")

    print("Get DTIs")
    PPI_dti_features = DTI_data_preparation.get_PPI_dti_feature_list(drug_list, protein_list)

    y_dti_data = DTI_data_preparation.get_DTIs(drug_list=drug_list, protein_list=protein_list)
    y_dti_data = y_dti_data.reshape((len(protein_list), len(drug_list)))
    print("Finished.\n")
    print("Finished loading data.\n")

    # building stellar graph
    print("Building Stellar graph data ...")
    PPI_graph = DTI_data_preparation.get_annotated_PPI_graph()
    G = sg.StellarGraph(PPI_graph, node_features='node_feature')
    print(G.info())

    embedding_batch_size = 50
    protein_node_embeddings = None
    generator = None
    if supervised:
        embedding_layer_sizes[-1] = embedding_output_size
        if embedding_method=='graphsage':
            generator = GraphSAGENodeGenerator(G, embedding_batch_size, num_samples)
        elif embedding_method == 'hinsage':
            generator = HinSAGENodeGenerator(G, embedding_batch_size, num_samples)
        elif embedding_method in ['gcn', 'chebyshev', 'gat']:
            generator = FullBatchNodeGenerator(G, method=embedding_method)
        elif embedding_method == 'sgc':
            generator = FullBatchNodeGenerator(G, method='sgc', k=2)
        elif embedding_method == 'ppnp':
            generator = FullBatchNodeGenerator(G,
                                               method="ppnp",
                                               sparse=False,
                                               teleport_probability=0.1
                                               )
        elif embedding_method == 'appnp':
            generator = FullBatchNodeGenerator(G, method="gcn", sparse=True)
        else:
            print("No valid embedding method chosen.")
            raise Exception
    else:
        if embedding_method == 'attr2vec':
            number_of_walks = 30
            length = 2
            unsupervised_samples = UnsupervisedSampler(G, nodes=list(G.nodes()), length=length, number_of_walks=number_of_walks)
            generator = Attri2VecLinkGenerator(G, embedding_batch_size).flow(unsupervised_samples)
        else:
            print("No valid embedding method chosen.")
            raise Exception

        train_gen = generator

        embedding_layer=None
        if embedding_method == 'attr2vec':
            embedding_layer = Attri2Vec(layer_sizes=embedding_layer_sizes,
                                        generator=train_gen,
                                        bias=True,
                                        normalize=None)
        x_inp, x_out = embedding_layer.build()
        prediction = link_classification(
            output_dim=1, output_act="sigmoid", edge_embedding_method='ip'
        )(x_out)

        model = models.Model(inputs=x_inp, outputs=prediction)

        model.compile(
            optimizer=optimizers.Adam(lr=1e-3),
            loss=losses.binary_crossentropy,
            metrics=['acc'],
        )

        history = model.fit_generator(
            generator,
            epochs=1,
            verbose=1,
            use_multiprocessing=True,
            workers=10000,
            shuffle=True
        )

        x_inp_src = x_inp[0]
        x_out_src = x_out[0]
        embedding_model = models.Model(inputs=x_inp_src, outputs=x_out_src)

        node_gen = Attri2VecNodeGenerator(G, batch_size).flow(protein_list)
        protein_node_embeddings = embedding_model.predict_generator(node_gen, workers=4, verbose=1)

    print("Finished.\n") # Fold parameters
    skf = KFold(n_splits=5, random_state=42)

    cv_scores = {'acc': [],
                 'auroc': [],
                 'f1-score': []}

    model = None
    embedding_model = None
    conf_matrix = None
    round = 1
    print("Starting folds ...")
    for train, test in skf.split(protein_list):
        print("Round", round)
        round += 1

        embedding_layer = None
        if supervised:
            train_gen = generator.flow(protein_list[train], PPI_dti_features[train])

            if embedding_method == 'graphsage':
                embedding_layer = GraphSAGE(
                    layer_sizes=embedding_layer_sizes,
                    generator=train_gen,
                    bias=True,
                    dropout=0.5
                )
            elif embedding_method == 'gcn':
                embedding_layer = GCN(
                    layer_sizes=embedding_layer_sizes,
                    activations=['relu', 'relu'],
                    generator=generator,
                    dropout=0.5
                )
            elif embedding_method == 'gat':
                embedding_layer = GAT(
                    layer_sizes=embedding_layer_sizes,
                    activations=["elu", "relu"],
                    attn_heads=8,
                    generator=generator,
                    in_dropout=0.5,
                    attn_dropout=0.5,
                    normalize=None
                )
            elif embedding_method == 'hinsage':
                embedding_layer = HinSAGE(
                    layer_sizes=embedding_layer_sizes,
                    generator=train_gen,
                    bias=True,
                    dropout=0.5
                )
            elif embedding_method == 'ppnp':
                embedding_layer = PPNP(
                    layer_sizes=embedding_layer_sizes,
                    activations=['relu']*len(embedding_layer_sizes),
                    generator=generator,
                    dropout=0.5,
                    kernel_regularizer=regularizers.l2(0.001)
                )
            elif embedding_method == 'appnp':
                embedding_layer = APPNP(
                    layer_sizes=embedding_layer_sizes,
                    activations=['relu'] * len(embedding_layer_sizes),
                    bias=True,
                    generator=generator,
                    teleport_probability=0.1,
                    dropout=0.5,
                    kernel_regularizer=regularizers.l2(0.001)
                )
            elif embedding_method == 'sgc':
                embedding_layer = GCN(
                    layer_sizes=[embedding_output_size],
                    generator=generator,
                    bias=True,
                    dropout=0.5,
                    activations=["softmax"],
                    kernel_regularizer=regularizers.l2(5e-4),
                )

            x_inp, x_out = embedding_layer.build() if embedding_method not in ['gcn', 'gat', 'ppnp', 'appnp', 'sgc'] else embedding_layer.node_model()

            prediction = layers.Dense(units=PPI_dti_features.shape[1], activation="linear")(x_out)

            encoder = models.Model(inputs=x_inp, outputs=x_out)
            embedding_model = models.Model(inputs=x_inp, outputs=prediction)

            embedding_model.compile(optimizer=optimizers.Adam(lr=0.005),
                                    loss=losses.binary_crossentropy,
                                    metrics=["acc"])

            embedding_model.summary()

            val_gen = generator.flow(protein_list[test], PPI_dti_features[test])

            history = embedding_model.fit_generator(
                train_gen,
                epochs=2,
                # validation_data=val_gen,
                verbose=1,
                shuffle=True
            )

            print("\nCalculate node embeddings ...")
            overall_generator = generator.flow(protein_list)

            protein_node_embeddings = encoder.predict_generator(overall_generator)
            print("Finished.\n")



        print(protein_node_embeddings.shape)

        # if plot:
            # dti_utils.plot_history(history)
        if embedding_method in ['gcn', 'gat', 'ppnp', 'appnp', 'sgc']:
            protein_node_embeddings = protein_node_embeddings.reshape(protein_node_embeddings.shape[1],
                                                                      protein_node_embeddings.shape[2])

        print("Scaling node embeddings ...")
        train_protein_node_embeddings = protein_node_embeddings[train]
        test_protein_node_embeddings = protein_node_embeddings[test]

        train_protein_node_embeddings = np.repeat(train_protein_node_embeddings, len(drug_list), axis=0)
        test_protein_node_embeddings = np.repeat(test_protein_node_embeddings, len(drug_list), axis=0)
        print("Finished.\n")

        print(protein_node_embeddings.shape)

        y_dti_train_data = y_dti_data[train].flatten()
        y_dti_test_data = y_dti_data[test].flatten()

        # Build actual dti model
        PPI_input = layers.Input(shape=(embedding_output_size,))

        # DDI_input = layers.Input(shape=(DDI_features.shape[2],))

        side_effect_input = layers.Input(shape=(side_effect_features.shape[2],))

        merge_1 = layers.Concatenate(axis=1)([PPI_input,
                                              # DDI_input,
                                              side_effect_input
                                              ])

        dense_1 = layers.Dense(500, activation='relu')(merge_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        dense_1 = layers.Dense(100, activation='relu')(dropout_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        dense_1 = layers.Dense(50, activation='relu')(dropout_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        dense_1 = layers.Dense(20, activation='relu')(dropout_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        output = layers.Dense(1, activation='sigmoid')(dropout_1)

        model = models.Model(inputs=[PPI_input,
                                     # DDI_input,
                                     side_effect_input
                                     ],
                             outputs=output)

        model.compile(loss=losses.binary_crossentropy,
                      optimizer=optimizers.Adam(lr=0.005),
                      metrics=[dti_utils.dti_f1_score,
                               dti_utils.dti_auroc,
                               'accuracy'])

        print(y_dti_train_data.sum())

        model.summary()

        imb_ratio = (len(y_dti_train_data)-y_dti_train_data.sum())/y_dti_train_data.sum()

        class_weight = {0: 1.,
                        1: imb_ratio}
        model.fit([train_protein_node_embeddings,
                   # DDI_features[train].reshape((len(drug_list)*len(train), len(drug_list))),
                   side_effect_features[train].reshape((len(drug_list)*len(train), side_effect_features.shape[-1]))
                   ],
                  y_dti_train_data,
                  batch_size=batch_size,
                  validation_data=([test_protein_node_embeddings,
                                   # DDI_features[test].reshape((len(drug_list)*len(test), len(drug_list))),
                                   side_effect_features[test].reshape((len(drug_list)*len(test), side_effect_features.shape[-1]))
                                    ],
                                   y_dti_test_data),
                  epochs=nb_epochs,
                  class_weight=class_weight,
                  shuffle=True,
                  )

        '''
        callbacks = [dti_utils.roc_callback(training_data=([train_protein_node_embeddings,
                                                            DDI_features[train].reshape(
                                                                (len(drug_list) * len(train), len(drug_list)))],
                                                           y_dti_train_data),
                                            validation_data=([test_protein_node_embeddings,
                                                              DDI_features[test].reshape(
                                                                  (len(drug_list) * len(test), len(drug_list)))],
                                                             y_dti_test_data))]
        '''

        y_pred = model.predict([test_protein_node_embeddings,
                                # DDI_features[test].reshape((len(drug_list)*len(test), len(drug_list))),
                                side_effect_features[test].reshape((len(drug_list)*len(test), side_effect_features.shape[-1]))
                                ])

        conf_matrix = metrics.confusion_matrix(y_true=y_dti_test_data,
                                               y_pred=(y_pred.reshape((len(y_pred))) >= 0.5).astype(int))

        print("Confusion matrix", conf_matrix)

        tn = conf_matrix[0, 0]
        tp = conf_matrix[1, 1]
        fp = conf_matrix[0, 1]
        fn = conf_matrix[1, 0]

        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        accuracy = (tp + tn) / (tp + tn + fp + fn)

        y_pred = (y_pred.reshape((y_pred.shape[0]))>0.5).astype(int)

        auroc = metrics.roc_auc_score(y_dti_test_data, y_pred)
        f1_score = metrics.f1_score(y_dti_test_data, y_pred)

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
        print("1" +"\t"+ "0" +"\t"+ "1" +"\t"+ "0" +"\t"+
              str(len(protein_list)) +"\t"+ str(len(drug_list)) +"\t"+
              str(nb_epochs) +"\t"+
              str(embedding_layer_sizes) + "\t" +
              str(supervised) + '\t' +
              str(embedding_method) + "\t" +
              str(np.mean(cv_scores['acc'])) +"\t"+ str(np.mean(cv_scores['auroc'])) +"\t"+ str(np.mean(cv_scores['f1-score'])),
              file=f)
    with open(file=results_filename, mode='a') as filehandler:
        print("DTI "+embedding_method+ " PREDICTION", file=filehandler)

        print("Including:", file=filehandler)
        print("- PPIs", file=filehandler)
        # print("- DDIs", file=filehandler)
        print("- side similarity scores", file=filehandler)
        # print("- HMM node features", file=filehandler)
        print("Number of targets:\t", len(protein_list), file=filehandler)
        print("Number of drugs:\t", len(drug_list), file=filehandler)

        print("Mean accuracy: {},\t Std: {}".format(np.mean(cv_scores['acc']), np.std(cv_scores['acc'])), file=filehandler)
        print("Mean auroc: {},\t Std: {}".format(np.mean(cv_scores['auroc']), np.std(cv_scores['auroc'])), file=filehandler)
        print("Mean f1: {},\t Std: {}".format(np.mean(cv_scores['f1-score']), np.std(cv_scores['f1-score'])), file=filehandler)

        print("Epochs: {}, Batch size: {}".format(nb_epochs, batch_size), file=filehandler)
        model.summary(print_fn=lambda x: filehandler.write(x + '\n'))

        # print confusion matrix
        print("Confusion matrix:",
              conf_matrix,
              file=filehandler)

        print("------------------------------------------------\n")

        # serialize model to JSON
        if plot:
            tf.keras.utils.plot_model(embedding_model,
                                      show_shapes=True,
                                      to_file='../models/'+embedding_method+'_model.png')
            tf.keras.utils.plot_model(model,
                                      show_shapes=True,
                                      to_file='../models/overall_model.png')



def GCN_missing_target_predictor():
    pass
    '''
    # build graph conv filters
    SYM_NORM = True
    num_filters = 32
    graph_conv_filters = preprocess_adj_tensor_with_identity(PPI_adj_mats[train], SYM_NORM)

    X_input = Input(shape=(PPI_node_features.shape[1], PPI_node_features.shape[2]))
    # A_input = Input(shape=(PPI_adj_mats.shape[1], PPI_adj_mats.shape[2]))
    graph_conv_filters_input = Input(shape=(graph_conv_filters.shape[1], graph_conv_filters.shape[2]))

    convolutional_1 = MultiGraphCNN(32, num_filters, graph_conv_filters, activation='elu')([X_input, graph_conv_filters_input])
    dropout_1 = layers.Dropout(0.4)(convolutional_1)
    convolutional_2 = MultiGraphCNN(64, num_filters, graph_conv_filters, activation='elu')([dropout_1, graph_conv_filters_input])
    dropout_2 = layers.Dropout(0.4)(convolutional_2)
    convolutional_3 = MultiGraphCNN(128, num_filters, graph_conv_filters, activation='relu')([dropout_2, graph_conv_filters_input])
    dropout_3 = layers.Dropout(0.2)(convolutional_3)

    graph_output = layers.Lambda(lambda x: K.mean(x, axis=1))(dropout_3)

    flatten_graph = layers.Flatten()(graph_output)

    DDI_input = Input(shape=(DDI_features.shape[1],))
    flatten_DDIs = layers.Flatten()(DDI_input)

    merge_1 = layers.Concatenate(axis=1)([flatten_graph, flatten_DDIs])

    dense_1 = layers.Dense(100, activation='relu')(merge_1)
    dropout_1 = layers.Dropout(0.5)(dense_1)
    dense_2 = layers.Dense(10, activation='relu')(dropout_1)
    dropout_2 = layers.Dropout(0.5)(dense_2)
    output = layers.Dense(1, activation='sigmoid')(dropout_2)


    model = Model(inputs=[X_input, graph_conv_filters_input, DDI_input], outputs=output)
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=[dti_auroc, 'acc'])

    model.fit([PPI_node_features[train],
               graph_conv_filters[train],
               DDI_features[train]],
              y_train_dti_data,
              batch_size=batch_size,
              validation_data=y_test_dti_data,
              verbose=1)
    '''

def pure_HMM_predictor():
    drug_to_Hmm_filtered_targets_dict = DTI_data_preparation.get_drug_to_HMM_filtered_targets_dict()
    drug_list = DTI_data_preparation.get_drug_list()
    protein_list = DTI_data_preparation.get_human_proteins()

    y_dtis = DTI_data_preparation.get_DTIs(drug_list, protein_list)

    y_pred = np.zeros(len(protein_list) * len(drug_list))

    print(len(set(drug_to_Hmm_filtered_targets_dict.keys()) & set(drug_list)))

    raise Exception

    for i in range(len(drug_list)):
        drug = drug_list[i]
        for j in range(len(protein_list)):
            protein = protein_list[j]

            if protein in drug_to_Hmm_filtered_targets_dict[drug]:
                y_pred[j * len(drug_list) + i] = 1

    conf_matrix = metrics.confusion_matrix(y_true=y_dtis, y_pred=y_pred)

    tn = conf_matrix[0, 0]
    tp = conf_matrix[1, 1]
    fp = conf_matrix[0, 1]
    fn = conf_matrix[1, 0]

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    accuracy = (tp + tn) / (tp + tn + fp + fn)

    auroc = metrics.roc_auc_score(y_dtis, y_pred)
    f1_score = metrics.f1_score(y_dtis, y_pred)

    print("accuracy", accuracy * 100)
    print("prec", precision * 100)
    print("recall", recall * 100)
    print("auroc", auroc * 100)
    print("f1-score", f1_score * 100)



if __name__ == '__main__':
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, num_samples=[50, 50], embedding_layer_sizes= [32, 32], embedding_output_size=32, embedding_method='attr2vec', supervised=False)
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, num_samples=[50, 50], embedding_layer_sizes= [32, 32], embedding_output_size=64)
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, num_samples=[100, 100], embedding_layer_sizes= [128, 64], embedding_output_size=64)
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, num_samples=[200, 100], embedding_layer_sizes= [128, 64], embedding_output_size=128)

    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, embedding_layer_sizes=[32,32], embedding_output_size=64, embedding_method='gcn')
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, embedding_layer_sizes=[32,32], embedding_output_size=64, embedding_method='gat')
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, embedding_layer_sizes=[32], num_samples=[100], embedding_output_size=64, embedding_method='hinsage')
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, embedding_layer_sizes=[64, 64, 64], embedding_output_size=128, embedding_method='ppnp')
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, embedding_layer_sizes=[64, 64, 64], embedding_output_size=128, embedding_method='appnp')
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, embedding_layer_sizes=[64, 64, 64], embedding_output_size=128, embedding_method='sgc')
    # missing_target_predictor(batch_size=10000, nb_epochs=20, plot=True, embedding_layer_sizes=[64, 64], embedding_output_size=128, embedding_method='graphsage')

    pure_HMM_predictor()

