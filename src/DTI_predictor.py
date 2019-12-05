import numpy as np
import networkx as nx
import pandas as pd

import DTI_data_preparation

from sklearn.model_selection import KFold
from sklearn import metrics

import matplotlib.pyplot as plt

import tensorflow as tf

import tensorflow.keras.models as models
import tensorflow.keras.layers as layers
import tensorflow.keras.backend as K
import tensorflow.keras.regularizers as regularizers
import tensorflow.keras.optimizers as optimizers
import tensorflow.keras.losses as losses


from keras_dgl.layers.graph_cnn_layer import GraphCNN
from keras_dgl.layers.graph_attention_cnn_layer import GraphAttentionCNN
from keras_dgl.layers.multi_graph_cnn_layer import MultiGraphCNN
from keras_dgl.layers.multi_graph_attention_cnn_layer import MultiGraphAttentionCNN
from keras_dgl.utils import *

import stellargraph as sg
from stellargraph import globalvar
from stellargraph.mapper import GraphSAGENodeGenerator
from stellargraph.layer import GraphSAGE

import dti_utils



def missing_target_predictor(results_filename='../results/results_log',
                             nb_epochs=3,
                             batch_size=500,
                             plot=False):

    print("Loading data ...")
    drug_list = np.array(DTI_data_preparation.get_drug_list())
    print("Get protein list ...")
    protein_list = np.array(DTI_data_preparation.get_human_proteins())[:1000]
    print("Finished.\n")

    # side_effect_features = DTI_data_preparation.get_side_effect_similarity_feature_list()
    print("Scaling data ...")
    DDI_features = np.tile(DTI_data_preparation.get_DDI_feature_list(drug_list), (len(protein_list),1))
    DDI_features = DDI_features.reshape((len(protein_list), len(drug_list), len(drug_list)))
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

    graphsage_batch_size = 50
    num_samples = [50, 20] # What do those values mean? [100,50,20]

    generator = GraphSAGENodeGenerator(G, graphsage_batch_size, num_samples)
    print("Finished.\n") # Fold parameters
    skf = KFold(n_splits=5, random_state=42)

    cv_scores = {'acc': [],
                 'auroc': [],
                 'f1-score': []}

    # model = None
    conf_matrix = None
    round = 0
    print("Starting folds ...")
    for train, test in skf.split(protein_list):
        print("Round", round)
        round += 1

        # parameters
        graphsage_output_size = 64

        train_gen = generator.flow(protein_list[train], PPI_dti_features[train], shuffle=True)

        graphsage_model_layer = GraphSAGE(layer_sizes=[32, graphsage_output_size],
                                    generator=generator,
                                    bias=True,
                                    dropout=0.5)

        x_inp, x_out = graphsage_model_layer.build()

        prediction = layers.Dense(units=PPI_dti_features.shape[1], activation="softmax")(x_out)

        encoder = models.Model(inputs=x_inp, outputs=x_out)
        graphsage_model = models.Model(inputs=x_inp, outputs=prediction)

        graphsage_model.compile(optimizer=optimizers.Adam(lr=0.005),
                                loss=losses.binary_crossentropy,
                                metrics=["acc"])

        graphsage_model.summary()

        # val_gen = generator.flow(protein_list[test], protein_list[test], shuffle=True)

        history = graphsage_model.fit_generator(
            train_gen,
            epochs=2,
            # validation_data=val_gen,
            verbose=1,
            shuffle=False
        )

        if plot:
            dti_utils.plot_history(history)

        overall_generator = generator.flow(protein_list)

        protein_node_embeddings = encoder.predict_generator(overall_generator)

        train_protein_node_embeddings = protein_node_embeddings[train]
        test_protein_node_embeddings = protein_node_embeddings[test]

        train_protein_node_embeddings = np.repeat(train_protein_node_embeddings, len(drug_list), axis=0)
        test_protein_node_embeddings = np.repeat(test_protein_node_embeddings, len(drug_list), axis=0)

        print(protein_node_embeddings.shape)

        y_dti_train_data = y_dti_data[train].flatten()
        y_dti_test_data = y_dti_data[test].flatten()

        # Build actual dti model
        PPI_input = layers.Input(shape=(graphsage_output_size,))

        DDI_input = layers.Input(shape=(DDI_features.shape[1],))

        merge_1 = layers.Concatenate(axis=1)([PPI_input, DDI_input])

        dense_1 = layers.Dense(100, activation='relu')(merge_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        dense_1 = layers.Dense(100, activation='relu')(dropout_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        dense_1 = layers.Dense(100, activation='relu')(dropout_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        dense_1 = layers.Dense(100, activation='relu')(dropout_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        output = layers.Dense(1, activation='sigmoid')(dropout_1)

        model = models.Model(inputs=[PPI_input, DDI_input],
                             outputs=output)

        model.compile(loss=losses.binary_crossentropy,
                      optimizers=optimizers.Adam,
                      metrics=[dti_utils.dti_auroc,
                               'accuracy'])

        model.summary()

        model.fit([train_protein_node_embeddings,
                   DDI_features[train].reshape((len(drug_list)*len(train), len(drug_list)))],
                  y_dti_train_data,
                  batch_size=batch_size,
                  validation_data=([test_protein_node_embeddings,
                                   DDI_features[test].reshape((len(drug_list)*len(test), len(drug_list)))],
                                   y_dti_test_data),
                  epochs=nb_epochs,
                  shuffle=True
                  )

        y_pred = model.predict([test_protein_node_embeddings,
                                DDI_features[test].reshape((len(drug_list)*len(test), len(drug_list)))])

        conf_matrix = metrics.confusion_matrix(y_true=y_dti_test_data,
                                               y_pred=(y_pred.reshape((len(y_pred))) >= 0.5).astype(int))

        print("Confusion matrix", conf_matrix)



        raise Exception

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

if __name__ == '__main__':
    missing_target_predictor()
