import numpy as np
import networkx as nx
import pandas as pd

import DTI_data_preparation

from sklearn.model_selection import KFold
from sklearn import metrics

import matplotlib.pyplot as plt

from keras.models import Input, Model
from keras.utils import plot_model
from keras import layers
import keras.backend as K
from keras.regularizers import *
from keras import optimizers, losses

import tensorflow as tf

from keras_dgl.layers.graph_cnn_layer import GraphCNN
from keras_dgl.layers.graph_attention_cnn_layer import GraphAttentionCNN
from keras_dgl.layers.multi_graph_cnn_layer import MultiGraphCNN
from keras_dgl.layers.multi_graph_attention_cnn_layer import MultiGraphAttentionCNN
from keras_dgl.utils import *

import stellargraph as sg
from stellargraph import globalvar
from stellargraph.mapper import GraphSAGENodeGenerator
from stellargraph.layer import GraphSAGE

import PPI_utils



def dti_auroc(y_true, y_pred):
    return tf.py_func(metrics.roc_auc_score, (y_true, y_pred), tf.double)

def dti_f1_score(y_true, y_pred):
    def recall(y_true, y_pred):
        """Recall metric.

        Only computes a batch-wise average of recall.

        Computes the recall, a metric for multi-label classification of
        how many relevant items are selected.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
        recall = true_positives / (possible_positives + K.epsilon())
        return recall

    def precision(y_true, y_pred):
        """Precision metric.

        Only computes a batch-wise average of precision.

        Computes the precision, a metric for multi-label classification of
        how many selected items are relevant.
        """
        true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
        predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
        precision = true_positives / (predicted_positives + K.epsilon())
        return precision
    precision = precision(y_true, y_pred)
    recall = recall(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))

def plot_history(history):
    metrics = sorted(history.history.keys())
    metrics = metrics[:len(metrics)//2]
    for m in metrics:
        # summarize history for metric m
        plt.plot(history.history[m])
        plt.plot(history.history['val_' + m])
        plt.title(m)
        plt.ylabel(m)
        plt.xlabel('epoch')
        plt.legend(['train', 'test'], loc='best')
        plt.show()



def missing_drug_predictor(results_filename='../results/results_log',
                           nb_epochs=3,
                           batch_size=500,
                           plot=False):

    print("Loading data ...")
    drug_list = DTI_data_preparation.get_drug_list()
    print("Get protein list ...")
    protein_list = DTI_data_preparation.get_human_proteins()
    print("Finished.\n")

    # side_effect_features = DTI_data_preparation.get_side_effect_similarity_feature_list()
    DDI_features = np.repeat(DTI_data_preparation.get_DDI_feature_list(), len(protein_list))
    # PPI_adj_mats = np.tile(DTI_data_preparation.get_PPI_adj_mat_list(), len(drug_list))
    PPI_node_features = DTI_data_preparation.get_PPI_node_feature_mat_list()
    # PPI_node_features = PPI_node_features.reshape(PPI_node_features.shape + (1,))

    PPI_dti_features = DTI_data_preparation.get_PPI_dti_feature_list()


    print("Finished.\n")

    skf = KFold(n_splits=5, random_state=42)

    cv_scores = {'acc': [],
                 'auroc': [],
                 'f1-score': []}

    df_node_features = pd.DataFrame(PPI_node_features, index=protein_list)
    PPI_graph = PPI_utils.get_PPI_graph()
    G = sg.StellarGraph(PPI_graph, node_features=df_node_features)

    print(G.info())

    graphsage_batch_size = 50
    num_samples = [10, 10] # What do those values mean?

    generator = GraphSAGENodeGenerator(G, graphsage_batch_size, num_samples)

    model = None
    conf_matrix = None
    round = 0
    for train, test in skf.split(protein_list):
        print("Round", round)
        round += 1

        train_gen = generator.flow(protein_list[train], PPI_dti_features[train], shuffle = True)

        graphsage_model = GraphSAGE(layer_sizes=[32, 32],
                                    generator=generator,
                                    bias=True,
                                    dropout=0.5)

        x_inp, x_out = graphsage_model.build()

        prediction = layers.Dense(units=PPI_dti_features.shape[1], activation="softmax")(x_out)

        graphsage_model = Model(inputs=x_inp, outputs=prediction)

        graphsage_model.compile(optimizer=optimizers.Adam(lr=0.005),
                                loss=losses.categorical_crossentropy,
                                metrics=["acc"])

        graphsage_model.summary()

        val_gen = generator.flow(protein_list[test], protein_list[test])

        history = model.fit_generator(
            train_gen,
            epochs=15,
            validation_data=val_gen,
            verbose=0,
            shuffle=False
        )

        if plot:
            plot_history(history)

        # encoder = Model(input=x_inp, output=)

        overall_generator = generator.flow(protein_list, PPI_dti_features)

        node_embeddings = graphsage_model.predict_generator()

        print(node_embeddings.shape)

        raise Exception


        '''
        # build graph conv filters
        SYM_NORM = True
        num_filters = 32
        graph_conv_filters = preprocess_adj_numpy(PPI_adj_mats, SYM_NORM)
        graph_conv_filters = np.concatenate([graph_conv_filters, np.matmul(graph_conv_filters, graph_conv_filters)], axis=0)
        graph_conv_filters = K.constant(graph_conv_filters)

        X_input = Input(shape=(PPI_node_features.shape[1], ))
        # graph_conv_filters_input = Input(shape=(graph_conv_filters.shape[1], graph_conv_filters.shape[2]))

        convolutional_1 = GraphCNN(32, num_filters, graph_conv_filters, activation='elu')(X_input)
        dropout_1 = layers.Dropout(0.4)(convolutional_1)
        convolutional_2 = GraphCNN(64, num_filters, graph_conv_filters, activation='elu')(dropout_1)
        dropout_2 = layers.Dropout(0.4)(convolutional_2)
        convolutional_3 = GraphCNN(128, num_filters, graph_conv_filters, activation='relu')(dropout_2)
        dropout_3 = layers.Dropout(0.2)(convolutional_3)

        flatten_graph = layers.Flatten()(dropout_3)

        DDI_input = Input(shape=(DDI_features.shape[1],))
        flatten_DDIs = layers.Flatten()(DDI_input)

        merge_1 = layers.Concatenate(axis=1)([flatten_graph, flatten_DDIs])

        dense_1 = layers.Dense(100, activation='relu')(merge_1)
        dropout_1 = layers.Dropout(0.5)(dense_1)
        dense_2 = layers.Dense(10, activation='relu')(dropout_1)
        dropout_2 = layers.Dropout(0.5)(dense_2)
        output = layers.Dense(1, activation='sigmoid')(dropout_2)


        model = Model(inputs=[X_input, DDI_input], outputs=output)
        model.compile(loss='binary_crossentropy', optimizer='adam', metrics=[dti_auroc, 'acc'])

        model.fit([PPI_adj_mats, DDI_features])

    '''









if __name__ == '__main__':
    missing_drug_predictor()
