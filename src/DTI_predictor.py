import numpy as np
import networkx as nx
import pandas as pd

import DTI_data_preparation

from sklearn.model_selection import KFold
from sklearn import metrics

from keras.models import Input, Model
from keras.utils import plot_model
from keras import layers
import keras.backend as K
from keras.regularizers import *

from keras_dgl.layers.graph_cnn_layer import GraphCNN
from keras_dgl.layers.graph_attention_cnn_layer import GraphAttentionCNN
from keras_dgl.utils import *

import stellargraph as sg
from stellargraph import globalvar
from stellargraph.mapper import GraphSAGENodeGenerator
from stellargraph.layer import GraphSAGE

import PPI_utils


def missing_drug_predictor(results_filename='../results/results_log',
                           nb_epochs=3,
                           batch_size=500,
                           plot=False):

    drug_list = DTI_data_preparation.get_drug_list()
    protein_list = DTI_data_preparation.get_human_proteins()

    # side_effect_features = DTI_data_preparation.get_side_effect_similarity_feature_list()
    # DDI_features = DTI_data_preparation.get_DDI_feature_list()
    # PPI_adj_mats = DTI_data_preparation.get_PPI_adj_mat_list()
    PPI_node_features = DTI_data_preparation.get_PPI_node_feature_mat_list()

    # skf = KFold(n_splits=5, random_state=42)

    cv_scores = {'acc': [],
                 'auroc': [],
                 'f1-score': []}

    df_node_features = pd.DataFrame(PPI_node_features, index=protein_list)
    PPI_graph = PPI_utils.get_PPI_graph()
    G = sg.StellarGraph(PPI_graph, node_features=df_node_features)

    print(G.info())

    batch_size = 50
    num_samples = [10, 10] # What do those values mean?

    generator = GraphSAGENodeGenerator(G, batch_size, num_samples)

    print("Generator", next(generator))








if __name__ == '__main__':
    missing_drug_predictor()
