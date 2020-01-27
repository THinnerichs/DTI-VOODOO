import numpy as np
import networkx as nx
from tqdm import tqdm

from sklearn.model_selection import KFold
from sklearn import metrics

from torch_dti_utils import *

import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv

class SimpleConvGCN(torch.nn.Module):
    def __init__(self, num_features):
        super(SimpleConvGCN, self).__init__()
        self.conv1 = GCNConv(num_features, 16, cached=False, normalize=True)


def enlightened_missing_target_predictor(results_filename='../results/torched_results_log',
                                         nb_epochs=3,
                                         plot=False,
                                         embedding_layer_sizes=[32, 64],
                                         embedding_method='gcn'):
    print("Loading data ...")
    dataset = FullNetworkDataset("../data/torch_raw/")
    # dataset is present in dimension (num_drugs * num_proteins)
    print("Finished.")

    # generate indices for proteins
    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    X = np.zeros((dataset.num_proteins,1))

    for train_index, test_index in kf.split(X):
        pass





if __name__ == '__main__':
    enlightened_missing_target_predictor()