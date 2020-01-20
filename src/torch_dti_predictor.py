import numpy as np
import networkx as nx
from tqdm import tqdm

from sklearn.model_selection import KFold
from sklearn import metrics

from torch_dti_utils import *

import torch
import torch.nn.functional as F

def enlightened_missing_target_predictor(results_filename='../results/torched_reults_log',
                                         nb_epochs=3,
                                         plot=False,
                                         embedding_layer_sizes=[32, 64],
                                         embedding_method='gcn'):
    print("Loading data ...")
    dataset = FullNetworkDataset("../data/torch_raw/")

    print("Finished.")

