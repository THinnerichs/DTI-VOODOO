import numpy as np
import networkx as nx
from tqdm import tqdm

from sklearn.model_selection import KFold
from sklearn import metrics

from torch_dti_utils import *
from torch_networks import *

import torch
import torch.nn.functional as F
from torch.nn import Linear
import torch_geometric.data as data


def enlightened_missing_target_predictor(results_filename='../results/torched_results_log',
                                         num_epochs=3,
                                         plot=False,
                                         embedding_layer_sizes=[32, 64],
                                         embedding_method='gcn'):
    # Learning parameters
    batch_size = 32
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    print("Loading data ...")
    dataset = FullNetworkDataset("../data/torch_raw/")
    # dataset is present in dimension (num_drugs * num_proteins)
    print("Finished.")

    # generate indices for proteins
    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    X = np.zeros((dataset.num_proteins,1))

    # build for help matrix for indices
    help_matrix = np.arange(dataset.num_drugs * dataset.num_proteins)
    help_matrix = help_matrix.reshape((dataset.num_drugs, dataset.num_proteins))

    val_losses, accs, durations = [], [], []
    fold = 0
    for train_protein_indices, test_protein_indices in kf.split(X):
        fold += 1
        print("Fold:", fold)


        # build train data over whole dataset with help matrix
        train_indices = help_matrix[:, train_protein_indices].flatten()
        test_indices = help_matrix[:, test_protein_indices].flatten()
        print(train_indices.shape, test_indices.shape)

        train_dataset = dataset[train_indices]
        test_dataset = dataset[test_indices]

        # build DataLoaders
        train_loader = data.DataLoader(train_dataset, batch_size, shuffle=True)
        test_loader = data.DataLoader(test_dataset, batch_size, shuffle=False)

        if torch.cuda.is_available():
            torch.cuda.synchronize(device=device)

        model = SimpleConvGCN(num_drugs=dataset.num_drugs,
                              num_features=dataset.num_PPI_features)
        optimizer = torch.optim.Adam([dict(params=model.reg_params, weight_decay=5e-4),
                                      dict(params=model.non_reg_params, weight_decay=0)],
                                     lr=0.01)

        for epoch in range(1, num_epochs + 1):
            train_loss = train(model=model, device=device, train_loader=train_loader, optimizer=optimizer, epoch=epoch)
            val_loss = eval_loss(model, test_loader, device)
            val_losses.append(val_loss)
            test_acc = eval_acc(model, test_loader, device)
            accs.append(test_acc)

            print('{:02d}/{:03d}: Val Loss: {:.4f}, Test Accuracy: {:.3f}'.format(fold, epoch, val_loss, test_acc))

        if torch.cuda.is_available():
            torch.cuda.synchronize(device=device)

    loss, acc = np.array(val_losses), np.array(accs)
    print("Acc mean:", acc.mean())

    print("Done.")













if __name__ == '__main__':
    enlightened_missing_target_predictor()