import pickle
import numpy as np

def write_predicted_DTIs(fold=3):
    filename = 'PPI_network_model_with_mol_features_fold_'+str(fold)+'_predictions.pkl'

    pred_list = None #[(drug, prot, label, pred)]
    with open(file=filename, mode='rb') as f:
        pred_list = pickle.load(f)

    pred_list = [tup for tup in pred_list if tup[label]]

    pred_list = sorted(pred_list, key=lambda tup: tup)