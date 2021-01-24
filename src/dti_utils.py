import numpy as np
# import tensorflow as tf
import sklearn.metrics as metrics
# import tensorflow.keras.backend as K
import matplotlib.pyplot as plt

import math


def dti_auroc(y_true, y_pred):
    if y_true.sum() == 0 or (1-y_true).sum() == 0:
        return metrics.accuracy_score(y_true=y_true, y_pred=y_pred)
    return metrics.roc_auc_score(y_true, y_pred)

def dti_auprc(y_true, y_pred):
    p,r, t = metrics.precision_recall_curve(y_true, y_pred)
    return metrics.auc(r, p)
    # return metrics.average_precision_score(y_true, y_pred, average='weighted')

def dti_mcc(y_true, y_pred):
    return metrics.matthews_corrcoef(y_true, y_pred)

def dti_f1_score(y_true, y_pred):
    return metrics.f1_score(y_true, y_pred)

def micro_AUC_per_prot(y_true, y_pred, num_drugs):
    y_true = y_true.reshape(num_drugs, -1)
    y_pred = y_pred.reshape(num_drugs, -1)

    num_prots = y_true.shape[1]
    prot_wise_auc = []
    for prot_index in range(num_prots):
        if y_true[:, prot_index].sum() == 0:
            pass
        prot_wise_auc.append(dti_auroc(y_true[:, prot_index], y_pred[:, prot_index]))

    prot_wise_auc = np.array(prot_wise_auc)

    return prot_wise_auc.mean()

def micro_AUC_per_drug(y_true, y_pred, num_drugs):
    y_true = y_true.reshape(num_drugs, -1)
    y_pred = y_pred.reshape(num_drugs, -1)

    num_prots = y_true.shape[1]
    drug_wise_auc = []
    for drug_index in range(num_drugs):
        if y_true[drug_index, :].sum() == 0:
            pass
        drug_wise_auc.append(dti_auroc(y_true[drug_index, :], y_pred[drug_index, :]))

    drug_wise_auc = np.array(drug_wise_auc)

    return drug_wise_auc.mean()

