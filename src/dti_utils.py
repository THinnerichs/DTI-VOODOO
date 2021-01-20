import numpy as np
# import tensorflow as tf
import sklearn.metrics as metrics
# import tensorflow.keras.backend as K
import matplotlib.pyplot as plt

import math


def dti_auroc(y_true, y_pred):
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
        prot_wise_auc.append(metrics.roc_auc_score(y_true[:, prot_index], y_pred[:, prot_index]))




