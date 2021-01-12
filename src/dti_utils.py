import numpy as np
# import tensorflow as tf
import sklearn.metrics as metrics
# import tensorflow.keras.backend as K
import matplotlib.pyplot as plt

import math


def dti_auroc(y_true, y_pred):
    if len(np.unique(y_true)) == 1:  # bug in roc_auc_score
        return 0.5
    if len(np.unique(y_pred)) == 0:
        print('y_pred.sum()', np.unique(y_pred), y_pred.sum())
    if y_pred.sum() == 0 or (1-y_pred).sum() == 0 or y_true.sum() == 0 or (1-y_true).sum()==0:
        return 0.5
    return metrics.roc_auc_score(y_true, y_pred)

def dti_auprc(y_true, y_pred):
    if len(np.unique(y_true)) == 1 or len(np.unique(y_pred)==1):  # bug in sklearn
        return 0.5
    if y_pred.sum() == 0 or (1-y_pred).sum() == 0 or y_true.sum() == 0 or (1-y_true).sum()==0:
        return 0.5
    return metrics.average_precision_score(y_true, y_pred)

def dti_mcc(y_true, y_pred):
    if len(np.unique(y_true)) == 1 or len(np.unique(y_pred)==1):  # bug in sklearn
        return 0.
    return metrics.matthews_corrcoef(y_true, y_pred)

def dti_f1_score(y_true, y_pred):
    if len(np.unique(y_true)) == 1 or len(np.unique(y_pred)==1):  # bug in sklearn
        return 0.
    return metrics.f1_score(y_true, y_pred)
