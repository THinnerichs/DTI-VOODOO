import numpy as np
import tensorflow as tf
import sklearn.metrics as metrics
# import tensorflow.keras.backend as K
import matplotlib.pyplot as plt

import math


def dti_auroc(y_true, y_pred):
    if len(np.unique(y_true)) == 1:  # bug in roc_auc_score
        return 0.5
    return metrics.roc_auc_score(y_true, y_pred)

def dti_f1_score(y_true, y_pred):
    if len(np.unique(y_true)) == 1:  # bug in roc_auc_score
        return 0.
    return metrics.f1_score(y_true, y_pred)

'''
def tf_dti_auroc(y_true, y_pred):
    return tf.py_func(dti_auroc_fix, (y_true, y_pred), tf.double)

def tf_dti_f1_score(y_true, y_pred):
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


'''
