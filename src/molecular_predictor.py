import numpy as np

import torch

from smiles_transformer.pretrain_trfm import TrfmSeq2seq
from smiles_transformer.pretrain_rnn import RNNSeq2Seq
from smiles_transformer.build_vocab import WordVocab
from smiles_transformer.utils import split

import DTI_data_preparation

import DDI_utils


def encode_drugs(drug_list,
                 mode='trfm'):
    drug_SMILES_dict = DTI_data_preparation.get_truncated_drug_to_SMILES_dict()

    pretrained_model_dir = '../models/drug_representation/'
    vocab = WordVocab.load_vocab(pretrained_model_dir + 'vocab.pkl')

    if mode=='trfm':
        trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 3)
        trfm.load_state_dict(torch.load(pretrained_model_dir+'trfm.pkl'))
        trfm.eval()

if __name__=='__main__':
    drug_list = DTI_data_preparation.get_drug_list()

    encode_drugs(drug_list)