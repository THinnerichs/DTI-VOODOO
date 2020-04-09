import numpy as np

import torch

from smiles_transformer.pretrain_trfm import TrfmSeq2seq
from smiles_transformer.pretrain_rnn import RNNSeq2Seq
from smiles_transformer.build_vocab import WordVocab
from smiles_transformer.utils import split


def encode_drugs(druglist,
                 mode='trfm'):
    pretrained_model_dir = '../models/drug_representation/'
    vocab = WordVocab.load_vocab(pretrained_model_dir + 'vocab.pkl')

    if mode=='trfm':
        trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 3)
        trfm.load_state_dict(torch.load(pretrained_model_dir+'trfm.pkl'))
        trfm.eval()
