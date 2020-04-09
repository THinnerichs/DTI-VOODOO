import numpy as np

import torch

from smiles_transformer.pretrain_trfm import TrfmSeq2seq
# from smiles_transformer.pretrain_rnn import RNNSeq2Seq
from smiles_transformer.build_vocab import WordVocab
from smiles_transformer.utils import split

import DTI_data_preparation

import DDI_utils


def encode_drugs(drug_list,
                 mode='trfm'):

    drug_SMILES_dict = DTI_data_preparation.get_truncated_drug_to_SMILES_dict()

    pad_index = 0
    unk_index = 1
    eos_index = 2
    sos_index = 3
    mask_index = 4

    # Adapted from https://github.com/DSPsleeporg/smiles-transformer
    pretrained_model_dir = '../models/drug_representation/'
    vocab = WordVocab.load_vocab(pretrained_model_dir + 'vocab.pkl')

    def get_inputs(sm):
        seq_len = 220
        sm = sm.split()
        if len(sm) > 218:
            print('SMILES is too long ({:d})'.format(len(sm)))
            sm = sm[:109] + sm[-109:]
        ids = [vocab.stoi.get(token, unk_index) for token in sm]
        ids = [sos_index] + ids + [eos_index]
        seg = [1] * len(ids)
        padding = [pad_index] * (seq_len - len(ids))
        ids.extend(padding), seg.extend(padding)
        return ids, seg

    def get_array(smiles):
        x_id, x_seg = [], []
        for sm in smiles:
            a, b = get_inputs(sm)
            x_id.append(a)
            x_seg.append(b)
        return torch.tensor(x_id), torch.tensor(x_seg)

    if mode=='trfm':
        trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
        trfm.load_state_dict(torch.load(pretrained_model_dir+'trfm.pkl'))
        trfm.eval()

        x_split = [split(sm) for sm in [drug_SMILES_dict[drug] for drug in drug_list]]
        xid, xseg = get_array(x_split)
        X = trfm.encode(torch.t(xid))

        print('Size:', X.size())

if __name__=='__main__':
    drug_list = DTI_data_preparation.get_drug_list()

    encode_drugs(drug_list)