import numpy as np
import gensim

import torch

import pickle

import DTI_data_preparation


def build_protein_embeddings(include_ppi=False):
    protein_list = DTI_data_preparation.get_human_proteins()

    DL2vec_path_prefix = '../data/DL2vec/DL2vec_embeddings/'
    uberon_model_filename = GO_model_filename = phenotype_model_filename = None

    if include_ppi:
        uberon_model_filename = DL2vec_path_prefix + 'uberon_intersection_ppi_embedding'
        GO_model_filename = DL2vec_path_prefix + 'go_intersection_ppi_embedding'
        phenotype_model_filename = DL2vec_path_prefix + 'mp_intersection_ppi_embedding'
    else:
        uberon_model_filename = DL2vec_path_prefix + 'uberon_intersection_embedding'
        GO_model_filename = DL2vec_path_prefix + 'go_intersection_embedding'
        phenotype_model_filename = DL2vec_path_prefix + 'mp_intersection_embedding'

    # load models
    uberon_model = gensim.models.Word2Vec.load(uberon_model_filename)
    GO_model = gensim.models.Word2Vec.load(GO_model_filename)
    phenotype_model = gensim.models.Word2Vec.load(phenotype_model_filename)

    # Build wordvector dicts
    uberon_model = uberon_model.wv
    GO_model = GO_model.wv
    phenotype_model = phenotype_model.wv

    uberon_embeddings = []
    GO_embeddings = []
    phenotype_embeddings = []
    for protein in protein_list:
        organism, protein_id = protein.strip().split('.')

        uberon_embeddings.append(uberon_model[protein_id])
        GO_embeddings.append(GO_model[protein_id])
        phenotype_embeddings.append(phenotype_model[protein_id])

    return np.array(uberon_embeddings), np.array(GO_embeddings), np.array(phenotype_embeddings)


if __name__ == '__main__':
    build_protein_embeddings()


