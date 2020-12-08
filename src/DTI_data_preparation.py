import numpy as np
import networkx as nx
import gensim

import torch

import pickle

import os
from tqdm import tqdm
import argparse

import PPI_utils
import DDI_utils
import similarity_measurement
import HPO_GO_similarity
import PhenomeNET_DL2vec_utils




def write_human_DTI_graph(min_score=700,
                          mode=''):
    filename = "../data/STITCH_data/9606.protein_chemical.links.transfer.v5.0.tsv"
    dti_graph = nx.Graph()

    # drug_set = similarity_measurement.get_SIDER_Boyce_Drubank_drug_intersection()
    # drug_set = get_drug_list()

    print("Loading chemical stereo to mono mapping...")
    stereo_mono_mapping = DDI_utils.get_chemical_stereo_to_normal_mapping()
    print("Done.\n")

    print("Parsing human drug-protein-links data ...")
    with open(file=filename, mode='r') as f:
        f.readline()
        for line in tqdm(f, total=15473940):
            split_line = line.split('\t')
            drug = split_line[0].strip()
            if 's' in drug:
                drug = stereo_mono_mapping.get(drug, None)
                if not drug:
                    continue
            # drug = split_line[0].replace('s','m')
            target = split_line[1]

            score = None
            if mode=='experimental':
                score = int((1- (1-int(split_line[2])/1000) * (1-int(split_line[3])/1000))*1000)
            elif mode=='database':
                score = int((1- (1-int(split_line[2])/1000) * (1-int(split_line[3])/1000) * (1-int(split_line[6])/1000) * (1-int(split_line[7])/1000))*1000)
            else:
                score = int(split_line[-1])
            # if not drug in drug_set:
                # continue

            if score >= min_score:
                dti_graph.add_node(drug)
                dti_graph.add_node(target)
                dti_graph.add_edge(drug, target, score=score)

    print("Finished.\n")

    print('num_nodes', len(dti_graph.nodes()))
    print('num_edges', len(dti_graph.edges()))

    print("Writing human only DTI-graph to disk ...")
    filename = "../data/STITCH_data/human_only_"+(mode+'_' if mode else '')+"DTI_graph"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(dti_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.\n".format(filename))

def get_human_DTI_graph(mode=''):
    filename = "../data/STITCH_data/human_only_"+(mode+'_' if mode else '')+"DTI_graph"
    with open(file=filename+'.pkl', mode='rb') as f:
        return pickle.load(f)


def write_human_protein_list(min_score=700,
                             mode=''):
    human_DTI_graph = get_human_DTI_graph(mode=mode)
    PPI_graph = PPI_utils.get_PPI_graph(min_score=min_score)

    print("Gathering proteins ...")
    protein_list = []
    for node in tqdm(human_DTI_graph.nodes()):
        if not node.startswith('CID') and node in list(PPI_graph.nodes()):
            protein_list.append(node)
    print("Finished.\n")

    # return sorted(PPI_utils.get_human_protein_list())

    filename = "../data/human_protein_list"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(sorted(protein_list), f, pickle.HIGHEST_PROTOCOL)

def get_human_proteins():
    filename = "../data/human_protein_list"
    with open(file=filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def write_human_prot_func_protein_list(mode=''):
    print('Loading standard proteins...')
    human_DTI_graph = get_human_DTI_graph(mode=mode)
    human_proteins = get_human_proteins()

    HPO_prots = HPO_GO_similarity.get_HPO_prot_list()

    DL2vec_path_prefix = '../data/DL2vec/DL2vec_embeddings/'
    uberon_model_filename = DL2vec_path_prefix + 'uberon_intersection_ppi_embedding'
    GO_model_filename = DL2vec_path_prefix + 'go_intersection_ppi_embedding'
    phenotype_model_filename = DL2vec_path_prefix + 'mp_intersection_ppi_embedding'

    # load model keys sets
    print('Loading models...')
    uberon_model = set(gensim.models.Word2Vec.load(uberon_model_filename).wv.vocab.keys())
    GO_model = set(gensim.models.Word2Vec.load(GO_model_filename).wv.vocab.keys())
    phenotype_model = set(gensim.models.Word2Vec.load(phenotype_model_filename).wv.vocab.keys())

    print("Gathering proteins ...")
    protein_list = []
    for node in tqdm(human_DTI_graph.nodes()):
        if not node.startswith('CID') and \
                node in human_proteins and \
                node in uberon_model and node in GO_model and node in phenotype_model and \
                node in HPO_prots:
            protein_list.append(node)
    print("Finished.\n")

    # return sorted(PPI_utils.get_human_protein_list())

    filename = "../data/human_prot_func_protein_list"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(sorted(protein_list), f, pickle.HIGHEST_PROTOCOL)

def get_human_prot_func_proteins():
    filename = "../data/human_prot_func_protein_list"
    with open(file=filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def write_human_PhenomeNET_proteins(mode=''):
    phenomeNET_proteins = PhenomeNET_DL2vec_utils.get_PhenomeNET_protein_list()
    human_DTI_graph = get_human_DTI_graph(mode=mode)
    human_proteins = get_human_proteins()

    return_protein_list = list(set(phenomeNET_proteins) & set(human_DTI_graph.nodes()) & set(human_proteins))
    return_protein_list = sorted(return_protein_list)

    filename = "../data/human_PhenomeNET_protein_list"
    with open(file=filename + '.pkl', mode='wb') as f:
        pickle.dump(return_protein_list, f, pickle.HIGHEST_PROTOCOL)

def get_human_PhenomeNET_proteins():
    filename = "../data/human_PhenomeNET_protein_list"
    with open(file=filename + '.pkl', mode='rb') as f:
        return pickle.load(f)


def get_drug_list(mode=''):
    human_DTI_graph = get_human_DTI_graph(mode=mode)
    # drug_list = sorted(list(similarity_measurement.get_SIDER_Boyce_Drubank_drug_intersection()))
    # drug_list = sorted(list(DDI_utils.get_DDI_Boyce_graph().nodes()))
    # drug_list = HPO_GO_similarity.get_HPO_SIDER_drug_list()

    drug_list = PhenomeNET_DL2vec_utils.get_PhenomeNET_drug_list()

    return np.array([drug for drug in drug_list if drug in human_DTI_graph.nodes()])


def get_side_effect_similarity_feature_list(intersect_drug_list):
    SIDER_drug_list = similarity_measurement.get_SIDER_drug_list()
    semsim_matrix = similarity_measurement.get_semantic_similarity_matrix()

    index_mapping = lambda drug: SIDER_drug_list.index(drug)

    indices = np.array(list(map(index_mapping, intersect_drug_list)))

    return semsim_matrix[indices,:][:,indices]

    # return np.array([semsim_matrix[index_mapping(drug),:] for drug in intersect_drug_list], dtype=np.float32)

def get_jaccard_side_effect_similarity_feature_list(intersect_drug_list):
    SIDER_drug_list = similarity_measurement.get_SIDER_drug_list()
    jaccard_feature_matrix = similarity_measurement.write_jaccard_se_similarity_graph()

    index_mapping = lambda drug: SIDER_drug_list.index(drug)
    indices = np.array(list(map(index_mapping, intersect_drug_list)))

    print('jaccard_feature_matrix.shape', jaccard_feature_matrix.shape)

    return jaccard_feature_matrix[indices,:][:,indices]

def get_DDI_feature_list(intersect_drug_list):
    # intersect_drug_list = get_drug_list()
    # merged_graph = DDI_utils.get_merged_DDI_graph()
    merged_graph = DDI_utils.get_DDI_Boyce_graph()

    feature_vec_list = []
    for drug in intersect_drug_list:
        feature_vector = np.zeros(len(intersect_drug_list))
        for neighbor in merged_graph.neighbors(drug):
            if neighbor in intersect_drug_list:
                feature_vector[list(intersect_drug_list).index(neighbor)] = 1
        feature_vec_list.append(feature_vector)

    return np.array(feature_vec_list)

def get_PPI_adj_mat_list(protein_list):
    # protein_list = get_human_proteins()
    protein_to_adj_mat_dict = PPI_utils.get_protein_to_adj_mat_dict()

    return np.array([protein_to_adj_mat_dict[protein] for protein in protein_list], dtype=np.int8)

def get_PPI_node_feature_mat_list(protein_list):
    # protein_list = get_human_proteins()

    protein_node_feature_dict = PPI_utils.get_protein_to_node_feature_dict()

    return np.array([protein_node_feature_dict[protein] for protein in protein_list], dtype=np.int8)

def get_PPI_dti_feature_list(drug_list, protein_list):
    human_dti_graph = get_human_DTI_graph()

    # drug_list = get_drug_list()
    # protein_list = get_human_proteins()

    protein_dti_mat = np.zeros((len(protein_list), len(drug_list)), dtype=np.int8)
    for protein_index in range(len(protein_list)):
        protein = protein_list[protein_index]
        neighbors = human_dti_graph.neighbors(protein)
        for drug_index in range(len(drug_list)):
            drug = drug_list[drug_index]
            if drug in neighbors:
                protein_dti_mat[protein_index, drug_index] = 1

    return protein_dti_mat


def get_DTIs(drug_list,
             protein_list,
             mode=''):
    DTI_graph = get_human_DTI_graph(mode=mode)

    y_data = np.zeros(len(drug_list)*len(protein_list))

    for i in range(len(protein_list)):
        protein = protein_list[i]
        if protein not in DTI_graph.nodes():
            continue
        for drug in DTI_graph.neighbors(protein):
            if drug not in drug_list:
                continue
            j = list(drug_list).index(drug)

            y_data[j * len(protein_list) + i] = 1

    return np.array(y_data, dtype=np.int8)

def get_annotated_PPI_graph():
    PPI_graph = PPI_utils.get_PPI_graph()
    node_feature_dict = PPI_utils.get_protein_to_node_feature_dict()

    # nx.set_node_attributes(PPI_graph, 'node_feature', node_feature_dict)
    for protein in PPI_graph.nodes():
        PPI_graph.node[protein]['node_feature'] = node_feature_dict[protein]

    return PPI_graph

def write_drug_to_HMM_filtered_targets_dict(rel_weight_method='wnone',
                                            alignment_method='mafft'):
    predicted_targets_dir = "../data/predicted_targets/"

    print("Parsing targets from Hmmer prediction ...")
    files = [file for file in os.listdir(predicted_targets_dir)
             if rel_weight_method in file and alignment_method in file]
    drug_filtered_targets_dict = {}
    for file in tqdm(files):
        drug = file.split('_')[0]
        target_list = []
        with open(file=predicted_targets_dir + file, mode='r') as f:
            first_line = f.readline()
            if '---' not in first_line:
                target_list.append(first_line.strip())
            for line in f:
                target_list.append(line.strip())

        drug_filtered_targets_dict[drug] = target_list

    print("Writing dict ...")
    filename = "../data/drug_to_HMM_filtered_"+alignment_method+"_"+rel_weight_method+"_targets_dict"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(drug_filtered_targets_dict, f, pickle.HIGHEST_PROTOCOL)
    print("Finished.\n")

def get_drug_to_HMM_filtered_targets_dict(rel_weight_method='wnone',
                                          alignment_method='mafft'):
    filename = "../data/drug_to_HMM_filtered_"+alignment_method+"_"+rel_weight_method+"_targets_dict"
    with open(file=filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def write_drug_protein_HMM_filtered_feature_matrix(drug_list,
                                                   protein_list,
                                                   alignment_method='mafft',
                                                   rel_weight_method='wnone'):
    return_matrix = np.zeros((len(protein_list),len(drug_list)))
    drug_HMM_filtered_targets_dict = get_drug_to_HMM_filtered_targets_dict(alignment_method=alignment_method,
                                                                           rel_weight_method=rel_weight_method)

    for protein_index, protein in enumerate(protein_list):
        for drug_index, drug in enumerate(drug_list):
            for predicted_target in drug_HMM_filtered_targets_dict[drug]:
                if predicted_target == protein:
                    return_matrix[protein_index, drug_index] = 1

    filename = "../data/HMM_filtered_"+alignment_method+"_"+rel_weight_method+"_"+"feature_matrix"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(return_matrix, f, pickle.HIGHEST_PROTOCOL)
    print("Finished.\n")

def get_drug_protein_HMM_filtered_feature_matrix(alignment_method='mafft',
                                                 rel_weight_method='wnone'):
    filename = "../data/HMM_filtered_"+alignment_method+"_"+rel_weight_method+"_"+"feature_matrix"
    with open(file=filename, mode='rb') as f:
        return pickle.load(f)

def get_protein_HMM_filtered_features(drug_list):
    # @TODO Implement this function for protein features
    # drug_toHMM_filtered_targets_dict = get_drug_HMM_filtered_features()

    print("Not implemented yet!")
    raise Exception

def get_PPI_DTI_graph_intersection():
    dti_graph = get_human_DTI_graph()
    ppi_graph = PPI_utils.get_PPI_graph()

    protein_set = set(dti_graph.nodes()) & set(ppi_graph.nodes())

    return ppi_graph.subgraph(protein_set)

def write_truncated_drug_to_SMILES_dict():
    drug_list = get_drug_list()
    drug_to_SMILES_dict = DDI_utils.get_drug_to_SMILES_dict()

    # sn_mapping = DDI_utils.get_chemical_stereo_to_normal_mapping()

    print("Writing truncated drug to SMILES dict...")

    return_dict = {drug: drug_to_SMILES_dict[drug] for drug in drug_list}

    filename = "../data/STITCH_data/truncated_drug_to_SMILES_dict"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(return_dict, f, pickle.HIGHEST_PROTOCOL)

    print("Done.")

def get_truncated_drug_to_SMILES_dict():
    filename = "../data/STITCH_data/truncated_drug_to_SMILES_dict"
    with open(file=filename + '.pkl', mode='rb') as f:
        return pickle.load(f)

def get_drughub_STRING_drug_intersection():
    drug_list = get_drug_list()

    print('Fetching drughub data...')
    drughub_drug_list = PPI_utils.get_drughub_drug_list()

    drug_intersect = set(drug_list) & set(drughub_drug_list)

    return list(drug_intersect)

def get_drughub_STRING_protein_intersection():
    protein_list = get_human_proteins()

    print('Fetching drughub data...')
    drughub_protein_list = PPI_utils.get_drughub_protein_list()

    protein_intersect = set(protein_list) & set(drughub_protein_list)

    return list(protein_intersect)

def get_DL2vec_features(entity_list):

    model_filename = "../data/PhenomeNET_data/embedding_model"
    # entities = gensim.models.Word2Vec.load(model_filename).wv.vocab.keys()
    # print('num present entities:', len(entities))


    vector_dict = gensim.models.Word2Vec.load(model_filename).wv

    return np.array([vector_dict[enti] for enti in entity_list])

def get_protein_function_embeddings(protein_list):
    DL2vec_path_prefix = '../data/DL2vec/DL2vec_embeddings/'

    uberon_model_filename = DL2vec_path_prefix + 'uberon_intersection_ppi_embedding'
    GO_model_filename = DL2vec_path_prefix + 'go_intersection_ppi_embedding'
    phenotype_model_filename = DL2vec_path_prefix + 'mp_intersection_ppi_embedding'

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
        # organism, protein_id = protein.strip().split('.')
        protein_id = protein

        uberon_embeddings.append(uberon_model[protein_id])
        GO_embeddings.append(GO_model[protein_id])
        phenotype_embeddings.append(phenotype_model[protein_id])

    uberon_embeddings = torch.Tensor(uberon_embeddings)
    GO_embeddings = torch.Tensor(GO_embeddings)
    phenotype_embeddings = torch.Tensor(phenotype_embeddings)

    return torch.stack([uberon_embeddings, GO_embeddings, phenotype_embeddings], dim=1)

def get_protein_degree_percentile(protein_list, n=10, PPI_min_score=700):
    node_degrees = PPI_utils.get_PPI_degree_for_proteins(protein_list=protein_list, PPI_min_score=PPI_min_score)


    max_degree = node_degrees.max()
    print('max_degree', max_degree)
    return_matrix = np.zeros((len(protein_list), n))

    for i, val in enumerate(np.percentile(node_degrees, np.linspace(0, 100, 100/n, endpoint=False))):
        return_matrix[:, i][node_degrees>val*max_degree/100]=1

    return return_matrix

def get_drug_to_name_mapping():
    drug_to_name_id_mapping = {} # only includes IUPAC, IUPHAR-DB, 816 (STITCH)

    filename='../data/STITCH_data/chemical.aliases.v5.0.tsv'
    print('Loading drug to IUPAC dict...')
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()

        for line in tqdm(f, total=174324327):
            drug_id, _, alias, source = line.strip().split('\t')
            if drug_id not in drug_to_name_id_mapping.keys():
                if 'IUPAC' in source or 'IUPHAR-DB' in source or '816' == source:
                    drug_to_name_id_mapping[drug_id.strip()] = alias.strip()

    return drug_to_name_id_mapping

def get_protein_to_EnsemblProtein_id():
    protein_to_Ensembl_protein_id_mapping = {}

    filename = '../data/STRING_data/9606.protein.aliases.v11.0.txt'
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()

        for line in tqdm(f, total=2224814):
            protein_id, alias, source = line.strip().split('\t')
            if source.strip() == 'Ensembl_protein_id':
                protein_to_Ensembl_protein_id_mapping[protein_id.strip()] = alias.strip()

    return protein_to_Ensembl_protein_id_mapping


def test():
    # print("DTI", len(get_human_proteins()))
    # print("PPI", len(PPI_utils.get_human_protein_list()))

    dti_graph = get_human_DTI_graph()
    print(len(dti_graph.nodes()))
    print(len(dti_graph.edges()))

    ppi_graph = PPI_utils.get_PPI_graph()
    print(len(set(dti_graph.nodes()) & set(ppi_graph.nodes())))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("--mode", type=str, default='')

    config = parser.parse_args()

    mode=config.mode
    # write_human_DTI_graph(300, mode='experimental')
    # write_human_DTI_graph(500, mode='experimental')
    # write_human_DTI_graph(700, mode=mode)

    # dti_graph = get_human_DTI_graph()
    # print("Nodes", len(dti_graph.nodes()))
    # print("Edges", len(dti_graph.edges()))

    # drug_list = get_drug_list()
    # print('Drugs:', len(drug_list))

    write_human_protein_list(mode=mode)


    # write_truncated_drug_to_SMILES_dict()
    # get_truncated_drug_to_SMILES_dict()

    # protein_list = get_human_proteins()
    # PPI_utils.write_protein_fasta(protein_list=protein_list)


    # test()

    # write_human_protein_list()
    # print('num_proteins', len(get_human_proteins()))

    # print(get_annotated_PPI_graph())

    # write_drug_to_HMM_filtered_targets_dict()


    '''
    drug_list = get_drug_list()
    protein_list = get_human_proteins()
    print(len(drug_list), len(protein_list))

    print('Fetching drughub data...')
    drughub_drug_list = PPI_utils.get_drughub_drug_list()
    drughub_protein_list = PPI_utils.get_drughub_protein_list()
    print(len(drughub_drug_list), len(drughub_protein_list))

    drug_intersect = set(drug_list) & set(drughub_drug_list)
    protein_intersect = set(protein_list) & set(drughub_protein_list)
    '''

    # write_human_prot_func_protein_list(mode=mode)
    write_human_PhenomeNET_proteins(mode=mode)
    # dicki = get_human_prot_func_proteins()
    # print('num_proteins', len(dicki))

    print('drugs', len(get_drug_list(mode=mode)))
    print('prots', len(get_human_PhenomeNET_proteins()))

    pass