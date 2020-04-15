import numpy as np
import networkx as nx

import pickle

import os
from tqdm import tqdm

import PPI_utils
import DDI_utils
import similarity_measurement



def write_human_DTI_graph(min_score=0):
    filename = "../data/STITCH_data/9606.protein_chemical.links.transfer.v5.0.tsv"
    dti_graph = nx.Graph()

    drug_set = similarity_measurement.get_SIDER_Boyce_Drubank_drug_intersection()

    print("Parsing human drug-protein-links data ...")
    with open(file=filename, mode='r') as f:
        f.readline()
        for line in f:
            split_line = line.split('\t')
            drug = split_line[0].replace('s','m')
            target = split_line[1]
            score = int(split_line[-1])

            if not drug in drug_set:
                continue

            dti_graph.add_node(drug)
            dti_graph.add_node(target)
            if score >= min_score:
                dti_graph.add_edge(drug, target, score=score)

    print("Finished.\n")

    print("Writing human only DTI-graph to disk ...")
    filename = "../data/STITCH_data/human_only_DTI_graph"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(dti_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.\n".format(filename))

def get_human_DTI_graph():
    filename = "../data/STITCH_data/human_only_DTI_graph"
    with open(file=filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def write_human_protein_list():
    human_DTI_graph = get_human_DTI_graph()
    protein_node_feature_dict = PPI_utils.get_protein_to_node_feature_dict()
    protein_adj_mat_dict = PPI_utils.get_protein_to_adj_mat_dict()

    print("Gathering proteins ...")
    protein_list = []
    for node in human_DTI_graph.nodes():
        if not node.startswith('CID') and \
                node in list(protein_node_feature_dict.keys()) and \
                node in list(protein_adj_mat_dict.keys()):
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

def get_drug_list():
    human_DTI_graph = get_human_DTI_graph()
    drug_list = sorted(list(similarity_measurement.get_SIDER_Boyce_Drubank_drug_intersection()))
    return np.array([drug for drug in drug_list if drug in human_DTI_graph.nodes()])


def get_side_effect_similarity_feature_list(intersect_drug_list):
    SIDER_drug_list = similarity_measurement.get_SIDER_drug_list()
    semsim_matrix = similarity_measurement.get_semantic_similarity_matrix()

    index_mapping = lambda drug: SIDER_drug_list.index(drug)

    return np.array([semsim_matrix[index_mapping(drug),:] for drug in intersect_drug_list], dtype=np.float32)

def get_DDI_feature_list(intersect_drug_list):
    # intersect_drug_list = get_drug_list()
    merged_graph = DDI_utils.get_merged_DDI_graph()

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


def get_DTIs(drug_list, protein_list):
    DTI_graph = get_human_DTI_graph()

    y_data = np.zeros(len(drug_list)*len(protein_list))

    for i in range(len(protein_list)):
        protein = protein_list[i]

        for drug in DTI_graph.neighbors(protein):
            if drug not in drug_list:
                continue
            j = list(drug_list).index(drug)

            y_data[i * len(drug_list) + j] = 1

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

    print("Writing truncated drug to SMILES dict...")
    return_dict = {drug: drug_to_SMILES_dict[drug.replace('m', 's')] for drug in drug_list}

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


def test():
    # print("DTI", len(get_human_proteins()))
    # print("PPI", len(PPI_utils.get_human_protein_list()))

    dti_graph = get_human_DTI_graph()
    print(len(dti_graph.nodes()))
    print(len(dti_graph.edges()))

    ppi_graph = PPI_utils.get_PPI_graph()
    print(len(set(dti_graph.nodes()) & set(ppi_graph.nodes())))


if __name__ == '__main__':
    # write_human_DTI_graph()

    # dti_graph = get_human_DTI_graph()
    # print("Nodes", len(dti_graph.nodes()))
    # print("Edges", len(dti_graph.edges()))

    # write_truncated_drug_to_SMILES_dict()
    # get_truncated_drug_to_SMILES_dict()

    # protein_list = get_human_proteins()
    # PPI_utils.write_protein_fasta(protein_list=protein_list)


    # test()

    # write_human_protein_list()

    # print(get_annotated_PPI_graph())

    # write_drug_to_HMM_filtered_targets_dict()


    drug_list = get_drug_list()
    protein_list = get_human_proteins()
    print(len(drug_list), len(protein_list))

    print('Fetching drughub data...')
    drughub_drug_list = PPI_utils.get_drughub_drug_list()
    drughub_protein_list = PPI_utils.get_drughub_protein_list()
    print(len(drughub_drug_list), len(drughub_protein_list))

    drug_intersect = set(drug_list) & set(drughub_drug_list)
    protein_intersect = set(protein_list) & set(drughub_protein_list)



    pass