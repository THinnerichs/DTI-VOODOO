import numpy as np
import networkx as nx

import pickle

import PPI_utils
import DDI_utils
import similarity_measurement



def write_human_DTI_graph():
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

    # intersect_drug_list = get_drug_list()

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

def test():
    # print("DTI", len(get_human_proteins()))
    # print("PPI", len(PPI_utils.get_human_protein_list()))
    dti_graph = get_human_DTI_graph()
    print(len(dti_graph.nodes()))
    print(len(dti_graph.edges()))




if __name__ == '__main__':
    # write_human_DTI_graph()

    test()

    # write_human_protein_list()

    # print(get_annotated_PPI_graph())
    pass