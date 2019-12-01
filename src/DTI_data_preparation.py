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

def get_human_proteins():
    human_DTI_graph = get_human_DTI_graph()
    protein_list = []
    for node in human_DTI_graph.nodes():
        if not node.startswith('CID'):
            protein_list.append(node)
    # return sorted(PPI_utils.get_human_protein_list())
    return sorted(protein_list)

def get_drug_list():
    return sorted(list(similarity_measurement.get_SIDER_Boyce_Drubank_drug_intersection()))

def get_side_effect_similarity_feature_list():
    SIDER_drug_list = similarity_measurement.get_SIDER_drug_list()
    semsim_matrix = similarity_measurement.get_semantic_similarity_matrix()

    intersect_drug_list = get_drug_list()

    index_mapping = lambda drug: SIDER_drug_list.index(drug)

    return np.array([semsim_matrix[index_mapping(drug),:] for drug in intersect_drug_list])

def get_DDI_feature_list():
    intersect_drug_list = get_drug_list()
    merged_graph = DDI_utils.get_merged_DDI_graph()

    feature_vec_list = []
    for drug in intersect_drug_list:
        feature_vector = np.zeros(len(intersect_drug_list))
        for neighbor in merged_graph.neighbors(drug):
            feature_vector[intersect_drug_list.index(neighbor)] = 1
        feature_vec_list.append(feature_vector)

    return np.array(feature_vec_list)

def get_PPI_adj_mat_list():
    protein_list = get_human_proteins()
    protein_to_adj_mat_dict = PPI_utils.get_protein_to_adj_mat_dict()

    return np.array([protein_to_adj_mat_dict[protein]
                     for protein in protein_list])

def get_PPI_node_feature_mat_list():
    protein_list = get_human_proteins()

    return np.array([1 for protein in protein_list])


def test():
    print("DTI", len(get_human_proteins()))
    print("PPI", len(PPI_utils.get_human_protein_list()))




if __name__ == '__main__':
    # write_human_DTI_graph()

    test()