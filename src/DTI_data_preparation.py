import numpy as np
import networkx as nx

import pickle

from similarity_measurement import get_SIDER_Boyce_Drubank_drug_intersection


def write_human_DTI_graph():
    filename = "../data/STITCH_data/9606.protein_chemical.links.transfer.v5.0.tsv"
    dti_graph = nx.Graph()

    drug_set = get_SIDER_Boyce_Drubank_drug_intersection()

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


if __name__ == '__main__':
    write_human_DTI_graph()