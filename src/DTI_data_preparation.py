import numpy as np
import networkx as nx


def write_human_DTIs():
    filename = "../data/STITCH_data/9606.protein_chemical.links.transfer.v5.0.tsv"
    dti_graph = nx.Graph()

    with open(file=filename, mode='r') as f:
        f.readline()
        for line in f:
            split_line = line.split('\t')
            drug = split_line[0]
            target = split_line[1]
            score = int(split_line[-1])

            dti_graph.add_node(drug)
            dti_graph.add_node(target)
            dti_graph.add_edge(drug, target, score=score)

