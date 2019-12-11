import networkx as nx
import pickle
import rdflib


def sumyyah_request():
    filename = "../data/STITCH_data/9606.protein_chemical.links.transfer.v5.0.tsv"
    dti_graph = nx.Graph()

    print("Parsing human drug-protein-links data ...")
    with open(file=filename, mode='r') as f:
        f.readline()
        for line in f:
            split_line = line.split('\t')
            drug = split_line[0].replace('s', 'm')
            target = split_line[1]
            score = int(split_line[-1])

            dti_graph.add_node(drug)
            dti_graph.add_node(target)
            dti_graph.add_edge(drug, target, score=score)

    print("Finished.\n")

