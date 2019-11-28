import numpy as np
import networkx as nx

from tqdm import tqdm
import pickle



def prune_protein_protein_db(min_score=700):

    filename = "../data/STRING_data/9606.protein.links.full.v11.0.txt"
    target_filename = "../data/STRING_data/9606.protein.links." + str(min_score) + "_min_score.v11.0.txt"

    print("Processing raw human protein links file ...")
    with open(file=filename, mode='r') as f, open(file=target_filename, mode='w') as targetfile:
        targetfile.write(f.readline())

        counter = 0

        for line in f:
            counter += 1
            if counter % 1000000 == 0:
                print("Processed lines:", counter)

            if int(line.strip().split(' ')[15]) < min_score:
                continue
            targetfile.write(line)
    print("Finished.")

def write_PPI_graph(min_score=700):
    pruned_PPI_file = "../data/STRING_data/9606.protein.links." + str(min_score) + "_min_score.v11.0.txt"

    print("Building PPI graph ...")
    PPI_graph = nx.Graph()
    num_lines = sum(1 for line in open(pruned_PPI_file, 'r'))
    with open(file=pruned_PPI_file, mode='r') as f:
        f.readline() # skip header
        for line in tqdm(f, total=num_lines):
            split_line = line.split(' ')

            node_1 = split_line[0]
            node_2 = split_line[1]

            PPI_graph.add_node(node_1)
            PPI_graph.add_node(node_2)
            PPI_graph.add_edge(node_1, node_2)
    print("Finished.")

    print("Writing PPI graph to disk ...")
    graph_filename = "../data/STRING_data/PPI_graph_"+str(min_score)+"_min_score"
    with open(file=graph_filename+'.pkl', mode='wb') as f:
        pickle.dump(PPI_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.\n".format(graph_filename))

def get_PPI_graph(min_score=700):
    filename = "../data/STRING_data/PPI_graph_" + str(min_score) + "_min_score.pkl"
    with open(file= filename, mode='rb') as f:
        return pickle.load(f)



if __name__ == '__main__':
    # prune_protein_protein_db(min_score=700)

    write_PPI_graph(min_score=700)
