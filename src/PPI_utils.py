import numpy as np
import math
import networkx as nx

from tqdm import tqdm
import pickle

from joblib import Parallel, delayed
import queue
import threading



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

def get_human_protein_list(min_score=700):
    return sorted(get_PPI_graph(min_score=min_score).nodes())

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
            score = int(split_line[-1])

            PPI_graph.add_node(node_1)
            PPI_graph.add_node(node_2)
            PPI_graph.add_edge(node_1, node_2, score=score)
    print("Finished.")

    print("Writing PPI graph to disk ...")
    graph_filename = "../data/PPI_data/PPI_graph_"+str(min_score)+"_min_score"
    with open(file=graph_filename+'.pkl', mode='wb') as f:
        pickle.dump(PPI_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.\n".format(graph_filename))

def get_PPI_graph(min_score=700):
    filename = "../data/PPI_data/PPI_graph_" + str(min_score) + "_min_score"
    with open(file= filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def write_protein_to_subgraph_dict(cutoff=0.7):
    PPI_graph = get_PPI_graph()
    for node1, node2 in PPI_graph.edges():
        PPI_graph[node1][node2]['score'] = math.log(PPI_graph[node1][node2]['score']/1000, cutoff)

    print("Build protein subgraph mapping ...")
    # ego_graph_wrapper = lambda prot: nx.ego_graph(PPI_graph, prot, radius=1, center=True, undirected=True, distance='score')
    # protein_list = sorted(PPI_graph.nodes())

    # result = Parallel(n_jobs=128)(delayed(ego_graph_wrapper)(prot) for prot in tqdm(protein_list))

    # protein_subgraph_dict = dict(zip(protein_list, result))
    '''
    for protein in tqdm(PPI_graph.nodes()):
        subgraph = nx.ego_graph(PPI_graph, protein, radius=1, center=True, undirected=True, distance='score')

        protein_subgraph_dict[protein] = subgraph
    '''

    ego_graph_wrapper = lambda prot: nx.ego_graph(PPI_graph, prot, radius=1, center=True, undirected=True, distance='score')

    protein_subgraph_dict = {}

    protein_list = sorted(PPI_graph.nodes())
    filename = "../data/PPI_data/protein_to_subgraph_dict"
    round = 1
    batch_size = 256
    workers = 64
    for batch in [protein_list[i:i+batch_size] for i in range(0, len(protein_list), batch_size)]:
        print("Round {} of {}".format(round, int(len(protein_list)/batch_size+1)))
        if round > 1:
            with open(file=filename + '.pkl', mode='rb') as f:
                protein_subgraph_dict = pickle.load(f)

        print(len(list(protein_subgraph_dict.keys())))
        round += 1

        # result = Parallel(n_jobs=32)(delayed(ego_graph_wrapper)(prot) for prot in tqdm(batch))

        # batch_dict = dict(zip(batch, result))

        '''
        for protein in tqdm(batch):
            batch_dict[protein] = nx.ego_graph(PPI_graph, protein, radius=1, center=True, undirected=True, distance='score')
        '''

        q = queue.Queue()

        for protein in batch:
            q.put(protein)

        batch_dict = {}
        def worker():
            while True:
                protein = q.get()
                if protein is None:  # EOF?
                    return
                batch_dict[protein] = nx.ego_graph(PPI_graph, protein, radius=1, center=True, undirected=True,
                                                   distance='score')

        threads = [threading.Thread(target=worker) for _i in range(workers)]
        for thread in threads:
            thread.start()
            q.put(None)  # one EOF marker for each thread

        with open(file=filename + '.pkl', mode='wb') as f:
            pickle.dump({**protein_subgraph_dict, **batch_dict}, f, pickle.HIGHEST_PROTOCOL)

    print("Finished.\n")

def get_protein_to_subgraph_dict():
    filename = "../data/PPI_data/protein_to_subgraph_dict"
    with open(file=filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def write_protein_to_adj_mat_dict():
    protein_to_subgraph_dict = get_protein_to_subgraph_dict()

    max_nodes = -1
    for protein, subgraph in protein_to_subgraph_dict.items():
        if len(subgraph.nodes()) > max_nodes:
            max_nodes = len(subgraph.nodes())
    print("Maximum nodes:", max_nodes)

    print("Calculating adjacency matrices ...")
    protein_to_adj_mat_dict = {}
    protein_to_node_feature_dict = {}
    for protein, subgraph in tqdm(protein_to_subgraph_dict.items()):
        adj_mat = np.zeros((max_nodes, max_nodes))
        help_mat = nx.adjacency_matrix(subgraph, weight=None)

        adj_mat[:help_mat[0], :help_mat[1]] = help_mat

        protein_to_adj_mat_dict[protein] = adj_mat

        node_feature_mat = np.zeros((max_nodes, 1))
        help_mat = np.transpose(np.ones(len(subgraph.nodes())))
        node_feature_mat[:help_mat[0], :help_mat[1]] = help_mat
        protein_to_node_feature_dict[protein] = node_feature_mat


    print("Finished.\n")

    print("Writing protein to adjacency matrix dict ...")
    filename = "../data/PPI_data/protein_to_adj_mat_dict"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(protein_to_adj_mat_dict, f, pickle.HIGHEST_PROTOCOL)
    print("Finished.\n")

    print("Writing protein to node feature matrix dict ...")
    filename = "../data/PPI_data/protein_to_node_features_dict"
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(protein_to_node_feature_dict, f, pickle.HIGHEST_PROTOCOL)
    print("Finished.\n")

def get_protein_to_adj_mat_dict():
    filename = "../data/PPI_data/protein_to_adj_mat_dict"
    with open(file=filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def get_protein_to_node_feature_dict():
    filename = "../data/PPI_data/protein_to_node_features_dict"
    with open(file=filename+'.pkl', mode='rb') as f:
        return pickle.load(f)



if __name__ == '__main__':
    # prune_protein_protein_db(min_score=700)

    # write_PPI_graph(min_score=700)
    write_protein_to_subgraph_dict()
    # write_protein_to_adj_mat_dict()
