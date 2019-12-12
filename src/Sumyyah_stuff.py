import networkx as nx
import pickle
import rdflib


def write_SIDER_only_graph():
    # Extract graph from raw SIDER2 data
    filename = "meddra_all_label_se.tsv"

    print("Extracting graph from raw data ...")
    G = nx.Graph()
    with open(file=filename, mode='r') as f:
        for line in f:
            _, flat_drug, stereo_drug, _, _, side_effect_id, side_effect_name = line.split('\t')

            flat_drug = flat_drug[:3] + "m" + flat_drug[4:]
            stereo_drug = stereo_drug[:3] + "s" + stereo_drug[4:]

            if flat_drug not in G.nodes():
                G.add_node(flat_drug)
                # G.add_node(stereo_drug)

            if side_effect_id not in G.nodes():
                G.add_node(side_effect_id)

            G.add_edge(flat_drug, side_effect_id)
            # G.add_edge(stereo_drug, side_effect_id)
    G.remove_node('')
    print("Finished.\n")

    print("Writing SIDER only graph to disc ...")
    filename = "bipartite_SIDER_only_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_SIDER_only_graph():
    print("Reading SIDER only graph ...\n")
    graph_filename = "bipartite_SIDER_only_graph"
    with open(graph_filename + '.pkl', 'rb') as f:
        return pickle.load(f)

def write_dti_graph(min_score=0):
    #
    filename = "9606.protein_chemical.links.transfer.v5.0.tsv"
    dti_graph = nx.Graph()

    print("Parsing human drug-protein-links data ...")
    with open(file=filename, mode='r') as f:
        f.readline() # skip header
        for line in f:
            split_line = line.split('\t')
            drug = split_line[0].replace('s', 'm')
            target = split_line[1]
            score = int(split_line[-1])

            dti_graph.add_node(drug)
            dti_graph.add_node(target)
            if score >= min_score:
                dti_graph.add_edge(drug, target, score=score)

    print("Finished.\n")

    print("Writing human only DTI-graph to disk ...")
    filename = "../data/STITCH_data/human_only_DTI_graph"
    with open(file=filename + '.pkl', mode='wb') as f:
        pickle.dump(dti_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.\n".format(filename))

def get_human_DTI_graph():
    filename = "../data/STITCH_data/human_only_DTI_graph"
    with open(file=filename + '.pkl', mode='rb') as f:
        return pickle.load(f)


