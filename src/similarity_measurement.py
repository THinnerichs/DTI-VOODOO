import numpy as np
from sklearn.metrics import jaccard_similarity_score
import networkx as nx
import rdflib

import subprocess
from joblib import Parallel, delayed
import time
import queue
import threading
import os
from random import shuffle
from tqdm import tqdm
import pickle

import itertools


def get_SIDER_drug_list():
    filename = "../data/SIDER_data/meddra_all_label_se.tsv"

    drug_list = set()
    with open(file=filename, mode='r') as f:
        for line in f:
            _, flat_drug, stereo_drug, _, _, side_effect_id, side_effect_name = line.split('\t')

            flat_drug = flat_drug[:3] + "m" + flat_drug[4:]
            stereo_drug = stereo_drug[:3] + "s" + stereo_drug[4:]

            drug_list.add(flat_drug)
            # drug_list.add(stereo_drug)

    return list(drug_list)

def get_SIDER_side_effect_list():
    filename = "../data/SIDER_data/meddra_all_label_se.tsv"

    side_effect_id_list = set()
    with open(file=filename, mode='r') as f:
        for line in f:
            _, flat_drug, stereo_drug, _, _, side_effect_id, side_effect_name = line.split('\t')

            side_effect_id_list.add(side_effect_id)

    return list(side_effect_id_list)

def write_SIDER_only_graph():
    # Extract graph from raw SIDER2 data
    filename = "../data/SIDER_data/meddra_all_label_se.tsv"

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
    filename = "../data/SIDER_data/bipartite_SIDER_only_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_SIDER_only_graph():
    print("Reading SIDER only graph ...\n")
    graph_filename = "../data/SIDER_data/bipartite_SIDER_only_graph"
    with open(graph_filename + '.pkl', 'rb') as f:
        return pickle.load(f)

def write_updated_MedDRA_label_SIDER_graph():
    SIDER_only_graph = get_SIDER_only_graph()

    # Intersection between all of the below is empty
    MedDRA_delete_list, MedDRA_merge_mapping_dict, MedDRA_simple_mapping_dict = get_MedDRA_mapping()

    print("Nodes to delete:", len(MedDRA_delete_list))
    print("Nodes to merge:", len(MedDRA_merge_mapping_dict))
    print("Nodes to relabel:", len(MedDRA_simple_mapping_dict.keys()), '\n')

    # Remove deleted
    print("Removing deprecated nodes ...")
    for cui in tqdm(MedDRA_delete_list):
        if cui in SIDER_only_graph.nodes():
            SIDER_only_graph.remove_node(cui)

    print("Relabeling nodes ...")
    SIDER_only_graph = nx.relabel_nodes(SIDER_only_graph, MedDRA_simple_mapping_dict)


    def merge_nodes(G, nodes, new_node):
        """
        Merges the selected `nodes` of the graph G into one `new_node`,
        meaning that all the edges that pointed to or from one of these
        `nodes` will point to or from the `new_node`.
        """

        G.add_node(new_node)  # Add the 'merged' node

        for n in nodes:
            for neighbor in G.neighbors(n):
                G.add_edge(neighbor, new_node)

        for n in nodes:  # remove the merged nodes
            G.remove_node(n)

    print("Merging nodes ...")
    for new_cui, old_cui_list in tqdm(MedDRA_merge_mapping_dict.items()):
        merge_nodes(SIDER_only_graph, old_cui_list, new_cui)
    print("Finished.\n")

    print("Writing updated MedDRA label SIDER graph to disc ...")
    filename = "../data/MedDRA_data/updated_MedDRA_label_SIDER_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(SIDER_only_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_updated_MedDRA_label_SIDER_graph():
    print("Reading updated MedDRA label SIDER only graph ...\n")
    graph_filename = "../data/SIDER_data/updated_MedDRA_label_SIDER_graph"
    with open(graph_filename + '.pkl', 'rb') as f:
        return pickle.load(f)

def write_jaccard_se_similarity_graph():
    # SIDER only graph
    SIDER_graph = get_SIDER_only_graph()
    drug_SIDER_list = get_SIDER_drug_list()

    print("Building jaccard similarity graph ...")
    similarity_graph = nx.Graph()
    similarity_graph.add_nodes_from(drug_SIDER_list)
    for drug1, drug2 in itertools.product(drug_SIDER_list, drug_SIDER_list):
        adj_set1 = set(SIDER_graph.neighbors(drug1))
        adj_set2 = set(SIDER_graph.neighbors(drug2))

        intersec_size = len(adj_set1 & adj_set2)
        union_size = len(adj_set1 | adj_set2)
        if union_size != 0:
            jaccard_similarity = intersec_size / union_size
            similarity_graph.add_edge(drug1, drug2, weight=jaccard_similarity)
        else:
            similarity_graph.add_edge(drug1, drug2, weight=0)
    print("Finished.\n")

    print("Writing jaccard side effect similarity graph to disc ...")
    filename = "../data/similarity_results/jaccard_similarity_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(similarity_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_jaccard_se_similarity_graph():
    print("Reading jaccard side effect graph ...\n")
    graph_filename = "../data/similarity_results/jaccard_similarity_graph"
    with open(graph_filename + '.pkl', 'rb') as f:
        return pickle.load(f)

def write_meddra_graph_to_disc():
    meddra_rdf_graph_filename = "../data/MedDRA_data/MEDDRA_RDF_original.ttl"
    meddra_graph = rdflib.Graph()
    result = meddra_graph.parse(meddra_rdf_graph_filename, format='n3')
    print(result)

    print("Writing meddra RDF graph to disc ...")
    filename = "../data/MedDRA_data/meddra_RDF_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(meddra_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def write_enriched_SIDER_graph():

    # read meddra RDF graph from disc
    print("Reading meddra RDF graph ...")
    graph_filename = "../data/MedDRA_data/meddra_RDF_graph"
    meddra_RDF_graph = None
    with open(graph_filename + '.pkl', 'rb') as f:
         meddra_RDF_graph = pickle.load(f)
    print("Finished.\n")

    # fetch mapping of MedDRA URIs to UMLS ids
    qres = meddra_RDF_graph.query(
        """SELECT DISTINCT ?UMLSid ?MedDRAid
           WHERE {
              ?MedDRAid umls:cui ?UMLSid .
           }""")

    UMLS_to_MedDRA_id_dict = {}
    for UMLSid, MedDRAid in qres:
        UMLS_to_MedDRA_id_dict[UMLSid.value] = MedDRAid

    updated_SIDER_graph = get_updated_MedDRA_label_SIDER_graph()
    drug_list = get_SIDER_drug_list()
    side_effect_list = get_SIDER_side_effect_list()

    # Add SIDER nodes to MedDRA RDF graph
    kaust_url = rdflib.Namespace("http://www.kaust_rdf.edu.sa/rdf_syntax#")
    counter = 0
    # build annotation graph with rdf labels
    annotation_graph = nx.Graph()
    for start_node, end_node in updated_SIDER_graph.edges():
        # Switch if end_node is drug
        if 'CID' in end_node:
            start_node, end_node = end_node, start_node

        subject = rdflib.term.URIRef(kaust_url+start_node)
        predicate = rdflib.namespace.RDF.type

        object = UMLS_to_MedDRA_id_dict.get(end_node, None)
        if object == None:
            continue
        counter += 1

        meddra_RDF_graph.add((subject, predicate, object))

        annotation_graph.add_node(subject)
        annotation_graph.add_node(object)
        annotation_graph.add_edge(subject, object)

        if counter % 10000 == 0:
            print("Added edges:", counter)

    # Write result to disc
    print("Writing meddra RDF graph to disc ...")
    target_filename = "../data/MedDRA_data/MedDRA_enriched_SIDER_RDF_graph.ttl"
    meddra_RDF_graph.serialize(destination=target_filename, format='turtle')
    print("Finished writing ", target_filename, '\n')

    # Write annotation graph to disc
    print("Writing meddra annotation graph to disc ...")
    filename = "../data/SIDER_data/SIDER_annotation_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(annotation_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_MedDRA_mapping():
    # Read MRCUI mapping file
    MRCUI_filename = "../data/MedDRA_data/MRCUI.RRF"
    simple_mapping_dict = {}
    merge_mapping_dict = {}
    delete_list = []
    with open(file=MRCUI_filename, mode='r') as f:
        for line in f:
            old_cui, database, mode, _, _, new_cui, _, _ = line.split('|')

            # database example: '2015AB' -> 2015, 'AB'
            database_year = int(database[:4])
            database_version = database[4:]

            # only interested in changes after
            if database_year < 2015:
                continue

            if mode in ['RB', 'RO', 'RN']:
                simple_mapping_dict[old_cui] = new_cui
            elif mode == 'SY':
                if merge_mapping_dict.get(new_cui, None):
                    merge_mapping_dict[new_cui].append(old_cui)
                else:
                    merge_mapping_dict[new_cui] = []
            elif mode == 'DEL':
                delete_list.append(old_cui)

    # Intersection between all of the below is empty
    return delete_list, merge_mapping_dict, simple_mapping_dict

def write_annotation_file():
    """
    Write annotation file for groovy script
    :return:
    """

    # read meddra RDF graph from disc
    print("Reading annotation graph ...")
    filename = "../data/SIDER_data/SIDER_annotation_graph"
    annotation_graph = None
    with open(filename + '.pkl', 'rb') as f:
         annotation_graph = pickle.load(f)
    print("Finished.\n")

    drug_list = get_SIDER_drug_list()

    print("Writing annotation file ...")
    annotation_file = "../data/annotation_file_for_groovy.tsv"
    with open(file=annotation_file, mode='w') as f:
        for drug in drug_list:
            present = False
            drug_rdf_term = None
            for rdf_id in annotation_graph.nodes():
                if drug in rdf_id:
                    present = True
                    drug_rdf_term = rdf_id
                    break
            if not present:
                continue
            neighbor_list = annotation_graph.neighbors(drug_rdf_term)
            f.write(drug+'\t'+'\t'.join(neighbor_list)+'\n')
    print("Finished.")

def execute_DD_semantic_similarity_matrix():
    """
    Wrapper for executing of external groovy script
    :return:
    """
    command = "groovy -cp './lib/*' SemanticSimilarity.groovy"
    print("Running {} ...".format(command))
    subprocess.call(command, shell=True)
    print("Finished.")

def get_semantic_similarity_matrix():

    drug_list = get_SIDER_drug_list()
    num_drugs = len(drug_list)
    semsim_matrix = np.zeros((len(drug_list), len(drug_list)), dtype=np.float)

    print("Parsing semantic similarity matrix ...")
    score_list_filename = "../data/similarity_results/semsim_drugs.txt"
    num_lines = sum(1 for line in open(score_list_filename, 'r'))
    with open(file=score_list_filename, mode='r') as f:
        counter = 0
        for line in tqdm(f, total=num_lines):
            i = int(counter/num_drugs)
            j = counter%num_drugs
            semsim_matrix[i,j] = float(line.strip())

            counter += 1
    print("Finished.\n")

    return semsim_matrix


if __name__ == '__main__':
    # write_SIDER_only_graph()
    # write_jaccard_se_similarity_graph()
    # get_jaccard_se_similarity_graph()

    # write_updated_MedDRA_label_SIDER_graph()
    # get_updated_MedDRA_label_SIDER_graph()

    # write_enriched_SIDER_graph()
    # write_annotation_file()

    execute_DD_semantic_similarity_matrix()

    # get_semantic_similarity_matrix()

