import json
import pickle
import os
from joblib import Parallel, delayed

import networkx as nx

from similarity_measurement import *

import pubchempy as pcp


def get_db_PubChem_id_mapping_dict():
    filename = "../data/DDI_data/db_Pubchem_mapping_data"
    db_PubChem_id_mapping_dict = {}
    with open(file=filename, mode='r') as f:
        for line in f:
            db_id, pubchem_id = line.split(',')
            pubchem_id = pubchem_id.strip()
            if len(pubchem_id) > 8:
                continue
            pubchem_id = 'CIDm' + (8-len(pubchem_id))*'0' + pubchem_id
            db_PubChem_id_mapping_dict[db_id] = pubchem_id

    return db_PubChem_id_mapping_dict

def get_DDI_Boyce_graph():

    # db := drugbank
    DDI_graph = nx.Graph()


    filename = "../data/pddi_data/CombinedDatasetNotConservative.csv"
    with open(file=filename, mode='r', encoding='utf-8-sig') as f:
        f.readline()
        for line in f:
            split_line = line.split('\t')

            if len(split_line) < 4:
                continue

            db_id_1 = split_line[0]
            db_id_2 = split_line[2]

            if 'http://bio2rdf.org/drugbank:DB' not in db_id_1 or \
                    'http://bio2rdf.org/drugbank:DB' not in db_id_2:
                continue

            db_id_1 = db_id_1[db_id_1.find('DB'):]
            db_id_2 = db_id_2[db_id_2.find('DB'):]

            DDI_graph.add_node(db_id_1)
            DDI_graph.add_node(db_id_2)

            DDI_graph.add_edge(db_id_1, db_id_2)

    return DDI_graph

    '''
    graph_filename = "../data/pddi_data/DDI_Boyce_graph"
    with open(file=graph_filename+'pkl', mode='wb') as f:
        pickle.dump(DDI_graph, f, pickle.HIGHEST_PROTOCOL)
    '''

def write_DDI_drugbank_graph():

    DDI_graph = nx.Graph()

    mapping_dict = get_db_PubChem_id_mapping_dict()

    filename = "../data/DDI_data/DDI_data_full.json"
    raw_json = None
    with open(file=filename, mode='r') as f:
        raw_json = json.load(f)

    hit_list = raw_json['hits']['hits']

    print("Parsing drugbank data ...")
    counter = 0
    for hit in hit_list:
        print(hit)

        db_id_1 = hit['_id']

        if not hit.get('_source', None) or not hit['_source'].get('drug-interactions', None):
            continue

        for interaction in hit['_score']['drug-interactions']:
            if not interaction.get('drugbank-id', None):
                continue
            db_id_2 = interaction['drugbank-id']

            pubchem_id_1 = mapping_dict[db_id_1]
            pubchem_id_2 = mapping_dict[db_id_2]
            DDI_graph.add_node(pubchem_id_1)
            DDI_graph.add_node(pubchem_id_2)
            DDI_graph.add_edge(pubchem_id_1, pubchem_id_2)
            counter += 1
            if counter % 1000 == 0:
                print("Added edges:", counter)
    print("Finished.")

    print("Writing DDI graph extracted from Drugbank ...")
    graph_filename = "../data/DDI_data/DDI_drugbank_graph"
    with open(file=graph_filename+'pkl', mode='wb') as f:
        pickle.dump(DDI_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.".format(graph_filename))

def get_DDI_drugbank_graph():
    graph_filename = "../data/DDI_data/DDI_drugbank_graph"
    DDI_drugbank_graph = None
    with open(file=graph_filename+'.pkl', mode='rb') as f:
        DDI_drugbank_graph = pickle.load(f)



def evaluate_dicts_and_graph():


    map_1 = get_db_PubChem_id_mapping_dict()
    map_2 = get_db_PubChem_id_mapping_dict_mahmud()
    map_3 = get_pddi_db_pubchem_mapping()

    key_set1 = set(map_1.keys())
    key_set2 = set(map_2.items())
    key_set3 = set(map_3.items())

    print(len(key_set1))
    print(len(key_set2))
    print(len(key_set3))

    print(len(key_set2 & key_set3))

    super_dict = dict(key_set2 | key_set2)

    print(len(list(super_dict.keys())))


    # SIDER_only_graph = get_SIDER_only_graph()
    # drug_set = set(get_SIDER_drug_list())

    # print(len(drug_set & key_set1))




if __name__ == '__main__':
    # get_DDI_Boyce_graph()

    # evaluate_dicts_and_graph()

    write_DDI_drugbank_graph()
    # get_DDI_drugbank_graph()
