import numpy as np

import json
import pickle
from tqdm import tqdm

import networkx as nx

import similarity_measurement



def get_db_PubChem_id_mapping_dict():
    STITCH_mapping = get_STITCH_db_Pubchem_mapping_dict()
    drugbank_mapping = get_drugbank_db_PubChem_id_mapping_dict()
    dhimmel_dict = get_dhimmel_db_to_Pubchem_mapping_dict()

    return {**STITCH_mapping, **drugbank_mapping, **dhimmel_dict}


def get_drugbank_db_PubChem_id_mapping_dict():
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

def write_STITCH_db_Pubchem_mapping_dict():
    filename = "../data/STITCH_data/chemical.aliases.v5.0.tsv"
    num_lines = 174324327
    db_pubchem_mapping_dict = {}
    with open(file=filename, mode='r') as f:
        for line in tqdm(f, total=num_lines):
            drug, stereo, alias, source = line.split('\t')
            if alias.startswith('DB') and not alias.startswith('DB-'):
                db_pubchem_mapping_dict[alias] = drug

    dict_filename = "../data/STITCH_data/STITCH_drugbank_pubchem_mapping_dict"
    with open(file=dict_filename+'.pkl', mode='wb') as f:
        pickle.dump(db_pubchem_mapping_dict, f, pickle.HIGHEST_PROTOCOL)

def get_dhimmel_db_to_Pubchem_mapping_dict():
    filename = '../data/DDI_utils/pubchem-mapping.tsv'

    return_dict = {}
    with open(file=filename, mode='r') as f:
        f.readline()
        for line in f:
            db_id, pubchem_id = line.strip().split('\t')
            pubchem_id = 'CIDm' + (8-len(pubchem_id))*'0' + pubchem_id
            return_dict[db_id] = pubchem_id

    return return_dict

def get_Yamanishi_db_to_PubChem_mapping_dict():
    filename = '../data/Yamanishi_data/drug_mapping.txt'
    stereo_to_mono_mapping = get_chemical_stereo_to_normal_mapping()

    return_dict = {}

    with open(file=filename, mode='r') as f:
        for line in f:
            if len(line) <= 9:
                continue
            db_id, pubchem_id = line.strip().split('\t')
            if len(pubchem_id) > 8:
                try:
                    pubchem_id = stereo_to_mono_mapping['CIDs' + pubchem_id [1:]]
                except:
                    continue
                    print(pubchem_id)
            else:
                pubchem_id = 'CIDm' + (8-len(pubchem_id))*'0' + pubchem_id

            return_dict[db_id] = pubchem_id

    return return_dict

def get_STITCH_db_Pubchem_mapping_dict():
    dict_filename = "../data/STITCH_data/STITCH_drugbank_pubchem_mapping_dict"
    with open(file=dict_filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def get_DDI_Boyce_graph():

    # db := drugbank
    DDI_graph = nx.Graph()

    mapping_dict = get_STITCH_db_Pubchem_mapping_dict()

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

            pubchem_id_1 = mapping_dict.get(db_id_1, None)
            pubchem_id_2 = mapping_dict.get(db_id_2, None)

            if not pubchem_id_1 or not pubchem_id_2:
                continue

            DDI_graph.add_node(pubchem_id_1)
            DDI_graph.add_node(pubchem_id_2)

            DDI_graph.add_edge(pubchem_id_1, pubchem_id_2)

    return DDI_graph

def write_DDI_drugbank_graph():

    DDI_graph = nx.Graph()

    mapping_dict = get_db_PubChem_id_mapping_dict()

    filename = "../data/DDI_data/DDI_data_full.json"
    raw_json = None
    print("Loading JSON ...")
    with open(file=filename, mode='r') as f:
        raw_json = json.load(f)
    print("Finished.\n")

    hit_list = raw_json['hits']['hits']

    print("Parsing drugbank data ...")
    counter = 0
    drugs_missing_in_dict_counter = 0
    for hit in hit_list:
        db_id_1 = hit['_id']

        if not hit.get('_source', None) or not hit['_source'].get('drug-interactions', None):
            continue

        for interaction in hit['_source']['drug-interactions']:
            if not interaction.get('drugbank-id', None):
                continue
            db_id_2 = interaction['drugbank-id']

            pubchem_id_1 = mapping_dict.get(db_id_1, None)
            pubchem_id_2 = mapping_dict.get(db_id_2, None)

            if not pubchem_id_1 or not pubchem_id_2:
                drugs_missing_in_dict_counter += 1
                continue

            DDI_graph.add_node(pubchem_id_1)
            DDI_graph.add_node(pubchem_id_2)
            DDI_graph.add_edge(pubchem_id_1, pubchem_id_2)
            counter += 1
            if counter % 100000 == 0:
                print("Added edges:", counter)
    print("Finished.")
    print("Skipped edges:", drugs_missing_in_dict_counter, '\n')

    print("Writing DDI graph extracted from Drugbank ...")
    graph_filename = "../data/DDI_data/DDI_drugbank_graph"
    with open(file=graph_filename+'.pkl', mode='wb') as f:
        pickle.dump(DDI_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.\n".format(graph_filename))

def get_DDI_drugbank_graph():
    graph_filename = "../data/DDI_data/DDI_drugbank_graph"
    with open(file=graph_filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def get_merged_DDI_graph():

    merged_graph = get_DDI_drugbank_graph()
    boyce_graph = get_DDI_Boyce_graph()

    merged_graph.add_nodes_from(boyce_graph.nodes())
    merged_graph.add_edges_from(boyce_graph.edges())

    return merged_graph

def write_drug_to_SMILES_dict():
    filename = '../data/STITCH_data/chemicals.v5.0.tsv'

    print("Writing drug to SMILES dict...")

    drug_to_smiles_dict = {}
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()

        for line in f:
            drug_id, _, _, drug_smiles_enc = line.strip().split('\t')
            drug_to_smiles_dict[drug_id] = drug_smiles_enc.strip()

    with open(file='../data/STITCH_data/drug_to_SMILES_dict.pkl', mode='wb') as f:
        pickle.dump(drug_to_smiles_dict, f, pickle.HIGHEST_PROTOCOL)

    print("Done.")

def get_drug_to_SMILES_dict():
    with open(file='../data/STITCH_data/drug_to_SMILES_dict.pkl', mode='rb') as f:
        return pickle.load(f)

def write_chemical_stereo_to_normal_mapping():
    filename = "../data/STITCH_data/chemical.aliases.v5.0.tsv"
    return_mapping_dict = {}
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()

        for line in tqdm(f, total=174324327):
            mono, stereo, _, _ = line.strip().split('\t')
            return_mapping_dict[stereo.strip()] = mono.strip()

    filename = '../data/STITCH_data/chemical_stereo_to_normal_mapping_dict'
    with open(file=filename, mode='wb') as f:
        pickle.dump(return_mapping_dict, f, pickle.HIGHEST_PROTOCOL)

def get_chemical_stereo_to_normal_mapping():
    filename = '../data/STITCH_data/chemical_stereo_to_normal_mapping_dict'
    with open(file=filename, mode='rb') as f:
        return pickle.load(f)

def get_yamanishi_drug_list():
    yamanishi_drug_mapping = get_Yamanishi_db_to_PubChem_mapping_dict()

    path = '../data/Yamanishi_data/'

    drug_list = []
    with open(file=path + 'drug.txt', mode='r') as f:
        for line in f:
            drug_list.append(line.strip())

    drug_list = np.array(list(map(lambda d: yamanishi_drug_mapping.get(d, None), drug_list)))

    return drug_list

def get_yamanishi_drug_side_effects():

    path = '../data/Yamanishi_data/'

    # build drug-side effect matrix
    filename = 'mat_drug_se.txt'
    drug_side_effect_matrix = []
    with open(file=path + filename, mode='r') as f:
        for line in f:
            drug_side_effect_matrix.append(list(map(int, list(line.strip().replace(' ', '')))))
    drug_side_effect_matrix = np.array(drug_side_effect_matrix)
    print('drug_side_effect_matrix.shape', drug_side_effect_matrix.shape)

    return drug_side_effect_matrix

def get_yamanishi_side_effect_annotations():
    path = '../data/Yamanishi_data/'

    # parse side effect HPO_term to name mapping
    mapping_dict = {}
    for db in ['hp','mp']:
        filename = f'{db}_for_synonyms.obo'
        with open(file=path+filename, mode='r') as f:

            id = None
            for line in f:
                line = line.strip()
                if line.startswith('id'):
                    id = line.split(' ')[-1]
                elif line.startswith('synonym'):
                    split_line = line.split('"')
                    name = split_line[1].lower()
                    mapping_dict[name] = id
                elif line.startswith('name'):
                    split_line = line.split(' ')

                    name = ' '.join(split_line[1:]).lower()
                    mapping_dict[name] = id

    # parse side effect names
    filename = 'se.txt'
    side_effect_list = []
    with open(file=path+filename, mode='r') as f:
        for line in f:
            side_effect_list.append(line.strip().lower())

    mapped_list = list(map(lambda se: mapping_dict.get(se, None), side_effect_list))

    print('total side effects:', len(mapped_list))
    print('Side effects with no mapping:', mapped_list.count(None))

    return np.array(mapped_list)







def evaluate_dicts_and_graph():

    '''
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
    '''

    drugbank_graph = get_DDI_drugbank_graph()
    boyce_graph = get_DDI_Boyce_graph()

    # print(len(drugbank_graph.nodes()))
    print('boyce', len(boyce_graph.nodes()), len(boyce_graph.edges()))
    print('db', len(drugbank_graph.nodes()), len(drugbank_graph.edges()))

    intersect = set(drugbank_graph.nodes()) & set(boyce_graph.nodes())
    print('intersect', len(intersect))

    merged_graph = get_merged_DDI_graph()
    print('merged', len(merged_graph.nodes()), len(merged_graph.edges()))

    SIDER_drugs = set(similarity_measurement.get_SIDER_drug_list())

    print(len(SIDER_drugs | intersect))


    # SIDER_only_graph = get_SIDER_only_graph()
    # drug_set = set(get_SIDER_drug_list())

    # print(len(drug_set & key_set1))





if __name__ == '__main__':
    # write_SITCH_db_Pubchem_mapping_dict()
    # get_STITCH_db_Pubchem_mapping_dict()

    # get_DDI_Boyce_graph()

    # write_DDI_drugbank_graph()
    # print(len(similarity_measurement.get_SIDER_Boyce_Drubank_drug_intersection()))

    evaluate_dicts_and_graph()

    # get_DDI_drugbank_graph()

    # print(get_merged_DDI_graph().nodes())

    # write_drug_to_SMILES_dict()

    pass
