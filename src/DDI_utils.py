import json
import pickle
import os

from similarity_measurement import *


def write_db_PubChem_id_mapping_dict():

    filename = "../data/DDI_data/DB_PubChem_mapping_data_full.json"
    raw_json = None
    with open(file=filename, mode='r') as f:
        raw_json = json.load(f)

    counter = 0
    db_PubChem_id_dict = {}
    for hit in raw_json['hits']['hits']:

        if hit['_type'] == 'drugbank':
            if not hit.get('_source', None):
                continue
            if not hit['_source'].get('external-identifiers', None):
                continue

            for ex_id in hit['_source']['external-identifiers']:
                if ex_id['resource'] == "PubChem Compound":
                    db_id = hit['_id']
                    pubchem_id = ex_id['identifier']
                    if len(pubchem_id) > 8:
                        continue
                    pubchem_id = "CIDm" + (8-len(pubchem_id)) * '0' + pubchem_id

                    db_PubChem_id_dict[db_id] = pubchem_id

    print("Writing  to disk ...")
    filename = "../data/DDI_data/db_pubchem_mapping_dict_old"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(db_PubChem_id_dict, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_db_PubChem_id_mapping_dict():
    print("Reading jaccard side effect graph ...\n")
    graph_filename = "../data/DDI_data/db_pubchem_mapping_dict_old"
    with open(graph_filename + '.pkl', 'rb') as f:
        return pickle.load(f)

def get_db_PubChem_id_mapping_dict_mahmud():
    filename = "../data/DDI_data/db_Pubchem_mapping_data_2"
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

def evaluate_dicts_and_graph():


    map_1 = get_db_PubChem_id_mapping_dict()
    map_2 = get_db_PubChem_id_mapping_dict_mahmud()

    key_set1 = set(map_1.keys())
    key_set2 = set(map_2.values())

    print(len(key_set1))
    print(len(key_set2))


    SIDER_only_graph = get_SIDER_only_graph()
    drug_set = set(get_SIDER_drug_list())

    print(len(drug_set & key_set1))

def read_pddi_Boyce_data():

    filename = "../data/pddi_data/CombinedDatasetNotConservative.csv"
    db_id_name_dict = {}
    with open(file=filename, mode='r', encoding='utf-8-sig') as f:
        print(f.readline())
        for line in f:
            split_line = line.split('\t')

            if len(split_line) < 4:
                continue

            db_id_1 = split_line.pop(0)
            db_name_1 = split_line.pop(0)
            db_id_2 = split_line.pop(0)
            db_name_2 = split_line.pop(0)

            if 'http://bio2rdf.org/drugbank:DB' not in db_id_1:
                continue
            db_id_name_dict[db_id_1] = db_name_1

            if 'http://bio2rdf.org/drugbank:DB' not in db_id_2:
                continue
            db_id_name_dict[db_id_2] = db_name_2
    print(len(list(db_id_name_dict.keys())))

    for 





def test_pubchempy_search():
    import pubchempy as pcp

    CID_list = pcp.get_cids('alprazolam', 'name', 'substance', list_return='flat')

    print(CID_list)





if __name__ == '__main__':
    # write_db_PubChem_id_mapping_dict()

    # get_db_PubChem_id_mapping_dict_mahmud()

    # evaluate_dicts_and_graph()

    # test_pubchempy_search()

    read_pddi_Boyce_data()

