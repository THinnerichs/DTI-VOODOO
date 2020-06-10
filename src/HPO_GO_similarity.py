import numpy as np
import networkx as nx

import pickle

import similarity_measurement


def get_SIDER_to_HPO_mapping_dict():
    filename = '../data/SIDER_data/MedGen_HPO_Mapping.txt'

    mapping_dict = {}
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()

        for line in f:
            split_line = line.strip().split('|')
            se_id, HPO_id = split_line[0], split_line[1]

            mapping_dict[se_id] = HPO_id

    return mapping_dict

def get_mapped_SIDER_graph():
    mapping_dict = get_SIDER_to_HPO_mapping_dict()
    SIDER_graph = similarity_measurement.get_updated_MedDRA_label_SIDER_graph()

    print('SIDER nodes')
    counter = 0
    while counter<10:
        for node in SIDER_graph.nodes():
            if node.startswith('CID'):
                continue
            counter +=1
            print(node)


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

    print('UMLS_to_MedDRA_dict', list(UMLS_to_MedDRA_id_dict.items())[:10])

    updated_SIDER_graph = SIDER_graph

    # Add SIDER nodes to MedDRA RDF graph
    # kaust_url = rdflib.Namespace("http://www.kaust_rdf.edu.sa/drugs#")
    counter = 0
    # build annotation graph with rdf labels
    annotation_graph = nx.Graph()
    for start_node, end_node in updated_SIDER_graph.edges():
        # Switch if end_node is drug
        if 'CID' in end_node:
            start_node, end_node = end_node, start_node

        # subject = rdflib.term.URIRef(kaust_url + start_node)
        # predicate = rdflib.namespace.RDF.type

        object = UMLS_to_MedDRA_id_dict.get(end_node, None)
        if object == None:
            continue

        print('start_node', start_node)
        print('end_node', end_node)
        print('object', object)
        counter += 1

    # perform partial mapping of node labels
    return_graph = nx.relabel_nodes(SIDER_graph, mapping_dict, copy=False)
    print('Affected nodes:', len(set(mapping_dict.keys()) & set(return_graph.nodes())))
    return return_graph

def get_UniProt_prot_to_GO_function_mapping():
    gene_go_feature = {}
    with open("../data/GO_data/goa_human.gaf", "r") as f:
        for line in f.readlines():
            if line.startswith('!'):
                continue
            data = line.strip().split("\t")
            prot_id, gene_id, go_id = data[1].strip(), data[2].strip(), data[4].strip()

            evidence_score = data[6].strip()
            if not evidence_score == 'IEA' or evidence_score == 'ND':
                gene_go_feature[prot_id] = go_id



            '''
            gene_name = data[2].strip()
            go_id = data[4].strip()
            if gene_name in geneName_to_id.keys():
                if not ((evidence_score == "IEA") or (evidence_score == "ND")):
                    human_gene = geneName_to_id[gene_name]
                    if (human_gene in mouse_to_human.values()):

                        try:
                            gene_go_feature[human_gene].append(go_id)
                        except:
                            gene_go_feature[human_gene] = [go_id]
            '''
            
