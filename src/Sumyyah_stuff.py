import networkx as nx
import pickle
import rdflib

from tqdm import tqdm


def write_SIDER_only_graph():
    # Extract graph from raw SIDER2 data
    filename = "/data/meddra_all_label_se.tsv"

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
    filename = "/data/bipartite_SIDER_only_graph"
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
    filename = "/data/9606.protein_chemical.links.transfer.v5.0.tsv"
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
    filename = "/data/human_only_DTI_graph"
    with open(file=filename + '.pkl', mode='wb') as f:
        pickle.dump(dti_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.\n".format(filename))

def get_human_DTI_graph():
    filename = "/data/human_only_DTI_graph"
    with open(file=filename + '.pkl', mode='rb') as f:
        return pickle.load(f)

def write_human_DTI_graph_to_RDF():
    dti_graph = get_human_DTI_graph()

    dti_RDF_graph = rdflib.Graph()
    kaust_drug_url = rdflib.Namespace("http://www.kaust_dti_rdf.edu.sa/drugs#")
    kaust_protein_url = rdflib.Namespace("http://www.kaust_dti_rdf.edu.sa/protein#")
    counter = 0
    for start_node, end_node in dti_graph.edges():
        # Switch if end_node is drug
        if 'CID' in end_node:
            start_node, end_node = end_node, start_node

        subject = rdflib.term.URIRef(kaust_drug_url+start_node)
        predicate = rdflib.term.URIRef(kaust_drug_url+"interacts_with")

        object = rdflib.term.URIRef(kaust_protein_url + end_node)
        if object == None:
            continue
        counter += 1

        dti_RDF_graph.add((subject, predicate, object))

        if counter % 10000 == 0:
            print("Added edges:", counter)

    print("Writing meddra RDF graph to disc ...")
    target_filename = "/data/SIDER_only_RDF.ttl"
    dti_RDF_graph.serialize(destination=target_filename, format='turtle')
    print("Finished writing ", target_filename, '\n')


def get_MedDRA_mapping():
    # Read MRCUI mapping file
    MRCUI_filename = "/data/MRCUI.RRF"
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
    filename = "/data/updated_MedDRA_label_SIDER_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(SIDER_only_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_updated_MedDRA_label_SIDER_graph():
    print("Reading updated MedDRA label SIDER only graph ...\n")
    graph_filename = "../data/SIDER_data/updated_MedDRA_label_SIDER_graph"
    with open(graph_filename + '.pkl', 'rb') as f:
        return pickle.load(f)

def write_meddra_graph_to_disk():
    meddra_rdf_graph_filename = "/data/MEDDRA_RDF_original.ttl"
    meddra_graph = rdflib.Graph()
    result = meddra_graph.parse(meddra_rdf_graph_filename, format='n3')

    print("Writing meddra RDF graph to disc ...")
    filename = "../data/MedDRA_data/meddra_RDF_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(meddra_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def write_enriched_SIDER_graph():

    # read meddra RDF graph from disc
    print("Reading meddra RDF graph ...")
    graph_filename = "/data/meddra_RDF_graph"
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

    # Add SIDER nodes to MedDRA RDF graph
    kaust_url = rdflib.Namespace("http://www.kaust_rdf.edu.sa/drugs#")
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
    target_filename = "/data/MedDRA_enriched_SIDER_RDF_graph.ttl"
    meddra_RDF_graph.serialize(destination=target_filename, format='turtle')
    print("Finished writing ", target_filename, '\n')

    # Write annotation graph to disc
    print("Writing meddra annotation graph to disc ...")
    filename = "/data/SIDER_annotation_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(annotation_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')


if __name__ == '__main__':
    write_SIDER_only_graph()
    write_dti_graph(min_score=0)
    write_human_DTI_graph_to_RDF()
    write_updated_MedDRA_label_SIDER_graph()
    write_meddra_graph_to_disk()
    write_enriched_SIDER_graph()
