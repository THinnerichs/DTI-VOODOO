import numpy as np
import networkx as nx

import pickle
import subprocess

import similarity_measurement

import torch



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
    for node in SIDER_graph.nodes():
        if node.startswith('CID'):
            continue
        counter +=1
        print(node)
        if counter >= 10:
            break


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

        if counter >10:
            break

    # perform partial mapping of node labels
    return_graph = nx.relabel_nodes(SIDER_graph, mapping_dict, copy=False)
    print('Affected nodes:', len(set(mapping_dict.keys()) & set(return_graph.nodes())))
    raise Exception

    return return_graph

def get_UniProt_prot_to_GO_function_mapping():
    # useless function
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

def get_prot_to_EntrezGene_mapping():
    filename = '../data/STRING_data/9606.protein.aliases.v11.0.txt'

    mapping_dict = {}
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()
        for line in f:
            prot_id, gene_id, rest = line.split('\t')
            if 'EntrezGene' in rest:
                mapping_dict[prot_id.strip()] = gene_id.strip()
    print('Prot to EntrezGene entries:', len(mapping_dict.keys()))
    return mapping_dict

def get_gene_HPO_class_associations():
    filename = '../data/HPO_data/genes_to_phenotype.txt'
    gene_to_HPO_class_dict = {}
    print('Parsing', filename)
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()
        counter = 0
        for line in f:
            split_line = line.split('\t')
            # build mapping for entrez ids
            gene_id, HPO_term = split_line[0].strip(), split_line[2].strip()
            if not gene_to_HPO_class_dict.get(gene_id, None):
                gene_to_HPO_class_dict[gene_id] = [HPO_term]
            else:
                gene_to_HPO_class_dict[gene_id] += [HPO_term]

            # repeat for Entrez symbols
            gene_id, HPO_term = split_line[1].strip(), split_line[2].strip()
            if not gene_to_HPO_class_dict.get(gene_id, None):
                gene_to_HPO_class_dict[gene_id] = [HPO_term]
            else:
                gene_to_HPO_class_dict[gene_id] += [HPO_term]

            counter += 1
    print('Done.\n')

    return gene_to_HPO_class_dict


def write_association_file():
    # Get drug associations
    SIDER_graph = similarity_measurement.get_updated_MedDRA_label_SIDER_graph()
    MedDRA_to_HPO_mapping = get_SIDER_to_HPO_mapping_dict()

    # prune dict so it doesn't contain any additional keys, otherwise throws error
    updated_mapping = {k:v for k,v in MedDRA_to_HPO_mapping.items() if k in list(SIDER_graph.nodes())}
    premapped_side_effects = list(set(updated_mapping.keys()) & set(SIDER_graph.nodes()))
    SIDER_graph = nx.relabel_nodes(SIDER_graph, updated_mapping, copy=False)
    print('SIDER original num nodes:', len(SIDER_graph.nodes()))

    # remove side effects that got no mapping to HPO
    # remove_nodes = [node for node in SIDER_graph if not node.startswith('CID') and node not in list(MedDRA_to_HPO_mapping.keys())]
    # SIDER_graph.remove_nodes_from(remove_nodes)

    num_SIDER_drugs = len([node for node in SIDER_graph if node.startswith('CID')])
    print('Num SIDER drugs:', num_SIDER_drugs)
    print('Num SIDER side effects:', len(SIDER_graph.nodes()) - num_SIDER_drugs)

    meddra_RDF_graph = similarity_measurement.get_meddra_graph()
    # fetch mapping of MedDRA URIs to UMLS ids
    qres = meddra_RDF_graph.query(
        """SELECT DISTINCT ?UMLSid ?MedDRAid ?MedDRAParentid ?UMLSParentid
           WHERE {
              ?MedDRAid umls:cui ?UMLSid .
              ?MedDRAid rdfs:subClassOf ?MedDRAParentid .
              ?MedDRAParentid umls:cui ?UMLSParentid .
           }""")

    UMLS_id_to_UMLS_parent_dict = {}
    for UMLS_id, MedDRA_id, MedDRAParent_id, UMLS_Parentid in qres:
        UMLS_id_to_UMLS_parent_dict[UMLS_id.value] = UMLS_Parentid.value

    side_effects = [node for node in SIDER_graph if not node.startswith('CID')]
    mapped_side_effects = list(set(side_effects) & set(UMLS_id_to_UMLS_parent_dict.keys()))

    # get drugs that have at least one side effect in the remaining set
    # thus build union over neighbours of side effects as graph is bipartite
    drugs = set()
    for side_effect in premapped_side_effects:
        drugs = drugs | set(SIDER_graph.neighbors(updated_mapping[side_effect]))
    for side_effect in side_effects:
        drugs = drugs | set(SIDER_graph.neighbors(side_effect))
    drugs = list(drugs)

    print('num drugs', len(drugs))

    # some analysis
    side_effects = [node for node in SIDER_graph if not node.startswith('CID')]
    pruned_side_effects = list(set(side_effects) & set(UMLS_id_to_UMLS_parent_dict.keys()))
    print('sidies', len(set(side_effects) & set(UMLS_id_to_UMLS_parent_dict.keys())))
    print('mappies', len(set(updated_mapping.keys()) & set(UMLS_id_to_UMLS_parent_dict.values())))
    SIDER_graph = SIDER_graph.subgraph(drugs + pruned_side_effects + [updated_mapping[se] for se in premapped_side_effects])


    # get protein associations
    prot_to_gene_mapping = get_prot_to_EntrezGene_mapping()
    gene_to_HPO_mapping = get_gene_HPO_class_associations()
    prot_list = [prot for prot in prot_to_gene_mapping.keys() if gene_to_HPO_mapping.get(prot_to_gene_mapping[prot])]

    print('prots', len(prot_list))

    filename = '../data/HPO_data/HPO_prots'
    with open(file=filename+'.pkl', mode='wb') as f:
        pickle.dump(prot_list, f, pickle.HIGHEST_PROTOCOL)

    updated_drugs = set()

    # Write association file
    print('Writing association file...')
    filename = '../data/HPO_data/association_file'
    prefix = '<http://purl.obolibrary.org/obo/'
    with open(file=filename, mode='w') as f:
        # write drug associations
        for node in SIDER_graph.nodes():
            if not node.startswith('CID') and not node.startswith('HP'):
                parent_id = UMLS_id_to_UMLS_parent_dict.get(node, None)
                if parent_id:
                    parent_HPO_class = MedDRA_to_HPO_mapping.get(parent_id, None)
                    if parent_HPO_class:
                        f.write(node+' '+prefix+parent_HPO_class.replace(':','_')+'>\n')

        for node in SIDER_graph.nodes():
            if not node.startswith('CID'):
                continue

            for neighbour in SIDER_graph.neighbors(node):
                if neighbour.startswith('HP'):
                    f.write(node+' '+prefix+neighbour.replace(':','_')+'>\n')
                else:
                    f.write(node+' '+neighbour+'\n')
                updated_drugs.add(node)

        # write protein associations
        for prot in prot_list:
            for HPO_class in gene_to_HPO_mapping[prot_to_gene_mapping[prot]]:
                f.write(prot+' '+ prefix+HPO_class.replace(':','_') + '>\n')

    filename = '../data/HPO_data/HPO_SIDER_drug_list'
    with open(file=filename + '.pkl', mode='wb') as f:
        pickle.dump(list(updated_drugs), f, pickle.HIGHEST_PROTOCOL)
    print('Eventual num of drugs:', len(updated_drugs))


    print('Done.\n')

def get_HPO_SIDER_drug_list():
    filename = '../data/HPO_data/HPO_SIDER_drug_list'
    with open(file=filename + '.pkl', mode='rb') as f:
        return pickle.load(f)

def get_HPO_prot_list():
    filename = '../data/HPO_data/HPO_prots'
    with open(file=filename + '.pkl', mode='rb') as f:
        return pickle.load(f)

def write_entity_list():
    drug_list =  get_HPO_SIDER_drug_list()
    prot_list = get_HPO_prot_list()

    print('Writing entity list ...')
    filename = '../data/HPO_data/entity_list'
    with open(file=filename, mode='w') as f:
        for drug in drug_list:
            f.write(drug+'\n')
        for protein in prot_list:
            f.write(protein+'\n')

    print('Done.')

def write_PhenomeNET_association_file():


def run_DL2vec_embeddings(embedsize=200):
    # This doesn't really work as groovy-DL2vec scripts cannot handle the path handling. Just copy the command and run it from the DL2vec dir

    print('Running DL2vec embedding generator ...')
    path = '../data/HPO_data/'
    onto = path+'hp.owl'
    asso = path+'association_file'
    outfile = path+'embedding_model'
    ents = path+'entity_list'
    command = 'python ./DL2vec/runDL2vec.py -embedsize {embedsize} -ontology {onto} -associations {asso} -outfile {outfile} -entity_list {ents}'.format(embedsize=embedsize,
                                                                                                                                                        onto=onto,
                                                                                                                                                        asso=asso,
                                                                                                                                                        outfile=outfile,
                                                                                                                                                        ents=ents)

    print('Command:', command)
    subprocess.call(command, shell=True)
    print('Done.')



