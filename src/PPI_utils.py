import numpy as np
import math
import networkx as nx

from tqdm import tqdm
import pickle

from joblib import Parallel, delayed
import queue
import threading
import sys
import os

from Bio import SeqIO



def prune_protein_protein_db(min_score=700,
                             mode=''):

    filename = "../data/STRING_data/9606.protein.links.full.v11.0.txt"
    target_filename = "../data/STRING_data/9606.protein.links." + str(min_score) + "_min_score.v11.0.txt"

    print("Processing raw human protein links file ...")
    p = 0.041 # see STRING documentation
    with open(file=filename, mode='r') as f, open(file=target_filename, mode='w') as targetfile:
        head = f.readline()
        targetfile.write(head)

        counter = 0

        for line in f:
            counter += 1
            if counter % 1000000 == 0:
                print("Processed lines:", counter)

            split_line = line.strip().split(' ')

            if mode=='experimental':
                experimental_score = (1-int(split_line[-6])/1000) * (1-int(split_line[-7])/1000)
                database_score = (1-int(split_line[-5])/1000) * (1-int(split_line[-4])/1000)
                experimental_score = int(1000 * (1-experimental_score * database_score))
                if experimental_score < min_score:
                    continue
                targetfile.write(split_line[0]+" "+ split_line[1]+" "+str(experimental_score)+'\n')
            else:
                total_score = int(split_line[15])/1000
                total_score_nop = (total_score-p)/(1-p)
                txt_score = int(split_line[14])/1000
                txt_score_nop = (txt_score - p)/(1-p)
                total_score_updated_nop = 1 - (1-total_score_nop)/(1-txt_score_nop)
                total_score_updated = total_score_updated_nop + p * (1-total_score_updated_nop)
                if total_score_updated * 1000 < min_score:
                    continue
                targetfile.write(split_line[0]+" "+ split_line[1]+" "+str(int(total_score_updated*1000))+'\n')
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

    print('nodes', len(PPI_graph.nodes()))
    print('edges', len(PPI_graph.edges()))

    print("Writing PPI graph to disk ...")
    graph_filename = "../data/PPI_data/PPI_graph_"+str(min_score)+"_min_score"
    with open(file=graph_filename+'.pkl', mode='wb') as f:
        pickle.dump(PPI_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing {}.\n".format(graph_filename))

def get_PPI_graph(min_score=700):
    filename = "../data/PPI_data/PPI_graph_" + str(min_score) + "_min_score"
    with open(file= filename+'.pkl', mode='rb') as f:
        return pickle.load(f)

def write_protein_fasta(protein_list):

    input_fasta_file = '../data/STITCH_data/9606.protein.sequences.v11.0.fa'

    return_sequences = []  # Setup an empty list
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        if record.id in protein_list:
            return_sequences.append(record)

    print("Found {} PPI protein sequences of {}".format(len(return_sequences), len(protein_list)))

    SeqIO.write(return_sequences, "../models/protein_representation/data/PPI_graph_protein_seqs.fasta", "fasta")

def write_drug_drughub_to_STRING_mapping():
    drughub_filename = '../data/drug_repurposing_hub/repurposing_drugs_20200324.txt'

    print('Reading drughub drugs...')
    drug_id_list = []
    with open(file=drughub_filename, mode='r') as f:
        # skip header
        for i in range(10):
            f.readline()
        for line in f:
            drug_id = line.split('\t')[0]
            drug_id_list.append(drug_id)
    print(len(drug_id_list), 'drugs found in drughub')

    alias_dict = {}
    print('Reading STITCH chemical aliases...')
    filename = '../data/STITCH_data/chemical.aliases.v5.0.tsv'
    counter = 0
    with open(file=filename, mode='r') as f:
        f.readline()
        for line in f:
            mono_id, _, alias, _ = line.strip().split('\t')
            try:
                pos = drug_id_list.index(alias)
                alias_dict[drug_id_list[pos]] = mono_id
            except:
                pass
            counter += 1
            if counter % 1000000 == 0:
                print('counter', counter, len(alias_dict.keys()))


    dict_filename = '../data/drug_repurposing_hub/drug_drughub_to_STRING_mapping'
    with open(file=dict_filename+'.pkl', mode='wb') as f:
        pickle.dump(alias_dict, f, pickle.HIGHEST_PROTOCOL)

def get_drug_drughub_to_STRING_mapping():
    dict_filename = '../data/drug_repurposing_hub/drug_drughub_to_STRING_mapping'
    with open(file=dict_filename + '.pkl', mode='rb') as f:
        return pickle.load(f)

def write_protein_drughub_to_STRING_mapping():
    filename = '../data/STRING_data/human.name_2_string.tsv'
    print("Writing drughub to STRING mapping for proteins...")
    mapping_dict = {}
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()
        for line in f:
            _, dh_id, STRING_id = line.split('\t')
            mapping_dict[dh_id] = STRING_id.strip()

    dict_filename = '../data/drug_repurposing_hub/protein_drughub_to_STRING_mapping'
    with open(file=dict_filename+'.pkl', mode='wb') as f:
        pickle.dump(mapping_dict, f, pickle.HIGHEST_PROTOCOL)

def get_protein_drughub_to_STRING_mapping():
    dict_filename = '../data/drug_repurposing_hub/protein_drughub_to_STRING_mapping'
    with open(file=dict_filename + '.pkl', mode='rb') as f:
        return pickle.load(f)

def write_drughub_dti_graph():
    drughub_filename = '../data/drug_repurposing_hub/repurposing_drugs_20200324.txt'
    drug_to_STRING_dict = get_drug_drughub_to_STRING_mapping()
    protein_to_STRING_dict = get_protein_drughub_to_STRING_mapping()

    print(len(drug_to_STRING_dict.keys()), 'drugs from drughub present.')
    print(len(protein_to_STRING_dict.keys()), 'proteins from drughub present.')

    drughub_dti_graph = nx.Graph()
    skipped_drugs = 0
    with open(file=drughub_filename, mode='r') as f:
        # Skip header
        for i in range(10):
            f.readline()
        for line in f:
            split_line = line.split('\t')
            drug, targets = split_line[0], split_line[3]
            drug = drug_to_STRING_dict.get(drug)
            if drug == None:
                skipped_drugs += 1
                continue
            for target in targets.strip().split('|'):
                target = protein_to_STRING_dict.get(target)
                if target == None:
                    continue
                drughub_dti_graph.add_edge(drug, target)

    print(skipped_drugs, 'drugs not present in mapping and skipped.')

    graph_filename = '../data/drug_repurposing_hub/drughub_dti_graph'
    with open(file=graph_filename, mode='wb') as f:
        pickle.dump(drughub_dti_graph, f, pickle.HIGHEST_PROTOCOL)

def get_drughub_dti_graph():
    graph_filename = '../data/drug_repurposing_hub/drughub_dti_graph'
    with open(file=graph_filename, mode='rb') as f:
        return pickle.load(f)

def get_drughub_drug_list():
    return list(get_drug_drughub_to_STRING_mapping().values())

def get_drughub_protein_list():
    # currently unused
    dti_graph = get_drughub_dti_graph()

    protein_list = []
    for node in dti_graph.nodes():
        if 'CID' not in node:
            protein_list.append(node)

    return protein_list

def get_PPI_degree_for_proteins(protein_list, PPI_min_score=700):
    PPI_graph = get_PPI_graph(min_score=PPI_min_score)


    return_list = []
    for protein in protein_list:
        return_list.append(len(list(PPI_graph.neighbors(protein))))

    return np.array(return_list)

def get_protein_Yamanishi_to_STITCH_mapping():
    # filename = '../data/Yamanishi_data/protein_dict_map.txt'
    filename = '../data/Yamanishi_data/uniprot_to_stitch_mapping.txt' #obtained from https://www.uniprot.org/help/uploadlists
    protein_yamanishi_to_Uniprot_mapping = {}
    print('Loading Yamanishi to Uniprot Mapping...')
    with open(file=filename, mode='r') as f:
        for line in f:
            # yamanishi_id, uniprot_id = line.strip().split(':')
            yamanishi_id, uniprot_id = line.strip().split('\t')
            protein_yamanishi_to_Uniprot_mapping[yamanishi_id] = uniprot_id
    '''
    filename = '../data/STRING_data/9606.protein.aliases.v11.0.txt'
    protein_to_STRING_mapping = {}
    with open(file=filename, mode='r') as f:
        # skip header
        f.readline()
        for line in f:
            protein_id, alias, source = line.strip().split('\t')
            if 'UniProt' in source or 'KEGG' in source:
                protein_to_STRING_mapping[alias] = protein_id
    '''

    print('yama_dict', len(protein_yamanishi_to_Uniprot_mapping))
    '''
    print('uniprot_dict', len(protein_to_STRING_mapping))
    return_dict = {prot: protein_to_STRING_mapping[protein_yamanishi_to_Uniprot_mapping[prot]] for prot
                   in protein_yamanishi_to_Uniprot_mapping.keys()
                   if protein_yamanishi_to_Uniprot_mapping[prot] in protein_to_STRING_mapping.keys()}
    '''
    return_dict = protein_yamanishi_to_Uniprot_mapping

    return return_dict


if __name__ == '__main__':
    # prune_protein_protein_db(min_score=700)

    # write_PPI_graph(min_score=700)

    # _, start = sys.argv
    # write_protein_to_subgraph_dict(cutoff=0.95)

    # write_protein_to_adj_mat_dict()
    # write_protein_to_node_feature_dict()


    # dicki = get_protein_to_node_feature_dict()
    # print(len(dicki))

    # write_drug_drughub_to_STRING_mapping()
    write_protein_drughub_to_STRING_mapping()

    write_drughub_dti_graph()

    pass
