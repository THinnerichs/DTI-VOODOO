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


def test_biopython_PairwiseAligner():
    # initialize and tune aligner
    from Bio import Align, SeqIO
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    # options:
    # STR, RISLER, RAO, PAM70, PAM30, PAM250, MDM78, MCLACHLAN, LEVIN, JONES, JOHNSON, GONNET1992, GENETIC,
    # FENG, DAYHOFF, BLOSUM,

    start_time = time.time()

    for matrix in ['STR', 'RISLER', 'RAO', 'PAM70', 'PAM30', 'PAM250', 'MDM78', 'MCLACHLAN', 'LEVIN', 'JONES', 'JOHNSON', 'GONNET1992', 'GENETIC',
                   'FENG', 'DAYHOFF', 'BLOSUM90', 'BENNER22', 'BENNER6', 'BENNER74', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80']:

        print("\n" + matrix)
        aligner.substitution_matrix = substitution_matrices.load(matrix)

        with open(file="../data/score_list", mode='a') as f:
            f.write("\nMATRIX: "+matrix+"\n")

        alignment_path = "../data/alignment_targets/"
        database_filename = "CIDm00000003_kalign_aligned_800_min_score.afa"
        query_filename = "CIDm00000006_kalign_aligned_800_min_score.afa"

        # make_blast_db
        # fasta_path = "../data/fasta_files/"
        # filename = "CIDm00000043_fasta_800_min_score.fasta"
        database_fasta_file = alignment_path + database_filename
        query_fasta_file = alignment_path + query_filename

        database_records = list(SeqIO.parse(database_fasta_file, 'fasta'))
        query_records = list(SeqIO.parse(query_fasta_file, 'fasta'))

        print(len(database_records))
        print(len(query_records))


        def help_func(doublet):
            with open(file="../data/score_list", mode='a') as f:
                seq1 = str(doublet[0].seq).replace('-', '-')
                seq2 = str(doublet[1].seq).replace('-', '-')
                f.write(str(aligner.score(seq1, seq2))+"\n")
        # help_func = lambda doublet: aligner.score(doublet[0].seq, doublet[1].seq)

        score_list = Parallel(n_jobs=40)(delayed(help_func)(doublet) for doublet in tqdm(list(itertools.product(database_records[:20], query_records[:20]))))


        # score_list = np.array(score_list)

        # np.save("../data/score_list_nparray.npy", score_list)

    print("This took {} seconds.".format(time.time()-start_time))

def replace_gap_symbols_in_alignment(samples=500):
    from Bio import SeqIO

    alignment_path = "../data/alignment_targets/"
    database_filename = "CIDm00000003_kalign_aligned_800_min_score.afa"
    query_filename = "CIDm00000006_kalign_aligned_800_min_score.afa"

    with open(file=alignment_path+database_filename[:-3]+"fasta", mode='w') as f:
        for record in list(SeqIO.parse(alignment_path+database_filename, 'fasta'))[:samples]:
            f.write(">"+str(record.id)+"\n")
            f.write(str(record.seq).replace('-', 'X')+"\n")

    with open(file=alignment_path + query_filename[:-3] + "fasta", mode='w') as f:
        for record in list(SeqIO.parse(alignment_path + query_filename, 'fasta'))[:samples]:
            f.write(">"+str(record.id)+"\n")
            f.write(str(record.seq).replace('-', 'X')+'\n')


def test_blast():

    # parameters
    threads = 40

    # make_blast_db
    # fasta_path = "../data/fasta_files/"
    # filename = "CIDm00000043_fasta_800_min_score.fasta"


    alignment_path = "../data/alignment_targets/"
    database_filename = "CIDm00000003_kalign_aligned_800_min_score.fasta"
    query_filename = "CIDm00000003_kalign_aligned_800_min_score.fasta"

    # subprocess.call("cp "+alignment_path+database_filename+" "+ alignment_path+database_filename[:-3]+"fasta")
    # subprocess.call("cp "+alignment_path+query_filename+" "+ alignment_path+query_filename[:-3]+"fasta")

    drug_name = database_filename.split("_")[0]
    database_name = "../data/" + drug_name + "_blast_db"

    # Build blast db
    print("Building database ...")
    command = "./makeblastdb -dbtype 'prot' "+\
              "-input_type 'fasta' "+\
              "-in " + alignment_path+database_filename + " "+\
              "-out " + database_name
    print(command)
    subprocess.call(command, shell=True)
    print("Finished.\n")


    print("Running query ...")
    results_filename = "../data/"+drug_name+"_blast_result.xml"

    blast_command = "./blastp "+\
                    "-task blastp-fast "+\
                    "-num_threads "+str(threads)+" "+\
                    "-query "+alignment_path+query_filename+" "+\
                    "-db "+database_name+" "+\
                    "-out "+results_filename+" "+\
                    "-outfmt 5"
                    # "-evalue 0.05 "+\
    print(blast_command)

    subprocess.call(blast_command, shell=True)
    print("Finished.\n")

def evaluate_Blast_XML():
    drug_name = "CIDm00000003"
    results_filename = "../data/"+drug_name+"_blast_result.xml"

    from Bio.Blast import NCBIXML
    E_VALUE_THRESH = 100

    '''
    num_lines = None
    fasta_path = "../data/fasta_files/"
    filename = "CIDm00000003_fasta_800_min_score.fasta"
    with open(file=fasta_path+filename, mode='r') as f:
        num_lines = sum(1 for line in f)

    print("Lines:", num_lines)
    '''



    start_time = time.time()
    print("Parsing similarity scores ...")
    # prot_prot_sim_dict = {}
    score_list = []
    for record in NCBIXML.parse(open(results_filename)):
        if record.alignments: # skip queries with no   matches
            # prot_prot_sim_dict[record.query] = {}
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        # prot_prot_sim_dict[record.query][alignment.title] = int(hsp.score)
                        score_list.append(int(hsp.score))

    print(score_list)
    # score_list = np.array(score_list)

    # print(score_list.mean(), score_list.std())

    print("Finished in {} seconds.".format(time.time() - start_time))

def run_similarity_pipeline(threads=8,
                            e_value_threshold=0.05,
                            min_score=800):

    alignment_path = "../data/alignment_targets/"

    def similarity_score(file1, file2):
        # build filenames
        database_alignment_file = alignment_path + file1
        query_alignment_file = alignment_path + file2

        database_drug_name = file1.split("_")[0]
        query_drug_name = file2.split("_")[0]
        database_name = alignment_path + database_drug_name + "_blast_db"

        # Build blast db
        print("Building database ...")
        command = "./makeblastdb -dbtype 'prot' " + \
                  "-in " + database_alignment_file + " " + \
                  "-out " + database_name
        print(command)
        subprocess.call(command, shell=True)
        print("Finished.\n")

        # run blast query
        print("Running query ...")
        results_filename = "" + database_drug_name + "_" + query_drug_name + "_blast_result.xml"

        blast_command = "./blastp " + \
                        "-task blastp-fast " + \
                        "-num_threads 32 " + \
                        "-query " + query_alignment_file + " " + \
                        "-db " + database_name + " " + \
                        "-out " + results_filename + " " + \
                        "-evalue " + str(e_value_threshold)+" "\
                        "-outfmt 5"
        print(blast_command)

        subprocess.call(blast_command, shell=True)
        print("Finished.\n")

        # Evaluate query
        from Bio.Blast import NCBIXML

        # Count lines in both files for normalization of results
        file1_num_lines = None
        database_fasta_file = "../data/fasta_files/" + database_drug_name + "_fasta_" + str(min_score) + "_min_score.afa"
        with open(file=database_fasta_file, mode='r') as f:
            file1_num_lines = sum(1 for line in f)

        file2_num_lines = None
        query_fasta_file = "../data/fasta_files/" + query_drug_name + "_fasta_" + str(min_score) + "_min_score.afa"
        with open(file=query_fasta_file, mode='r') as f:
            file2_num_lines = sum(1 for line in f)

        # Parse xml file
        start_time = time.time()
        print("Parsing similarity scores ...")
        score_sum = 0
        for record in NCBIXML.parse(open(results_filename)):
            if record.alignments:  # skip queries with no   matches
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < e_value_threshold:
                            score_sum += int(hsp.score)

        print("Finished in {} seconds.".format(time.time() - start_time))
        final_score = score_sum / (file1_num_lines/2 + file2_num_lines/2)

        evaluation_results_filename = "../data/similarity_results"
        with open(file=evaluation_results_filename, mode='a') as f:
            f.write(query_drug_name+"\t"+database_drug_name+"\t"+str(final_score))

    q = queue.Queue()

    help_file_doublets = itertools.product(os.listdir(alignment_path), os.listdir(alignment_path))

    # omit symmetric doublets
    file_doublets=[]
    for doublet in help_file_doublets:
        file1, file2 = doublet
        if (file2, file2) not in file_doublets and not file1 == file2:
            file_doublets.append(doublet)

    # initialize queue
    for doublet in file_doublets:
        q.put(doublet)

    # define worker
    def worker():
        while True:
            doublet = q.get()
            if doublet is None:  # EOF?
                return
            file1, file2 = doublet

            similarity_score(file1, file2)

    # Start the workers
    max_overall_threads = 50
    num_threads = [threading.Thread(target=worker) for _i in range(int(max_overall_threads/threads))]
    for thread in num_threads:
        thread.start()
        q.put(None)  # one EOF marker for each thread

def get_SIDER_drug_list():
    filename = "../data/meddra_all_label_se.tsv"

    drug_list = set()
    with open(file=filename, mode='r') as f:
        for line in f:
            _, flat_drug, stereo_drug, _, _, side_effect_id, side_effect_name = line.split('\t')

            flat_drug = flat_drug[:3] + "m" + flat_drug[4:]
            stereo_drug = stereo_drug[:3] + "s" + stereo_drug[4:]

            drug_list.add(flat_drug)
            drug_list.add(stereo_drug)

    return list(drug_list)

def get_SIDER_side_effect_list():
    filename = "../data/meddra_all_label_se.tsv"

    side_effect_id_list = set()
    with open(file=filename, mode='r') as f:
        for line in f:
            _, flat_drug, stereo_drug, _, _, side_effect_id, side_effect_name = line.split('\t')

            side_effect_id_list.add(side_effect_id)

    return list(side_effect_id_list)

def write_SIDER_only_graph():
    # Extract graph from raw SIDER2 data
    filename = "../data/meddra_all_label_se.tsv"

    print("Extracting graph from raw data ...")
    G = nx.Graph()
    with open(file=filename, mode='r') as f:
        for line in f:
            _, flat_drug, stereo_drug, _, _, side_effect_id, side_effect_name = line.split('\t')

            flat_drug = flat_drug[:3] + "m" + flat_drug[4:]
            stereo_drug = stereo_drug[:3] + "s" + stereo_drug[4:]

            if flat_drug not in G.nodes():
                G.add_node(flat_drug)
                G.add_node(stereo_drug)

            if side_effect_id not in G.nodes():
                G.add_node(side_effect_id)

            G.add_edge(flat_drug, side_effect_id)
            G.add_edge(stereo_drug, side_effect_id)
    G.remove_node('')
    print("Finished.\n")

    print("Writing SIDER only graph to disc ...")
    filename = "../data/bipartite_SIDER_only_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(G, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_SIDER_only_graph():
    print("Reading SIDER only graph ...\n")
    graph_filename = "../data/bipartite_SIDER_only_graph"
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
    filename = "../data/updated_MedDRA_label_SIDER_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(SIDER_only_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_updated_MedDRA_label_SIDER_graph():
    print("Reading updated MedDRA label SIDER only graph ...\n")
    graph_filename = "../data/updated_MedDRA_label_SIDER_graph"
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
    filename = "../data/jaccard_similarity_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(similarity_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_jaccard_se_similarity_graph():
    print("Reading jaccard side effect graph ...\n")
    graph_filename = "../data/jaccard_similarity_graph"
    with open(graph_filename + '.pkl', 'rb') as f:
        return pickle.load(f)

def write_meddra_graph_to_disc():
    meddra_rdf_graph_filename = "../data/MEDDRA_RDF_original.ttl"
    meddra_graph = rdflib.Graph()
    result = meddra_graph.parse(meddra_rdf_graph_filename, format='n3')
    print(result)

    print("Writing meddra RDF graph to disc ...")
    filename = "../data/meddra_RDF_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(meddra_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def write_enriched_SIDER_graph():

    # read meddra RDF graph from disc
    print("Reading meddra RDF graph ...")
    graph_filename = "../data/meddra_RDF_graph"
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
    target_filename = "../data/MedDRA_enriched_SIDER_RDF_graph.ttl"
    meddra_RDF_graph.serialize(destination=target_filename, format='turtle')
    print("Finished writing ", target_filename, '\n')

    # Write annotation graph to disc
    print("Writing meddra annotation graph to disc ...")
    filename = "../data/SIDER_annotation_graph"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(annotation_graph, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename, '\n')

def get_MedDRA_mapping():
    # Read MRCUI mapping file
    MRCUI_filename = "../data/MRCUI.RRF"
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
    filename = "../data/SIDER_annotation_graph"
    annotation_graph = None
    with open(filename + '.pkl', 'rb') as f:
         annotation_graph = pickle.load(f)
    print("Finished.\n")

    print(annotation_graph.nodes()[:100])

    raise Exception

    drug_list = get_SIDER_drug_list()

    print("Writing annotation file ...")
    annotation_file = "../data/annotation_file_for_groovy.tsv"
    with open(file=annotation_file, mode='w') as f:
        for drug in drug_list:
            if drug not in annotation_graph.nodes():
                continue
            neighbor_list = annotation_graph.neighbors(drug)
            f.write(drug+'\t'+'\t'.join(neighbor_list)+'\n')
    print("Finished.")





if __name__ == '__main__':
    # replace_gap_symbols_in_alignment()
    # test_blast()
    # evaluate_Blast_XML()
    # run_similarity_pipeline(threads=8)

    # test_biopython_PairwiseAligner()

    # write_SIDER_only_graph()
    # write_jaccard_se_similarity_graph()
    # get_jaccard_se_similarity_graph()

    # write_updated_MedDRA_label_SIDER_graph()
    # get_updated_MedDRA_label_SIDER_graph()

    # write_enriched_SIDER_graph()
    write_annotation_file()

