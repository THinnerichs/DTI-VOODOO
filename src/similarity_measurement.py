import numpy as np

import subprocess
from joblib import Parallel, delayed
import time
import queue
import threading
import os
from random import shuffle
from tqdm import tqdm

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

    for matrix in ['STR', 'RISLER', 'RAO', 'PAM70', 'PAM30', 'PAM250', 'MDM78', 'MCLACHLAN', 'LEVIN', 'JONES', 'JOHNSON', 'GONNET1992', 'GENETIC',
                   'FENG', 'DAYHOFF', 'BLOSUM90', 'BENNER22', 'BENNER6', 'BENNER74', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80']:

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

        start_time = time.time()

        def help_func(doublet):
            with open(file="../data/score_list", mode='a') as f:
                seq1 = str(doublet[0].seq).replace('-', 'X')
                seq2 = str(doublet[1].seq).replace('-', 'X')
                f.write(str(aligner.score(seq1, seq2))+"\n")
        # help_func = lambda doublet: aligner.score(doublet[0].seq, doublet[1].seq)

        score_list = Parallel(n_jobs=40)(delayed(help_func)(doublet) for doublet in tqdm(itertools.product(database_records, query_records)[:200]))


        # score_list = np.array(score_list)

        # np.save("../data/score_list_nparray.npy", score_list)

    print("This took {} seconds.".format(time.time()-start_time))

def test_blast():

    # parameters
    threads = 8

    # make_blast_db
    fasta_path = "../data/fasta_files/"
    filename = "CIDm00000043_fasta_800_min_score.fasta"
    database_fasta_file = fasta_path + filename
    query_fasta_file = fasta_path + filename

    drug_name = filename.split("_")[0]
    database_name = fasta_path + drug_name + "_blast_db"

    # Build blast db
    print("Building database ...")
    command = "./makeblastdb -dbtype 'prot' "+\
              "-in " + database_fasta_file + " "+\
              "-out " + database_name
    print(command)
    subprocess.call(command, shell=True)
    print("Finished.\n")


    print("Running query ...")
    results_filename = ""+drug_name+"_blast_result.xml"

    blast_command = "./blastp "+\
                    "-task blastp-fast "+\
                    "-num_threads 32 "+\
                    "-query "+query_fasta_file+" "+\
                    "-db "+database_name+" "+\
                    "-out "+results_filename+" "+\
                    "-evalue 1e-20 "+\
                    "-outfmt 5"
    print(blast_command)

    subprocess.call(blast_command, shell=True)
    print("Finished.\n")

def evaluate_Blast_XML():
    drug_name = "CIDm00000043"
    results_filename = ""+drug_name+"_blast_result.xml"

    from Bio.Blast import NCBIXML
    E_VALUE_THRESH = 0.05

    num_lines = None
    fasta_path = "../data/fasta_files/"
    filename = "CIDm00000043_fasta_800_min_score.fasta"
    with open(file=fasta_path+filename, mode='r') as f:
        num_lines = sum(1 for line in f)

    print("Lines:", num_lines)

    raise Exception


    start_time = time.time()
    print("Parsing similarity scores ...")
    prot_prot_sim_dict = {}
    for record in NCBIXML.parse(open(results_filename)):
        if record.alignments: # skip queries with no   matches
            prot_prot_sim_dict[record.query] = {}
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        prot_prot_sim_dict[record.query][alignment.title] = int(hsp.score)

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



if __name__ == '__main__':
    # test_blast()
    # evaluate_Blast_XML()
    # run_similarity_pipeline(threads=8)

    test_biopython_PairwiseAligner()
