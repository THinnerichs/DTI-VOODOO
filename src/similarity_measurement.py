import numpy as np

import subprocess
from joblib import Parallel, delayed
import time



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



if __name__ == '__main__':
    # test_blast()
    evaluate_Blast_XML()