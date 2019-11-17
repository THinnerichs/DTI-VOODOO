import numpy as np

import subprocess
from joblib import Parallel, delayed



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
    results_filename = ""+drug_name+"_blast_result"

    blast_command = "./blastp "+\
                    "-task blastp-fast "+\
                    "-num_threads 8 "+\
                    "-query "+query_fasta_file+" "+\
                    "-db "+database_name+" "+\
                    "-out "+results_filename+" "+\
                    "-outfmt 0"
    print(blast_command)

    subprocess.call(blast_command, shell=True)
    print("Finished.\n")


if __name__ == '__main__':
    test_blast()