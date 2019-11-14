import numpy as np

import subprocess
from joblib import Parallel, delayed



def test_blast():
    # make_blast_db
    from Bio.Blast.Applications import NcbiblastpCommandline

    fasta_path = "../data/fasta_files/"
    filename = "CIDm00000043_fasta_800_min_score.fasta"
    database_fasta_file = fasta_path + filename
    query_fasta_file = fasta_path + filename

    drug_name = filename.split("_")[0]
    database_name = fasta_path + drug_name + "_blast_db"

    # Build blast db
    command = "./makeblastdb -dbtype 'prot' "+\
              "-in " + database_fasta_file + " "+\
              "-out " + database_name
    subprocess.call(command, shell=True)
    print(command)

    

    # cline = NcbiblastpCommandline(query=)




if __name__ == '__main__':
    test_blast()