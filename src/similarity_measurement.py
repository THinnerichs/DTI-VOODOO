import numpy as np

import subprocess
from joblib import Parallel, delayed



def test_blast():
    # make_blast_db
    from Bio.Blast.Applications import NcbiblastpCommandline
    database_fasta_file = "../data/fasta_files/CIDm00000043_fasta_800_min_score.fasta"
    query_fasta_file = "../data/fasta_files/CIDm00000043_fasta_800_min_score.fasta"

    # Build blast db
    command = "./makeblastdb -dbtype 'prot' "+\
              "-in " + database_fasta_file
    subprocess.call(command, shell=True)
    print(command)

    raise Exception

    # cline = NcbiblastpCommandline(query=)




if __name__ == '__main__':
    test_blast()