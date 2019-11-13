import numpy as np

from Bio import SeqIO
import subprocess
import time
import os
from joblib import Parallel, delayed

import pickle

from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline, MafftCommandline, MSAProbsCommandline, TCoffeeCommandline

def separate_prots_to_files(file_min_score=400,
                            min_score=400):

    # process drug-protein-interaction file
    print("Processing drug protein data...")
    protein_chemical_links_filename = "../data/protein_chemical.links.min_score_"+str(file_min_score)+".tsv"


    counter = 0
    with open(file=protein_chemical_links_filename, mode='r') as f:
        f.readline()
        for line in f:
            counter += 1

            if counter%5000000==0:
                # print("Progress: %.2f %%\r"%(counter*100/total_lines))
                print("Processed lines:", counter)

            split_line = line.strip().split('\t')
            if int(split_line[10]) < min_score:
                continue

            drug = split_line[0]
            split_protein = split_line[1].strip().split('.')
            organism = split_protein.pop(0)
            protein = ".".join(split_protein)
            score = int(split_line[10])

            with open(file="../data/drug_target_relations/"+drug+"_targets", mode='a') as drug_handler:
                drug_handler.write(organism+'.'+protein + '\t' + str(score) + '\n')
    print("Finished.")

def merge_drug_files():
    path = '../data/drug_target_relations/'
    files = os.listdir(path)

    m_file = None
    for i in range(100000000):
        if 'm' in files[i]:
            m_file = files[i]
            break

    print("m_file:", m_file)

    m_targets = []
    with open(path+m_file, mode='r') as f:
        for line in f:
            m_targets.append(line.split('\t')[0])

    s_targets = []
    with open(path+'/'+m_file.replace('m','s'), mode='r') as f:
        for line in f:
            s_targets.append(line.split('\t')[0])

    print("length m_targets:", len(m_targets))
    print("length s_targets:", len(s_targets))

    print("intersection:", len(set(m_targets) & set(s_targets)))

def create_fasta_files(min_score=800):
    path = "../data/drug_target_relations/"
    files = os.listdir(path)

    for file in files:
        drug_name = file.strip()[:-8]
        fasta_filename = "../data/fasta_files/"+drug_name+"_fasta_" + str(min_score) + "_min_score.fasta"

        with open(file=path+file, mode='r') as filehandler, open(file=fasta_filename, mode='w') as fasta_filehandler:
            for line in filehandler:
                protein, aa_seq, score = line.split('\t')
                if int(score) < min_score:
                    continue

                fasta_filehandler.write(">"+protein+'\n')
                fasta_filehandler.write(aa_seq+'\n')

def run_MSA(min_score=800,
            alignment_method='mafft'):


    fasta_path = '../data/fasta_files/'
    target_path = '../data/alignment_targets/'

    def msa(file):
        drug_name = file.strip()[:-8]
        target_file = target_path + drug_name + "_"+alignment_method+"_aligned_"+str(min_score)+"_min_score.afa"
        if not os.path.exists(target_file):
            fasta_file = "../data/fasta_files/"+drug_name+"_fasta_" + str(min_score) + "_min_score.fasta"

            start_time = time.time()

            command = None
            if alignment_method == 'clustalo':
                # literally takes forever and thus disqualifies itself for this usecase
                command = ClustalOmegaCommandline(infile=fasta_file,
                                                  outfile=target_file,
                                                  verbose=True,
                                                  auto=True)
                command = "./" + str(command)
            elif alignment_method == 'muscle':
                command = MuscleCommandline(input=fasta_file,
                                            out=target_file)
                command = "./" + str(command)
            elif alignment_method == 'mafft':
                command = MafftCommandline(input=fasta_file)
                command = "./mafft-linux64/mafft.bat --anysymbol --auto " + ' '.join(str(command).split(' ')[1:]) + " > " + target_file

                print(command)

                print("Starting {} alignment ...".format(alignment_method))
                subprocess.call(str(command), shell=True)

                print("Finished in {} sec.".format(time.time()-start_time))

            elif alignment_method == 'msaprobs':
                print("MSAProbs not supported yet.")
                raise Exception

                command = MSAProbsCommandline(infile=fasta_file,
                                              outfile=target_file)
                command = "./MSAProbs-0.9.7/MSAProbs/msaprobs " + ' '.join(str(command).split(' ')[1:])
            elif alignment_method == 'kalign':
                command = "./kalign2/kalign -i " + fasta_file + " -o " + target_file
            else:
                print("No valid alignment method selected.")
                raise Exception

            '''
            elif alignment_method == 'tcoffee':
                command = TCoffeeCommandline(infile=in_file,
                                             outfile=out_file,
                                             verbose=True,
                                             auto=True)
            '''

            print(command)

            print("Starting {} alignment ...".format(alignment_method))
            subprocess.call(str(command), shell=True)
            print("Finished in {} sec.\n".format(time.time()-start_time))

            return drug_name


    Parallel(n_jobs=16)(delayed(msa)(filename) for filename in os.listdir(fasta_path))

def write_predicted_targets(min_score=800,
                            alignment_method='mafft',
                            sym_frac=0.5,
                            frag_thresh=0.5,
                            rel_weight_method='wpb',
                            cores=2):

    fasta_path = '../data/fasta_files/'
    target_path = '../data/alignment_targets/'

    def hmm_build(file):
        drug_name = file

        fasta_file = "../data/fasta_files/"+drug_name+"_fasta_" + str(min_score) + "_min_score.fasta"
        target_file = target_path + drug_name + "_"+alignment_method+"_aligned_"+str(min_score)+"_min_score.afa"

        print("Building Hidden Markov Model ...")
        command = "hmmbuild --amino "\
                  "--cpu "+str(cores)+" "+\
                  "--symfrac "+str(sym_frac)+" "+\
                  "--fragthresh " + str(frag_thresh) +" "+\
                  "--"+rel_weight_method+" "+\
                  out_file+" "+\
                  alignment_file
        print(command)
        subprocess.call(command, shell=True)
        print("Finished.\n")




if __name__ == '__main__':
    separate_prots_to_files(file_min_score=400,
                            min_score=400)
