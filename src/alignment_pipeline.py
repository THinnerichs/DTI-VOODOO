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
        if not('m' in files[i] and os.path.exists(path+'/'+files[i].replace('m','s'))):
            continue

        m_file = files[i]
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

        print("intersection:", len(set(m_targets) & set(s_targets)), '\n')

def create_fasta_files(min_score=800):

    # Create protein amino acid sequence from fasta
    print("Reading pruned dict...")
    dict_filename = "../data/prot_aa_seq_dict"
    protein_aa_seq_dict = None
    with open(dict_filename + '.pkl', 'rb') as f:
        protein_aa_seq_dict = pickle.load(f)
    print("Finished.")

    path = "../data/drug_target_relations/"
    files = os.listdir(path)

    for i in range(len(files)):
        if i % 2000 == 0:
            print("Processed %.2f %% of files"%(i*100/len(files)))
        file = files[i]
        drug_name = file.strip()[:-8]
        fasta_filename = "../data/fasta_files/"+drug_name+"_fasta_" + str(min_score) + "_min_score.fasta"

        with open(file=path+file, mode='r') as filehandler, open(file=fasta_filename, mode='w') as fasta_filehandler:
            for line in filehandler:
                protein, score = line.split('\t')
                if int(score) < min_score:
                    continue

                split_protein = protein.strip().split('.')
                organism = split_protein.pop(0)
                protein = ".".join(split_protein)

                aa_seq = protein_aa_seq_dict[organism][protein]

                fasta_filehandler.write(">"+organism+'.'+protein+'\n')
                fasta_filehandler.write(aa_seq +'\n')

def run_MSA(min_score=800,
            alignment_method='mafft'):


    fasta_path = '../data/fasta_files/'
    target_path = '../data/alignment_targets/'

    def msa(file):
        # Check whether right min_score is present
        if str(min_score) not in file:
            return
        # Check whether file is empty for speedup
        if os.stat(fasta_path+file).st_size == 0:
            return

        drug_name = file.split("_")[0].strip()

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


    Parallel(n_jobs=50)(delayed(msa)(filename) for filename in os.listdir(fasta_path))

def write_predicted_targets(min_score=800,
                            alignment_method='kalign'):

    alignment_path = '../data/alignment_targets/'

    def hmm_pipeline(file):

        # hmm build to perform on precomputed fasta files
        hmmbuild_target_path = '../data/hmm_builds/'

        # hmmbuild params
        sym_frac = 0.5
        frag_thresh = 0.5
        rel_weight_method = 'wpb'
        cores = 2

        drug_name = file.split("_")[0]

        alignment_file = alignment_path + drug_name + "_"+alignment_method+"_aligned_"+str(min_score)+"_min_score.afa"
        hmmbuild_file = hmmbuild_target_path + drug_name + "_"+alignment_method+"_aligned_"+str(min_score)+"_min_score.hmm"

        print("Building Hidden Markov Model ...")
        command = "hmmbuild --amino "\
                  "--cpu "+str(cores)+" "+\
                  "--symfrac "+str(sym_frac)+" "+\
                  "--fragthresh " + str(frag_thresh) +" "+\
                  "--"+rel_weight_method+" "+\
                  hmmbuild_file+" "+\
                  alignment_file
        print(command)
        subprocess.call(command, shell=True)
        print("Finished.\n")

        # Parallel(n_jobs=10)(delayed(hmm_build)(filename) for filename in [file for file in os.listdir(alignment_path) if file.startswith("CID")])

        # run hmm search on previously computed hidden markov model over all sequences
        # hmmsearch params
        max_flag = False
        cores = 1

        all_prots_fasta_filename = "../data/protein.sequences.v10.fa"

        hmm_search_results_path = '../data/hmm_search_results/'
        hmmsearch_file = hmm_search_results_path + drug_name + "_" + alignment_method + "_aligned_" + str(min_score) + "_min_score.out"

        command = "hmmsearch " + \
                  ("--max " if max_flag else "") + \
                  "--nonli " + \
                  "--nontextw "+\
                  "--cpu " + str(cores) + " " + \
                  hmmbuild_file + " " + all_prots_fasta_filename + " > " + hmmsearch_file

        start_time = time.time()
        print("Querying all protein sequences ...")
        print(command)
        subprocess.call(command, shell=True)
        print("Finished in {} seconds.\n".format(time.time() - start_time))

        # Evaluate the results of the search
        predicted_targets_path = "../data/predicted_targets/"
        predicted_targets_file = predicted_targets_path + drug_name + "_predicted_targets"

        # extract actual predicted targets from hmmsearch output file and write them to a new file
        with open(file=hmmsearch_file, mode='r') as hmmsearch_filehandler, open(file=predicted_targets_file, mode='w') as predicted_targets_filehandler:
            # Skip first 14 lines
            for i in range(14):
                hmmsearch_filehandler.readline()
            for line in hmmsearch_filehandler:
                if line.strip() == "":
                    print("NO MATCHES FOUND!")
                    break
                if "inclusion threshold" in line:
                    break
                split_list = []
                for ele in line.split(' '):
                    if not ele == '':
                        split_list.append(ele)
                protein_id = split_list[8].strip()
                predicted_targets_filehandler.write(protein_id+'\n')

    Parallel(n_jobs=8)(delayed(hmm_pipeline)(filename) for filename in os.listdir(alignment_path))




if __name__ == '__main__':
    '''
    separate_prots_to_files(file_min_score=400,
                            min_score=400)
    '''

    # merge_drug_files()

    create_fasta_files(min_score=700)

    run_MSA(min_score=800,
            alignment_method='kalign')



