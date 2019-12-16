import numpy as np
from random import shuffle

import subprocess
import queue
import threading
import time
import os
import sys
from joblib import Parallel, delayed
import tqdm

import pickle

from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline, MafftCommandline, MSAProbsCommandline, TCoffeeCommandline

from similarity_measurement import *
import DDI_utils


def write_pruned_SeqIO_fasta_dict():
    # Create protein amino acid sequence from fasta
    print("Reading fasta file to dict...")
    fasta_filename = "../data/STITCH_data/9606.protein.sequences.v11.0.fa"
    protein_aa_seq_dict = SeqIO.index(fasta_filename, 'fasta')
    print("Finished.")

    pruned_dict = {}
    for key, seqio_seq in protein_aa_seq_dict.items():
        split_list = key.split('.')

        organism = split_list.pop(0)
        protein = '.'.join(split_list)
        result = pruned_dict.get(organism, None)

        if not result:
            pruned_dict[organism] = {}

        pruned_dict[organism][protein] = str(seqio_seq.seq)
    print("Finished.")

    print("Writing pruned dict to disc")
    filename = "../data/prot_aa_seq_dict"
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(pruned_dict, f, pickle.HIGHEST_PROTOCOL)
    print("Finished writing ", filename)

def prune_drug_protein_db(min_score=700):

    filename = "../data/STITCH_data/protein_chemical.links.transfer.v5.0.tsv"
    target_filename = "../data/STITCH_data/protein_chemical.links.min_score_" + str(min_score) + ".tsv"

    print("Processing huge file ...")
    with open(file=filename, mode='r') as f, open(file=target_filename, mode='w') as targetfile:
        targetfile.write(f.readline())

        counter = 0

        for line in f:
            counter += 1
            if counter % 10000000 == 0:
                print("Processed lines:", counter)

            if int(line.strip().split('\t')[10]) < min_score:
                continue
            targetfile.write(line)
    print("Finished.")

def separate_prots_to_files(file_min_score=400,
                            min_score=400):
    """
    Iterate over large drug-target-interaction file and filter drugs and their targets to separate files.
    :param file_min_score:              min_score given in filename
    :param min_score:                   min_score that is really applied,
    :return:
    """
    # process drug-protein-interaction file
    print("Processing drug-protein-links data...")
    protein_chemical_links_filename = "../data/STITCH_data/protein_chemical.links.min_score_"+str(file_min_score)+".tsv"

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

            # merge CIDm and s
            drug = split_line[0].replace('s', 'm')
            split_protein = split_line[1].strip().split('.')
            organism = split_protein.pop(0)
            if organism == "9606": # if organism is human
                continue
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
    """
    Create fasta files for all filtered files from separate_prots_to_file() using a previously created dict, mapping
    targets to their amino acid sequence.
    :param min_score:                   min_score both applied and in filename
    :return:
    """
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
            alignment_method='mafft',
            workers=16,
            threads_per_process=2,
            overwrite=False,
            start='',
            end='',
            human_only=False):
    """
    A wrapper for the different MSA techniques that are eventually executed in parallel.

    :param min_score:                       min_score on which pipeline shall be executed
    :param alignment_method:                alignment method which shall be used for the msa, see elif cascade below for options
    :param workers:                         number of workers that each will execute the pipeline, watch memory limitations
    :return:
    """


    fasta_path = '../data/fasta_files/'
    target_path = '../data/alignment_targets/'

    def msa(file):
        # Check whether right min_score is present
        if str(min_score) not in file:
            return

        drug_name = file.split("_")[0].strip()
        # Check whether fasta is existent
        if not os.path.exists(fasta_path+file):
            with open(file="../data/non_existent_fastas", mode='a') as f:
                f.write(drug_name+'\n')
            return
        # Check whether file is empty for speedup
        if os.stat(fasta_path+file).st_size == 0:
            with open(file="../data/empty_fastas", mode='a') as f:
                f.write(drug_name+'\n')
            return

        drug_name = file.split("_")[0].strip()

        # if int(drug_name[4:]) > 1200:
            # return

        target_file = target_path + drug_name + "_"+alignment_method+"_aligned_"+str(min_score)+"_min_score.afa"
        if overwrite or not os.path.exists(target_file):
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
                command = "./mafft-linux64/mafft.bat " \
                          "--anysymbol " \
                          "--auto " + \
                          "--thread " + str(threads_per_process)+" "+\
                          fasta_file + " > " + target_file
                print(command)

            elif alignment_method == 'msaprobs':
                print("MSAProbs not supported yet.")
                raise Exception

                command = MSAProbsCommandline(infile=fasta_file,
                                              outfile=target_file)
                command = "./MSAProbs-0.9.7/MSAProbs/msaprobs " + ' '.join(str(command).split(' ')[1:])
            elif alignment_method == 'kalign':
                command = "./kalign2/kalign -i " + fasta_file + " -o " + target_file
            elif alignment_method == 'famsa':
                command = "./famsa-1.2.5-linux " +\
                          "-t "+str(threads_per_process)+" "+\
                          fasta_file + " " + target_file
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


    # Parallel(n_jobs=50)(delayed(msa)(filename) for filename in os.listdir(fasta_path))

    # for filename in os.listdir(fasta_path):
        # msa(filename)

    q = queue.Queue()

    files = None
    if human_only:
        merged_DDI_graph = DDI_utils.get_merged_DDI_graph()
        drug_list = set(merged_DDI_graph.nodes()) | set(get_SIDER_drug_list())
        # drug_list = get_SIDER_Boyce_Drubank_drug_intersection()
        # drug_list = get_SIDER_drug_list()
        files = [drug_name+"_fasta_" + str(min_score) + "_min_score.fasta" for drug_name in drug_list]
        shuffle(files)
    else:
        files = [file for file in os.listdir(fasta_path) if (str(min_score) in file and os.stat(fasta_path+file).st_size!=0)]

    # shuffle(files)

    if start:
        files = files[int(start):]
    if end:
        files = files[:int(end)-int(start)]

    for fileName in files:
        q.put(fileName)

    def worker():
        while True:
            fileName = q.get()
            if fileName is None:  # EOF?
                return
            msa(fileName)

    threads = [threading.Thread(target=worker) for _i in range(workers)]
    for thread in threads:
        thread.start()
        q.put(None)  # one EOF marker for each thread


def run_hmm_build_pipeline(min_score=700,
                           alignment_method='famsa',
                           workers=20,
                           threads_per_worker=4,
                           rel_weight_method='wpb'
                           ):
    """
    A wrapper for the hmm pipeline that builds the Hidden Markov Model for each multi sequence alignment and also searches
    for it against the whole amount of amino acid sequences.
    :param min_score:                       min_score for both filename
    :param alignment_method:                alignment method that was used in the previous msa step
    :param workers:                         number of workers that will each execute the pipeline for a single alignment
    :param threads_per_worker:              Number of threads used per worker
    :param rel_weight_method:               options: wpb, wgsc, wblosum, wnone; Parameter for hmmbuild
    :return:
    """

    alignment_path = '../data/alignment_targets/'

    def hmm_build_pipeline(file):
        # Check whether right min_score is present
        if str(min_score) not in file:
            return
        # Check whether file is empty for speedup
        if os.stat(alignment_path+file).st_size == 0:
            return


        # hmm build to perform on precomputed fasta files
        hmmbuild_target_path = '../data/hmm_builds/'

        # hmmbuild params
        sym_frac = 0.5
        frag_thresh = 0.5

        drug_name = file.split("_")[0]

        alignment_file = alignment_path + drug_name + "_"+alignment_method+"_aligned_"+str(min_score)+"_min_score.afa"
        hmmbuild_file = hmmbuild_target_path + drug_name + "_" + rel_weight_method + "_"+alignment_method+"_aligned_"+str(min_score)+"_min_score.hmm"

        if os.path.exists(hmmbuild_file) and os.stat(hmmbuild_file).st_size != 0:
            return

        print("Building Hidden Markov Model ...")
        command = "hmmbuild --amino "\
                  "--cpu "+str(threads_per_worker)+" "+\
                  "--" + rel_weight_method + " " +\
                  hmmbuild_file+" "+\
                  alignment_file
        # "--symfrac "+str(sym_frac)+" "+\
        # "--fragthresh " + str(frag_thresh) +" "+\
        # "--"+rel_weight_method+" "+\

        print(command)
        subprocess.call(command, shell=True)
        print("Finished.\n")

    q = queue.Queue()
    files = os.listdir(alignment_path)
    shuffle(files)

    for fileName in files:
        q.put(fileName)

    def worker():
        while True:
            fileName = q.get()
            if fileName is None:  # EOF?
                return
            hmm_build_pipeline(fileName)

    threads = [threading.Thread(target=worker) for _i in range(workers)]
    for thread in threads:
        thread.start()
        q.put(None)  # one EOF marker for each thread

def run_hmm_search_pipeline(min_score=700,
                            alignment_method='famsa',
                            workers=20,
                            threads_per_worker=4,
                            rel_weight_method='wpb'
                            ):
    hmmbuild_target_path = '../data/hmm_builds/'

    def hmm_search_pipeline(file):

        drug_name = file.split("_")[0]
        hmmbuild_file = hmmbuild_target_path + file

        # Check whether right min_score is present
        if str(min_score) not in file:
            return
        # Only calculate builds for certain alignment method
        if alignment_method not in file:
            return
        # Only run if weighting method is correct
        if rel_weight_method not in file:
            return
        # Check whether hmm_build file exists at all
        if os.path.exists(hmmbuild_file) and os.stat(hmmbuild_file).st_size == 0:
            return
        # Check whether file is empty for speedup

        # run hmm search on previously computed hidden markov model over all sequences
        # hmmsearch params
        max_flag = False

        all_prots_fasta_filename = "../data/STITCH_data/9606.protein.sequences.v11.0.fa"

        hmm_search_results_path = '../data/hmm_search_results/'
        hmmsearch_file = hmm_search_results_path + drug_name + "_" + rel_weight_method + "_" + alignment_method + "_aligned_" + str(min_score) + "_min_score.out"

        command = "hmmsearch " + \
                  "--cpu " + str(threads_per_worker) + " " + \
                  "--nonull2 --nobias " +\
                  hmmbuild_file + " " + all_prots_fasta_filename + " > " + hmmsearch_file

        start_time = time.time()
        print("Querying all protein sequences ...")
        print(command)
        subprocess.call(command, shell=True)
        print("Finished in {} seconds.\n".format(time.time() - start_time))

        # Evaluate the results of the search
        predicted_targets_path = "../data/predicted_targets/"
        predicted_targets_file = predicted_targets_path + drug_name + "_" + rel_weight_method + "_" + alignment_method+"_predicted_targets"

        # extract actual predicted targets from hmmsearch output file and write them to a new file
        with open(file=hmmsearch_file, mode='r') as hmmsearch_filehandler, \
                open(file=predicted_targets_file, mode='w') as predicted_targets_filehandler:
            # Skip first 15 lines
            for i in range(15):
                hmmsearch_filehandler.readline()
            for line in hmmsearch_filehandler:
                if line.strip() == "":
                    print("NO MATCHES FOUND!")
                    break
                if "inclusion threshold" in line:
                    break
                split_list = list(filter(None, line.split(' ')))
                e_value = split_list[0]
                protein_id = split_list[8]
                predicted_targets_filehandler.write(protein_id+'\t'+e_value+'\n')

    q = queue.Queue()
    files = os.listdir(hmmbuild_target_path)
    shuffle(files)

    for fileName in files:
        q.put(fileName)

    def worker():
        while True:
            fileName = q.get()
            if fileName is None:  # EOF?
                return
            hmm_search_pipeline(fileName)

    threads = [threading.Thread(target=worker) for _i in range(workers)]
    for thread in threads:
        thread.start()
        q.put(None)  # one EOF marker for each thread

def write_predicted_targets(alignment_method='famsa'):
    hmm_search_path = "../data/hmm_search_results/"
    files = [filename for filename in os.listdir(hmm_search_path)
             if alignment_method in filename]

    for file in files:
        drug_name = file.split('_')[0]
        with open(file=hmm_search_path + file, mode='r') as f:
            for i in range(15):
                f.readline()
            for line in f:
                split_line = list(filter(None, line.split(' ')))  # remove empty strings from list
                e_value = float(split_line[0])









if __name__ == '__main__':
    '''
    separate_prots_to_files(file_min_score=400,
                            min_score=400)
    '''

    # create_fasta_files(min_score=700)

    '''
    args = sys.argv + ['', '']
    start = args[1]
    end = args[2]

    run_MSA(min_score=700,
            alignment_method='mafft',
            workers=2,
            threads_per_process=10,
            start=start,
            end=end,
            human_only=True)
    '''

    # famsa builds
    '''
    run_hmm_build_pipeline(min_score=700,
                           alignment_method='famsa',
                           workers=2,
                           threads_per_worker=10,
                           rel_weight_method='wpb')
    run_hmm_build_pipeline(min_score=700,
                           alignment_method='famsa',
                           workers=2,
                           threads_per_worker=10,
                           rel_weight_method='wgsc')
    run_hmm_build_pipeline(min_score=700,
                           alignment_method='famsa',
                           workers=2,
                           threads_per_worker=10,
                           rel_weight_method='wblosum')
    run_hmm_build_pipeline(min_score=700,
                           alignment_method='famsa',
                           workers=2,
                           threads_per_worker=10,
                           rel_weight_method='wnone')
    '''

    # mafft builds
    run_hmm_build_pipeline(min_score=700,
                           alignment_method='mafft',
                           workers=2,
                           threads_per_worker=10,
                           rel_weight_method='wpb')
    run_hmm_build_pipeline(min_score=700,
                           alignment_method='mafft',
                           workers=2,
                           threads_per_worker=10,
                           rel_weight_method='wgsc')
    run_hmm_build_pipeline(min_score=700,
                           alignment_method='mafft',
                           workers=2,
                           threads_per_worker=10,
                           rel_weight_method='wblosum')
    run_hmm_build_pipeline(min_score=700,
                           alignment_method='mafft',
                           workers=2,
                           threads_per_worker=10,
                           rel_weight_method='wnone')

    '''
    # hmm search
    run_hmm_search_pipeline(min_score=700,
                            alignment_method='famsa',
                            workers=5,
                            threads_per_worker=4,
                            rel_weight_method='wpb')
    run_hmm_search_pipeline(min_score=700,
                            alignment_method='famsa',
                            workers=5,
                            threads_per_worker=4,
                            rel_weight_method='wgsc')
    run_hmm_search_pipeline(min_score=700,
                            alignment_method='famsa',
                            workers=5,
                            threads_per_worker=4,
                            rel_weight_method='wlosum')
    run_hmm_search_pipeline(min_score=700,
                            alignment_method='famsa',
                            workers=5,
                            threads_per_worker=4,
                            rel_weight_method='wnone')
    '''

