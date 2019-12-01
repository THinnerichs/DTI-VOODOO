import numpy as np

from Bio import SeqIO
import subprocess
import time

import pickle


def write_pruned_SeqIO_fasta_dict():
    # Create protein amino acid sequence from fasta
    print("Reading fasta file to dict...")
    fasta_filename = "../data/protein.sequences.v10.fa"
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
    filename = "../data/protein_chemical.links.transfer.v5.0.tsv"
    target_filename = "../data/protein_chemical.links.min_score_" + str(min_score) + ".tsv"

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


def write_paracetamol_prots_to_file(file_min_score=400,
                                    min_score=700,
                                    mol_name='para'):
    # Create protein amino acid sequence from fasta
    print("Reading pruned dict...")
    dict_filename = "../data/prot_aa_seq_dict"
    protein_aa_seq_dict = None
    with open(dict_filename + '.pkl', 'rb') as f:
        protein_aa_seq_dict = pickle.load(f)
    print("Finished.")

    # process drug-protein-interaction file
    print("Processing drug protein data...")
    protein_chemical_links_filename = "../data/protein_chemical.links.min_score_" + str(file_min_score) + ".tsv"
    # paracetamol_id = 'CIDs00001983'
    rofecoxib_id = 'CIDs00005090'
    counter = 0
    # total_lines = 7019540873
    paracetamol_prots_count = 0
    with open(file=protein_chemical_links_filename, mode='r') as f:
        for line in f:
            counter += 1
            paracetamol_prots_count += 1

            if counter % 5000000 == 0:
                # print("Progress: %.2f %%\r"%(counter*100/total_lines))
                print("Processed lines:", counter)

            if rofecoxib_id not in line:
                continue

            split_line = line.strip().split('\t')
            if int(split_line[10]) < min_score:
                continue

            drug = split_line[0]
            split_protein = split_line[1].strip().split('.')
            organism = split_protein.pop(0)
            protein = ".".join(split_protein)
            score = int(split_line[10])

            aa_seq = protein_aa_seq_dict[organism][protein]

            with open(file="../data/" + mol_name + "_targets", mode='a') as para_handler:
                para_handler.write(organism + '.' + protein + "\t" + aa_seq + '\t' + str(score) + "\n")
    print("Finished.")

    print("Total line count:", counter)
    print("Paracetamol count:", paracetamol_prots_count)


def write_paracetamol_prots_to_fasta(min_score=400,
                                     mol_name='para'):
    para_filename = "../data/" + mol_name + "_targets"
    para_fasta_filename = "../data/" + mol_name + "_fasta_" + str(min_score) + "_min_score.fasta"
    print("Processing {} and writing {} ...".format(para_filename, para_fasta_filename))
    with open(file=para_filename, mode='r') as para_file, open(file=para_fasta_filename, mode='w') as fasta_file:
        for line in para_file:
            protein, aa_seq, score = line.split('\t')
            if int(score) < min_score:
                continue

            fasta_file.write(">" + protein + '\n')
            fasta_file.write(aa_seq + '\n')

    print("Finished.")


def run_stitch_db_query():
    import requests
    import sys

    # build url
    paracetamol_id = "CID000001983"
    # paracetamol_id = "CID3676"
    required_score = 400
    url = "http://stitch.embl.de/api/tsv/interactors?identifier=" + paracetamol_id + \
          "&required_score=" + str(required_score) + "%0A" \
                                                     "species=9606"
    url2 = "http://stitch.embl.de/api/tsv/resolve?identifier=" + paracetamol_id
    print(url2)

    # Run query with exception handling
    r = requests.get(url)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    seq = r.text

    print(seq)


def eliminate_para_target_duplicates():
    filename = "../data/para_targets"
    protein_set = set()
    with open(file=filename, mode='r') as f:
        for line in f:
            org_prot, _, score = line.split('\t')
            if int(score) < 700:
                continue
            split_prot = org_prot.split('.')
            organism = split_prot.pop(0)
            protein = ".".join(split_prot)

            protein_set.add(protein)

    print("set size:", len(protein_set))

    with open(file="../data/setified_para_proteins_700", mode="w") as f:
        for ele in protein_set:
            f.write(ele + '\n')


def run_para_multi_sequence_alignment(min_score=700,
                                      alignment_method='clustalo',
                                      mol_name='para'):
    from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline, MafftCommandline, \
        MSAProbsCommandline, TCoffeeCommandline

    in_file = "../data/" + mol_name + "_fasta_" + str(min_score) + "_min_score.fasta"
    out_file = "../data/" + mol_name + "_" + alignment_method + "_aligned_" + str(min_score) + "_min_score.afa"
    start_time = time.time()

    command = None
    if alignment_method == 'clustalo':
        command = ClustalOmegaCommandline(infile=in_file,
                                          outfile=out_file,
                                          verbose=True,
                                          auto=True)
        command = "./" + str(command)
    elif alignment_method == 'muscle':
        command = MuscleCommandline(input=in_file,
                                    out=out_file)
        command = "./" + str(command)
    elif alignment_method == 'mafft':
        command = MafftCommandline(input=in_file)
        command = "./mafft-linux64/mafft.bat --anysymbol --auto " + ' '.join(
            str(command).split(' ')[1:]) + " > " + out_file

        print(command)

        print("Starting {} alignment ...".format(alignment_method))
        subprocess.call(str(command), shell=True)

        print("Finished in {} sec.".format(time.time() - start_time))

    elif alignment_method == 'msaprobs':
        command = MSAProbsCommandline(infile=in_file,
                                      outfile=out_file)
        command = "./MSAProbs-0.9.7/MSAProbs/msaprobs " + ' '.join(str(command).split(' ')[1:])
    elif alignment_method == 'kalign':
        command = "./kalign2/kalign -i " + in_file + " -o " + out_file
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
    print("Finished in {} sec.".format(time.time() - start_time))


def para_PWM_from_alignment(min_score=700,
                            alignment_method='mafft',
                            mol_name='para'):
    from Bio.Alphabet import IUPAC, Gapped
    from Bio import AlignIO, Alphabet, motifs, SeqIO

    alphabet = Gapped(IUPAC.extended_protein)

    filename = "../data/" + mol_name + "_" + alignment_method + "_aligned_" + str(min_score) + "_min_score.afa"

    print("Reading alignment file ...")
    alignment = AlignIO.read(filename,
                             'fasta',
                             alphabet=alphabet)
    print("Finished.")

    import math
    in_file = "../data/" + mol_name + "_fasta_" + str(min_score) + "_min_score.fasta"
    min_seq_length = math.inf
    for record in SeqIO.parse(in_file, 'fasta'):
        if min_seq_length > len(str(record.seq)):
            min_seq_length = len(str(record.seq))

    print(min_seq_length)

    raise Exception

    print("Creating motifs ...")
    m = motifs.create([x.seq for x in alignment], alphabet=alphabet)
    print("Finished.")

    # m.weblogo(fname="../results/weblogo_"+mol_name+"_"+alignment_method+"_"+str(min_score)+"min_score.png")

    # print(m.consensus)
    # print(m.counts)

    # See http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc213

    # Usually, pseudocounts are added to each position before normalizing. This avoids overfitting of the position-weight
    # matrix to the limited number of motif instances in the alignment, and can also prevent probabilities from becoming zero.
    pwm = m.counts.normalize(pseudocounts=0.5)
    # print(pwm)

    pssm = pwm.log_odds()

    pssm_string = pssm.__str__()
    pssm_lines = pssm_string.split('\n')

    print(pssm_lines[1])

    # pssm_mat = [[word.strip() for word in line.split('\t')] for line in pssm_string.split('\n')]

    raise Exception

    pssm = np.array(pssm.__str__())

    # print(pssm)
    print(pssm.shape)

    raise Exception

    print("Calculating distribution ...")
    all_fastas_filename = "../data/protein.sequences.v10.fa"
    # distribution = pssm.distribution(precision=10 ** 4)
    print("Finished.")

    # Obtain threshold from distribution
    # We can specify the requested false-positive rate (probability of “finding” a motif instance in background generated sequence)
    # threshold = distribution.threshold_fpr(0.01)
    # or the false-negative rate (probability of “not finding” an instance generated from the motif)
    # threshold = distribution.threshold_fnr(0.1)
    # or a threshold (approximately) satisfying some relation between the false-positive rate and the false-negative rate (fnr/fpr = t)
    # threshold = distribution.threshold_balanced(1000)
    # or a threshold satisfying (roughly) the equality between the −log of the false-positive rate and the information content (as used in patser software by Hertz and Stormo)
    # threshold = distribution.threshold_patser()

    threshold = 5

    print("Analyzing fasta all sequences with distribution ...")
    counter = 0
    hit_counter = 0
    for record in SeqIO.parse(all_fastas_filename, 'fasta'):
        # Set threshold constant to 3.0, or according to distribution
        print("COUNTER:", counter)
        for position, score in pssm.search(record.seq, threshold=threshold, both=False):
            print("hit_counter: {},\tposition: {},\tscore: {}".format(hit_counter, position, score))
            hit_counter += 1

        counter += 1
        if hit_counter == 100:
            break
    print("Finished.")


def run_HMMER_build(min_score=700,
                    alignment_method='mafft',
                    mol_name='para',
                    sym_frac=0.5,
                    frag_thresh=0.5,
                    rel_weight_method='wpb',
                    cores=2):
    """
    :param min_score:
    :param alignment_method:
    :param mol_name:
    :param sym_frac:
    :param frag_thresh:
    :param rel_weight_method:           options: wpb, wgsc, wblosum
    :return:
    """
    alignment_file = "../data/" + mol_name + "_" + alignment_method + "_aligned_" + str(min_score) + "_min_score.afa"
    out_file = "../data/hmm_build_" + mol_name + "_" + alignment_method + "_aligned_" + str(
        min_score) + "_min_score.hmm"

    print("Building Hidden Markov Model ...")
    command = "hmmbuild --amino " \
              "--cpu " + str(cores) + " " + \
              "--symfrac " + str(sym_frac) + " " + \
              "--fragthresh " + str(frag_thresh) + " " + \
              "--" + rel_weight_method + " " + \
              out_file + " " + \
              alignment_file
    print(command)
    subprocess.call(command, shell=True)
    print("Finished.\n")


def run_HMMER_search(min_score=700,
                     alignment_method='mafft',
                     mol_name='para',
                     max_flag=False,
                     cores=2):
    hmm_file = "../data/hmm_build_" + mol_name + "_" + alignment_method + "_aligned_" + str(
        min_score) + "_min_score.hmm"
    all_prots_fasta_filename = "../data/protein.sequences.v10.fa"
    out_file = "../data/hmm_search_" + mol_name + "_" + alignment_method + "_aligned_" + str(
        min_score) + "_min_score.out"

    command = "hmmsearch " + \
              ("--max " if max_flag else "") + \
              "--cpu " + str(cores) + " " + \
              hmm_file + " " + all_prots_fasta_filename + " > " + out_file

    start_time = time.time()
    print("Querying all protein sequences ...")
    print(command)
    subprocess.call(command, shell=True)
    print("Finished in {} seconds.\n".format(time.time() - start_time))


def create_HMMER_logo(min_score=700,
                      alignment_method='mafft',
                      mol_name='para'):
    hmm_file = "../data/hmm_build_" + mol_name + "_" + alignment_method + "_aligned_" + str(
        min_score) + "_min_score.hmm"
    logo_file = "../results/hmm_logo_" + mol_name + "_" + alignment_method + "_aligned_" + str(min_score) + "_min_score"

    print("Building HMM Logo ...")
    command = "hmmlogo " + hmm_file + " > " + logo_file
    print("Finished.")


def evaluate_HMMER_search(min_score=700,
                          alignment_method='mafft',
                          mol_name='para'):
    out_file = "../data/hmm_search_" + mol_name + "_" + alignment_method + "_aligned_" + str(
        min_score) + "_min_score.out"
    print("Evaluating {}".format(out_file))

    protein_id_list = []
    with open(file=out_file, mode='r') as f:
        # Skip first 14 lines
        for i in range(14):
            f.readline()
        for line in f:
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
            protein_id_list.append(protein_id)

    print("Contains {} proteins below threshold.".format(len(protein_id_list)))

    # Get all proteins from mol_name targets
    targets = []
    with open(file="../data/" + mol_name + "_targets") as f:
        for line in f:
            protein = line.split('\t')[0].strip()
            targets.append(protein)

    intersection = set(protein_id_list) & set(targets)
    others = set(protein_id_list) - intersection

    print("Detected {} of {}".format(len(intersection), len(targets)))
    print("And {} others\n".format(len(others)))


from joblib import Parallel, delayed
import queue
import os
import itertools
import threading


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

    for matrix in ['STR', 'RISLER', 'RAO', 'PAM70', 'PAM30', 'PAM250', 'MDM78', 'MCLACHLAN', 'LEVIN', 'JONES',
                   'JOHNSON', 'GONNET1992', 'GENETIC',
                   'FENG', 'DAYHOFF', 'BLOSUM90', 'BENNER22', 'BENNER6', 'BENNER74', 'BLOSUM45', 'BLOSUM50', 'BLOSUM62',
                   'BLOSUM80']:
        print("\n" + matrix)
        aligner.substitution_matrix = substitution_matrices.load(matrix)

        with open(file="../data/score_list", mode='a') as f:
            f.write("\nMATRIX: " + matrix + "\n")

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
                f.write(str(aligner.score(seq1, seq2)) + "\n")

        # help_func = lambda doublet: aligner.score(doublet[0].seq, doublet[1].seq)

        score_list = Parallel(n_jobs=40)(delayed(help_func)(doublet) for doublet in
                                         tqdm(list(itertools.product(database_records[:20], query_records[:20]))))

        # score_list = np.array(score_list)

        # np.save("../data/score_list_nparray.npy", score_list)

    print("This took {} seconds.".format(time.time() - start_time))


def replace_gap_symbols_in_alignment(samples=500):
    from Bio import SeqIO

    alignment_path = "../data/alignment_targets/"
    database_filename = "CIDm00000003_kalign_aligned_800_min_score.afa"
    query_filename = "CIDm00000006_kalign_aligned_800_min_score.afa"

    with open(file=alignment_path + database_filename[:-3] + "fasta", mode='w') as f:
        for record in list(SeqIO.parse(alignment_path + database_filename, 'fasta'))[:samples]:
            f.write(">" + str(record.id) + "\n")
            f.write(str(record.seq).replace('-', 'X') + "\n")

    with open(file=alignment_path + query_filename[:-3] + "fasta", mode='w') as f:
        for record in list(SeqIO.parse(alignment_path + query_filename, 'fasta'))[:samples]:
            f.write(">" + str(record.id) + "\n")
            f.write(str(record.seq).replace('-', 'X') + '\n')


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
    command = "./makeblastdb -dbtype 'prot' " + \
              "-input_type 'fasta' " + \
              "-in " + alignment_path + database_filename + " " + \
              "-out " + database_name
    print(command)
    subprocess.call(command, shell=True)
    print("Finished.\n")

    print("Running query ...")
    results_filename = "../data/" + drug_name + "_blast_result.xml"

    blast_command = "./blastp " + \
                    "-task blastp-fast " + \
                    "-num_threads " + str(threads) + " " + \
                    "-query " + alignment_path + query_filename + " " + \
                    "-db " + database_name + " " + \
                    "-out " + results_filename + " " + \
                    "-outfmt 5"
    # "-evalue 0.05 "+\
    print(blast_command)

    subprocess.call(blast_command, shell=True)
    print("Finished.\n")


def evaluate_Blast_XML():
    drug_name = "CIDm00000003"
    results_filename = "../data/" + drug_name + "_blast_result.xml"

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
        if record.alignments:  # skip queries with no   matches
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
                        "-evalue " + str(e_value_threshold) + " " \
                                                              "-outfmt 5"
        print(blast_command)

        subprocess.call(blast_command, shell=True)
        print("Finished.\n")

        # Evaluate query
        from Bio.Blast import NCBIXML

        # Count lines in both files for normalization of results
        file1_num_lines = None
        database_fasta_file = "../data/fasta_files/" + database_drug_name + "_fasta_" + str(
            min_score) + "_min_score.afa"
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
        final_score = score_sum / (file1_num_lines / 2 + file2_num_lines / 2)

        evaluation_results_filename = "../data/similarity_results"
        with open(file=evaluation_results_filename, mode='a') as f:
            f.write(query_drug_name + "\t" + database_drug_name + "\t" + str(final_score))

    q = queue.Queue()

    help_file_doublets = itertools.product(os.listdir(alignment_path), os.listdir(alignment_path))

    # omit symmetric doublets
    file_doublets = []
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
    num_threads = [threading.Thread(target=worker) for _i in range(int(max_overall_threads / threads))]
    for thread in num_threads:
        thread.start()
        q.put(None)  # one EOF marker for each thread


if __name__ == '__main__':
    # write_pruned_SeqIO_fasta_dict()
    # write_paracetamol_prots_to_file()

    # eliminate_para_target_duplicates()
    # run_stitch_db_query()
    # prune_drug_protein_db(min_score=400)

    # write_paracetamol_prots_to_file(file_min_score=400, min_score=700)

    # write_paracetamol_prots_to_fasta(min_score=700)

    '''
    for mn in ['para', 'rofec']:
        for ms in [700, 800]:
            for am in ['mafft', 'kalign']:
                run_para_multi_sequence_alignment(min_score=700,
                                                  mol_name='para',
                                                  alignment_method='mafft')
    '''

    '''
    for mn in ['para', 'rofec']:
        run_para_multi_sequence_alignment(min_score=800,
                                          mol_name=mn,
                                          alignment_method='clustalo')
    '''

    '''
    para_PWM_from_alignment(min_score=800,
                            alignment_method='mafft',
                            mol_name='para')
    '''

    '''
    for m_s in [700, 800]:
        for a_m in ['mafft', 'kalign']:
            for mol in ['para', 'rofec']:
                # run_HMMER_build(min_score=m_s, alignment_method=a_m, mol_name=mol)
                # run_HMMER_search(min_score=m_s, alignment_method=a_m, mol_name=mol)
                evaluate_HMMER_search(min_score=m_s, alignment_method=a_m, mol_name=mol)
    '''

    for sym_frac in [i / 10.0 for i in range(3, 6)]:
        for frag_thresh in [i / 10.0 for i in range(3, 6)]:
            print("SYM_FRAC", sym_frac, "FRAG_THRESH", frag_thresh)
            run_HMMER_build(min_score=700, alignment_method='mafft', mol_name='para', sym_frac=sym_frac,
                            frag_thresh=frag_thresh, cores=16)
            run_HMMER_search(min_score=700, alignment_method='mafft', mol_name='para', cores=16)
            evaluate_HMMER_search(min_score=700, alignment_method='mafft', mol_name='para')