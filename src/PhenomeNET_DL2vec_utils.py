import numpy as np

import pickle as pkl
import rdflib
from tqdm import tqdm

import DDI_utils


def parse_drug_indications():
    path = '../data/drug_indications/'
    filename = 'drug_indication_links.csv'

    indication_name_to_UMLS_mapping = {}
    drug_name_list = []
    drug_indication_links = []

    with open(file=path+filename, mode='r') as f:
        for _ in range(3):
            f.readline()

        for line in f:
            line = line.strip()
            term, term_alteration = line.split('\t')[42:44]
            umls_term, umls_phenotype_term, umls_id = line.split('\t')[46:49]

            umls_id = umls_id.strip()

            # Add drug names to drug_list
            drug = line.split('\t')[3]
            drug_synonym = line.split('\t')[15]

            drug_name_list.append(drug)
            drug_name_list.append(drug_synonym)


            for name in [term, term_alteration, umls_term, umls_phenotype_term]:
                if name and umls_id:
                    name = name.strip().replace('"', '').replace("'", "")
                    indication_name_to_UMLS_mapping[name.strip().replace('"','').lower()] = umls_id

            drug_indication_links.append((drug.lower(), umls_id))
    # remove duplicates
    drug_indication_links = list(set(drug_indication_links))

    filename = 'meddra.tsv'
    with open(file=path+filename, mode='r') as f:
        for line in f:
            line = line.strip()
            umls_id, _, _, indication_name = line.split('\t')

            indication_name_to_UMLS_mapping[indication_name.lower()] = umls_id

    # write drug names
    print('Writing drug names to file ...')
    filename = 'drug_names'
    drug_name_list = list(set(drug_name_list))
    with open(file=path+filename, mode='w') as f:
        for drug in drug_name_list:
            if drug:
                f.write(drug+'\n')

    # parse drug name to id mapping
    print('Parsing drug name to id mapping ...')
    filename = 'drug_name_to_ID_mapping'
    drug_name_to_id_mapping = {}
    with open(file=path+filename, mode='r') as f:
        for line in f:
            drug, pubchem_id = line.split('\t')

            pubchem_id = pubchem_id.strip()

            if not pubchem_id or len(pubchem_id) >= 9:
                continue

            pubchem_id = 'CIDm' + (8-len(pubchem_id))*'0' + pubchem_id
            drug_name_to_id_mapping[drug.lower()] = pubchem_id

    mapped_drug_indications = [(drug_name_to_id_mapping[drug], umls_id) for drug, umls_id in drug_indication_links if drug in drug_name_to_id_mapping.keys()]
    print('Num drug indication pairs:', len(mapped_drug_indications))

    # parse drug indication links
    yamanishi_path = '../data/Yamanishi_data/'
    filename = 'disease.txt'

    # parse diseases from yamanishi_dataset
    disease_list = []
    with open(file=yamanishi_path+filename, mode='r') as f:
        for line in f:
            line = line.strip()
            disease_list.append(line.lower().replace('&apos;', ''))

    mapped_disease_list = list(map(lambda d: indication_name_to_UMLS_mapping.get(d, None), disease_list))
    disease_indices = [i for i, disease in enumerate(mapped_disease_list) if disease!=None and disease.startswith('C')]
    print('indications matched:', len(disease_indices), 'of', len(disease_list))

    disease_list = np.array(mapped_disease_list)

    # parse drug_indication matrix provided by https://github.com/luoyunan/DTINet/
    dii_matrix = []
    with open(file=yamanishi_path + 'mat_drug_disease.txt', mode='r') as f:
        for line in f:
            dii_matrix.append(list(map(int, list(line.strip().replace(' ', '')))))
    dii_matrix = np.array(dii_matrix)
    print('dii_matrix.shape', dii_matrix.shape)

    return disease_list[disease_indices], dii_matrix[:, disease_indices], mapped_drug_indications

def write_UMLS_NET_files():
    path_prefix = "../data/drug_indications/"
    onto_prefix = "<http://phenomebrowser.net/umls#{entity}>"

    phenomenet_drugs = get_PhenomeNET_drug_list()

    disease_list, dii_matrix, sharp_drug_indications = parse_drug_indications()
    drug_list = DDI_utils.get_yamanishi_drug_list()

    drug_indices = [list(drug_list).index(drug) for drug in drug_list if drug!=None and not drug.startswith('1')]
    drug_list = drug_list[drug_indices]
    dii_matrix = dii_matrix[drug_indices, :]

    sharp_drug_list = list(set(map(lambda tup: tup[0], sharp_drug_indications)))

    print('drug_list:', drug_list.shape)
    print('pruned dii_list:', dii_matrix.shape)

    print('num drugs:', len((set(sharp_drug_list) | set(drug_list)) & set(phenomenet_drugs)))
    print(set(drug_list) <= set(phenomenet_drugs))

    drug_indication_pairs = []
    for drug_index in range(len(drug_list)):
        interactors = disease_list[dii_matrix[drug_index, :]==1]
        for disease in interactors:
            disease = onto_prefix.format(entity=disease)
            drug_indication_pairs.append((drug_list[drug_index], disease))

    print('drug_indications:', len(drug_indication_pairs) + len(sharp_drug_indications))

    # write association file
    asso_filename = 'drug_indication_association_file'
    print('Writing association file...')
    with open(file=path_prefix+asso_filename, mode='w') as f:
        for drug, term in drug_indication_pairs:
            if drug in phenomenet_drugs:
                f.write(drug+' '+term+'\n')
        for drug, term in sharp_drug_indications:
            if drug in phenomenet_drugs:
                f.write(drug+' '+onto_prefix.format(entity=term)+'\n')

    # write entity list
    ent_filename = 'drug_indication_entity_list'
    print('Writing entitiy list...')
    with open(file=path_prefix + ent_filename, mode='w') as f:
        for drug in set(drug_list) & set(sharp_drug_list):
            if drug in phenomenet_drugs:
                f.write(drug+'\n')

    path_prefix = '../' + path_prefix
    print('Run this in src/DL2vec (with suitable number of workers):')
    command = f"python runDL2vec.py -embedsize 200 -ontology {path_prefix+'UMLS.owl'} -associations {path_prefix+asso_filename} -outfile {path_prefix+'embedding_model'} -entity_list {path_prefix+ent_filename} -num_workers {64}"
    print(command)


def write_PhenomeNET_files(mode='all'):
    path_prefix = "../data/PhenomeNET_data/"
    onto_prefix = "<http://purl.obolibrary.org/obo/{entity}>"


    # parse drug HPO annotations and build updated drug list
    # filename = "drug_SIDER_HPO_annotations.csv"

    drug_list = []
    drug_HPO_pairs = []
    if mode=='drug' or mode=='all':
        UMLS_to_phenomeNET_mapping = DDI_utils.get_UMLS_to_phenomeNET_mapping()
        drug_list, drug_se_pairs = DDI_utils.parse_SIDER()

        drug_HPO_pairs = [(drug, onto_prefix.format(se.replace(':','_'))) for drug, se in drug_se_pairs]

        print('Num drug-HPO-pairs:', len(drug_HPO_pairs))
        print('Num present drugs', len(drug_list))

    drug_HPO_pairs = list(set(drug_HPO_pairs))

    # parse GO and MP annotations for proteins and update protein list
    filename = "final_GO_ProteinID_human.txt"

    protein_GO_list = []
    protein_GO_term_pairs = []
    if mode=='GO' or mode=='all':
        print('Parsing protein-GO associations...')
        with open(file=path_prefix + filename, mode='r') as f:
            for line in f:
                GO_term, protein = line.strip().split('\t')
                GO_term = GO_term.replace(':', '_')
                GO_term = onto_prefix.format(entity=GO_term)

                protein_GO_list.append(protein)
                protein_GO_term_pairs.append((protein, GO_term))

        print('Num protein-GO-associations:', len(protein_GO_term_pairs))

    filename = "final_MP_ProteinID_human.txt"

    protein_MP_list = []
    protein_MP_term_pairs = []
    if mode=='MP' or mode=='all':
        print('Parsing protein-MP associations...')
        with open(file=path_prefix + filename, mode='r') as f:
            for line in f:
                MP_term, protein = line.strip().split('\t')
                MP_term = MP_term.replace(':', '_')
                MP_term = onto_prefix.format(entity=MP_term)

                protein_MP_list.append(protein)
                protein_MP_term_pairs.append((protein, MP_term))
        print('Num protein-MP-associations:', len(protein_MP_term_pairs))

    filename = "final_uberon_ProteinID_human.txt"

    protein_uberon_list = []
    protein_uberon_term_pairs = []
    if mode=='uberon' or mode=='all':
        print('Parsing protein-Uberon associations...')
        with open(file=path_prefix + filename, mode='r') as f:
            for line in f:
                uberon_term, protein = line.strip().split('\t')
                uberon_term = uberon_term.replace(':', '_')
                uberon_term = onto_prefix.format(entity=uberon_term)

                protein_uberon_list.append(protein)
                protein_uberon_term_pairs.append((protein, uberon_term))
        print('Num protein-uberon-associations:', len(protein_uberon_term_pairs))

    # build distinct drugs and proteins lists
    drug_list = sorted(list(set(drug_list)))
    protein_list = None
    if mode=='all':
        protein_list = sorted(list(set(protein_GO_list) & set(protein_MP_list) & set(protein_uberon_list)))
    else:
        protein_list = sorted(list(set(protein_GO_list) | set(protein_MP_list) | set(protein_uberon_list)))


    # write drug and protein lists
    filename = 'PhenomeNET_drug_list'
    if mode=='drug' or mode=='all':
        print('Writing drug and protein list...')
        with open(file=path_prefix+filename+'.pkl', mode='wb') as f:
            pkl.dump(drug_list, f, pkl.HIGHEST_PROTOCOL)

    filename = mode+'_protein_list' if mode in ['GO','MP','uberon','drug'] else 'PhenomeNET_protein_list'
    with open(file=path_prefix + filename + '.pkl', mode='wb') as f:
        pkl.dump(protein_list, f, pkl.HIGHEST_PROTOCOL)

    print('drugs/proteins present:', len(drug_list), len(protein_list))

    # write drug-HPO-annotations and others to annotations file
    filename = mode+'_association_file' if mode in ['GO','MP','uberon','drug'] else 'association_file'
    print('Writing association file...')
    with open(file=path_prefix+filename, mode='w') as f:
        if mode in ['drug','all']:
            for drug, HPO_term in drug_HPO_pairs:
                f.write(drug+' '+HPO_term+'\n')

        if mode in ['GO','all']:
            for protein, GO_term in protein_GO_term_pairs:
                f.write(protein+' '+GO_term+'\n')

        if mode in ['MP','all']:
            for protein, MP_term in protein_MP_term_pairs:
                f.write(protein+' '+MP_term+'\n')

        if mode in ['uberon','all']:
            for protein, uberon_term in protein_uberon_term_pairs:
                f.write(protein+' '+uberon_term+'\n')

    # write entity list
    filename = mode+'_entity_list' if mode in ['GO','MP','uberon','drug'] else 'entity_list'
    print('Writing entitiy list...')
    with open(file=path_prefix + filename, mode='w') as f:
        if mode =='drug' or mode=='all':
            for drug in drug_list:
                f.write(drug+'\n')
        if mode in ['GO','uberon', 'MP']:
            for protein in protein_list:
                f.write(protein+'\n')

def output_example_DL2vec_command(workers=48, embedsize=200):
    print('Running DL2vec embedding generator ...')
    path = '../../data/PhenomeNET_data/'
    onto = path + 'phenomenet.owl'
    asso = path + 'association_file'
    outfile = path + 'embedding_model'
    ents = path + 'entity_list'
    command = 'python runDL2vec.py -embedsize {embedsize} -ontology {onto} -associations {asso} -outfile {outfile} -entity_list {ents} -num_workers {num_workers}'.format(
        embedsize=embedsize,
        onto=onto,
        asso=asso,
        outfile=outfile,
        ents=ents,
        num_workers=workers)

    print('Command:', command)

def get_PhenomeNET_drug_list():
    filename = '../data/PhenomeNET_data/PhenomeNET_drug_list'
    with open(file=filename+'.pkl', mode='rb') as f:
        return pkl.load(f)

def get_PhenomeNET_protein_list(mode='all'):
    path = '../data/PhenomeNET_data/'
    filename = mode+'_protein_list' if mode in ['GO','MP','uberon'] else 'PhenomeNET_protein_list'
    with open(file=path+filename+'.pkl', mode='rb') as f:
        return pkl.load(f)

def write_query_drugpheno_rdf_graph():
    drugpheno_rdf_graph_filename = "../data/PhenomeNET_data/data-2020-12-07/drugphenotype.rdf"
    drugpheno_graph = rdflib.ConjunctiveGraph()
    print('Parsing drugpheno RDF graph...')
    result = drugpheno_graph.parse(drugpheno_rdf_graph_filename, format='xml')

    print('Query drugpheno RDF graph...')
    qres = drugpheno_graph.query(
        """SELECT ?concept ?phenotype
WHERE {
    ?association rdf:type rdf:Statement .
    ?association rdf:predicate obo:RO_0002200 .
    ?association rdf:object ?phenotype .
    ?association rdf:subject ?concept .
    ?concept rdf:type <http://phenomebrowser.net/Drug> .
    ?association dc:provenance ?prov .
    ?prov dc:creator ?creator .
    ?prov dcterms:source ?source .
}
""")

    outfile = '../data/PhenomeNET_data/drugpheno_query_results.tsv'

    print('Writing results...')

    results = []
    counter = 0
    for d,p in qres:
        results.append((d,p))
        counter += 1
        if counter % 1000 == 0:
            print(counter)
    print(len(results))
    with open(file=outfile, mode='w') as f:
        for drug, phenotype in tqdm(results):
            drug = drug.strip().split('/')[-1]
            phenotype = phenotype.strip().split('/')[-1]
            print(drug+'\t'+phenotype, file=f)
    print('Done.', len(results), 'drugpheno pairs.')


if __name__ == '__main__':
    write_PhenomeNET_files('drug')
    # write_query_drugpheno_rdf_graph()
