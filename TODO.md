## TODO

Also write down failed approach?

Options:
- GCN vs GAT vs SAGE vs others
- Protein vs drug split

- molecular drug representation:
  - SMILES transformer?
  - proper GCN over molecular structure? (What node features? See other papers for that)
  - take model from them and enhance it with PPI graph?
- Model for drug ? For each protein?


Draw images for Networks
- PPI & DDI & Side effects data
- Molecular predictor
- Graph networks
- Images for splits

Robert talk:
  - 10 years ago side effect similarity paper in Science  (0.7 AUC)
  - locality in PPI network coincides with phenotypes
    - See paper in Skype

  - function phenotypes

- What is usage of/hypothesis for 
  - GO / protein functions
  - anatomical site of expression
  - loss of function phenotypes
    - [link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933875/)

- remove self loops? 
- Resnik Similarity drug-protein-similarity

- Remove self-loops from GCN


Data for experimental only:
- 300: num_nodes 1828, num_edges 4494
- 500: num_nodes 1725, num_edges 3620
- 700: num_nodes 1578, num_edges 3043

Robert talk 28.05.:
- side effects and phenotype similarity
  - see mouse model features 
  - use Siamese network

- semsim out of order? values at identity are not 1 
  - Is it algorithm or my dumbness?
  - Get maximum of each row and its index

Data for PPI graph reduction:
- no reduction: nodes: 16042, edges: 401752
- reduction: 
  - 100: nodes: 17982, edges: 1711038 | 17978/1705657
  - 300: nodes: 15647, edges: 220167 | 15616/218487
  - 500: nodes: 12120, edges: 80325 | 11377/71528
  - 700: nodes: 8290, edges: 48558 | 7963/45906
  - 900: nodes: 5055, edges: 23601 | 3396/14859

Robert 04.06.:
- side effect and protein function similarity HPO (Human phenotype ontology)
- Interpro motive 

- phenotypes and functions localize on PPI network

- take example: 
  - e.g. Diclophenac 

------------------------------

1. map SIDER side effects to HPO
2. See what DL2vec has actually done
  2.1. Ontology comes in in which form?
  2.2. How to build GO function according to this
  2.3. How to build HPO according to this
  2.4. How to build the actual model?!
3. get GO functions for proteins
4. Build ontologies according to DL2vec
5. Rebuild their model

------------------------------

- Rebuild prot to gene mapping with UniProt mapping stuff

- add Mouse phenotypes to ontology 
  - rebuild association file with updated mapping

- Build plain (non-graph conv.) network (siamese) for similarity

- take example and test hypothesis for Diclophenac
  - DDI:
    - take DDIs and analyze their targets within PPI graph
      - does it localize? how to measure? 
        - #edges between them vs avg edge count?
      - anything obvious in there?
  - side effects:
    - does side effect similarity work? does it localize? 
      - take highest values from similarity matrix and check
  - DL2vec:
    - 

- PPI not only experimental:
  - also take:
    - from other databases

- run model with updated walks

- Update side effect similarity pipeline with proper mapping 
  - and obviously also calculate it 

- Check what part of DL2vec (@Jun) yielded that high performance
  - Just check the different embeddings obtained from different embeddings 
- 

-----
What I have done since last meeting 25.06.:
- DL2vec walks were stuck in 'HasAssociation'
- non-graph nn over DL2vec -> no significant difference
- fixed some errors in DL2vec walks
- Manual check
  - 1.4% in DDI-targest vs. 0.15% for rof
  - 1.7% for paracetamol
  - 5.9% for dapagliflozin
  BUT those edges aren't specific at all 16000 overlapping edges, mainly due to protein hubs
  -> same for side effect similarity, BUT has much higher specificity

- diameter of graph:
```
len([len(c) for c in sorted(nx.connected_components(ppi_graph), key=len, reverse=True)]) = 134
```
  


-----
Robert Talk 25.06.:
- remove normalization from GCN layer
1. Does network for my 3 drugs work?
2. train classification in new manner
  - swap negative samples or train only on targets
3. add GO functions, phenotypes, diseases as node features

------

Plan:
1. Search for degree aware GCN methods
2. if not possible build one of the following methods:
  2.1. decentile (!)
  2.2. full encoding
  2.3. normalized number
3. run stuff

Network:
- build embedding from proteins and combine with drugs afterwards
  - what to use as node features?
    - DL2vec?
- retest removal of normalization

DL2Vec TODO:
- test quality of my DL2vec features
  - use my drug features as features for prot\_func
PhenomeNET:
- SIDER included?
Training: 
- put mask over results and randomly select negativ samples

What I have done since last meeting:
- rerun without HasAssociation limitation
- test quality of my DL2vec features
  - use my drug features as features for prot\_func


- Test pathways manually
- cosine similarity of target and drug as single features 
- create protein embeddings over GO, UBERON and HP
  
- apply weighting with respect to node degree?

- Just test on the their non-PPI data -> too few proteins in non-PPI data
- check graph data building whether its correct

- protein union instead of intersection

- predict PPIs by DL2vec embeddings (GO) as sanity test for embeddings

## A new hope (27.10.2020)
- test GCN with different approaches:
  - pre similarity
  - post similarity
  - plain GCN
  - different layer types
  - different intermediate embedding sizes
- save model from HPO pred and run introduce it to GCN stuff 

- discrepancy between len(protein\_list) on kw60407 and glogin!!
  - size of dti graph
- fix molecular predictor issue

Results:
- Flat distr intersec:  81.0 %AUROC
- Flat distr union:     83.1 %AUROC
- Flat combi intersec:  80.1 %AUROC
- GCN distr union:      78.3 %AUROC

Discussion:
- Sigmoid -> GCN -> Sigmoid 
  - output GCNConv weight for this issue
  - remove last Sigmoid
  - reread GCN paper


## TODO:

- utilize edge weight
- increase PPI graph min score

- use molecular features
  - pretrain
  - add to GCN model
- test other layer types
- DenseConvolution
- normalization
  - Graph norm
  - batch norm

1. Check whether GCN doesn't put anything on themselves
2. Fix molecular predictor
3. check whether it localizes on graph
4. What is pooling exactly

- GENConv
  - MsgNorm = True
  - learn msgscale = True

## Some questions
- is it feasible to have the test entities already in the graph already?


## Additional stuff
- test on Saras driver gene cancer
- test for gene OMIM stuff

- find and test a gold standard 

## TODO 
1. use drug indication file for mappings from indications to UMLS id
2. fuse yamanishi drug-indications and M's stuff
3. execute umls2owl 
4. start DL2vec on these indications
  4.1. add as new mode to not completely wipe the existing embeddings
5. integrate onto model
 



- Analysis of subclasses
- cherry pick drug and discuss that very drug 

- intro, method results discussion

- methods:
  - datasets
    - mappings
      - how much data is lost
  - hyperparameter tuning 
    - why is model looking like it is
      - formulas
    - grid search
    - split implementation
    - size of dataset
    - evaluation
    - loss
  - training process
    - split implementation
  - evaluation method
    1. how did I evaluate (evaluation metrics)
    2. comparative evaluation (comparative evaluation/comparison)
    - how long for training
    - why did i stop
    - size of dataset
    - evaluation
    - loss
  - discussion drug-driver genes

- results
  - description of method (high level) + figure for workflow
    - da kommen diese features her, da die anderen, methode generisch
  - experiments
    - "I took this dataset and did this."
    - Table(ROC-Curve)
      - first features plain
      - then molecular
      - then phenotype
      - then with network
    - what did I learn?
  - comparison
    - Yamanishi benchmark
    - comparison table 

- discussion
  - what is new
    - reached goal of combining graph with ..
  - driver gene stuff
    - cherry pick
  - technical discussion
    - what did I learn
      - softmax
    - which features any why
    - many other features possible (What could be done different)
      - molecular features
        - why as chosen
      - ontology features
        - many other ontology methods
          - OWL2vec
          - GNN (end-to-end)
- conclusion (summary)
          


