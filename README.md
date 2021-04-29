# DTI-Voodoo: 
### Machine learning over interaction networks and ontology-based background knowledge predicts drug--target interactions

**Motivation:** *In silico* drug--target
  interaction (DTI) prediction is important for drug discovery and
  drug repurposing.  Approaches to predict DTIs can proceed
  indirectly, top-down, using phenotypic effects of drugs to identify
  potential drug targets, or they can be direct, bottom-up and use
  molecular information to directly predict binding potentials.  Both
  approaches can be combined with information about interaction
  networks.

**Results:** We developed DTI-Voodoo as a computational method
  that combines molecular features and ontology-encoded phenotypic
  effects of drugs with protein--protein interaction networks, and
  uses a graph neural networks to predict DTIs.  We demonstrate that
  drug effect features can exploit information in the interaction
  network whereas molecular features do not.  DTI-Voodoo is designed to
  predict candidate drugs for a given protein; we use this formulation
  to show that common DTI datasets contain intrinsic biases with major
  affects on performance evaluation and comparison of DTI prediction
  methods. Using a modified evaluation scheme, we demonstrate that
  DTI-Voodoo improves substantially over state of the art DTI prediction
  methods.

## Requirements
`Python 3.7` packages:
```
- pytorch 1.6+
- pytorch-geometric 1.6+ (please refer to https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html)
- numpy 1.19+
- scikit-learn 
- networkx
- gensim
- rdflib
- BioPython
- tqdm (for better evaluation and sanity preservation)
```

Others:
```
- Groovy (Groovy Version: 2.4.10 JVM: 1.8.0_121) with Grape for dependency management (http://docs.groovy-lang.org/latest/html/documentation/grape.html) for DL2vec axiom generation
- diamond & blastp (for DeepGO feature generation)
```

## How to run

Example commands for all used datasets are provided in the `src/` folder.

## Datasets
We provide all used datasets and download instructions in the `data/` and associated data preparation methods in the `src/` directory. 

