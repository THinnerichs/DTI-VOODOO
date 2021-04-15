# DTI-Voodoo
**Motivation:** *In silico* drug--target interaction (DTI) prediction is crucial for drug discovery and down-stream tasks, such as drug indication and drug--disease prediction. However, DTI datasets are prone to certain biases, skewing the effective predictive expressiveness over those datasets. Additionally, only few approaches combine phenotypical and molecular features for DTI prediction.

**Results:** We present DTI-Voodoo to combine molecular features and functional information with protein--protein interaction networks under usage of graph convolutional neural networks, showing that phenotypic information localizes on the interaction network. We hereby focus on protein-centric evaluation, predicting drugs that may target specific targets. We further show that common drug--target interaction datasets contain intrinsic biases, majorly affecting predictive performance of compared models. Additionally, we propose an modified evaluation scheme and metric to circumvent such skews and to better reflect the problem DTI prediction methods aim to solve. 

## Requirements
`Python 3.7` packages:
```python
- pytorch 1.6+
- pytorch-geometric 1.6+
- numpy 1.19+
- scikit-learn 
- rdflib
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

