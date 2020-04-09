## TODO

1. Maybe check for some other chemical features for nodes for PPI graph 
2. Check some better fitting data for drugs -> e.g. SMILES and learn from them
3. Add best filter feature to node features
4. Add upsampling to data handling
5. Write draft for paper -> Put all ideas in there ask Robert for more
6. Check performance of filters -> Copy approach from dti\_predictor.py


Questions for Robert:
- Models
  - node classification vs graph classification
  - inductive vs. transductive: what exactly?
  - discuss models in detail
  - weakness of filters

- Data
  - What could be node features fitting better for bottom up approach?
  - Compare against results of other approaches on what dataset? Which ones got dataset? Union sufficient? 
    - HAVE LOOK AT SOTA PAPERS!!!!
  - Maybe use better features for drugs?
 


---------------------------------------------------

1. Perform hmm\_search pipeline and extract features
2. Build dummy features and get torch network to run
3. put in real features

---------------------------------------------------
- Get other method to run on my dataset (literature!!)
  - Have to get SMILES and stuff too -> meh
  - How can this be made proof? Model is tweaked on this dataset while foreign methods are not

----------------------------------------------------
- extract predictions from search results

- update targets based on new alignments
  - update after new hmm\_search went trough :)

---------------------------------------------------
- DTI\_graph:
  - min\_score:
    - > 700: 18987
    - > 400: 33549
    - all: 171168

---------------------------------------------------

- Find SOTA for DTI
  - take intersection of data and perform on that
  - Maybe one top-down and one bottom up?

--------------------------------------------------

Talk with Robert on 09.01.:

Work until now:
  - performed mafft and FAMSA alignment with search -> poor results
    1. tweak this! because its new!
    2. take other bottom up approach to fix this
  - Tweaking GCN, GAT model -> not really sufficient, performs horrible 
    - heavy problems getting that to run -> decompose stellargraph methods
    - other library?! -> you had time over christmas
      - Top k pooling? pooling layers in general?
      - Pytorch geometric

  - artificially reduce PPI > 700 score

  - Take Proteins without links into model at the end, not for tuning
  
  - used DTIs accurate? currently just taking all with score aboth 700 (just captures cooccurence)
  - used PPIs accurate? Currently taking all of them


## Ideas for presentation

- include numbers of interactions
  - histograms?
- Other approaches lack ability to generalize -> guilt by association

- histogram over scores of STRING and STITCH
- show image of network for STRING and STITCH
- Meaning of score 
- Functions to slides B)


## Other ideas

- Calculate similarity between drugs and metabolites...
	... over structural similarity
	... over alignment of e.g. SMILES

- get DD similarity trough ATC code/hierarchy?

## On performance of alignments

This [link](https://www.ebi.ac.uk/Tools/msa/)

## Other sources for MSAs

- [wiki article on MSAs](https://en.wikipedia.org/wiki/List_of_sequence_alignment_software#Multiple_sequence_alignment)

## Questions

- What is the actual hypothesis? 
- What could be some meat to add here?

## Databases

- CTD

Tray Eidecke

