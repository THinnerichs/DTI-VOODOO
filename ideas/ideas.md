## Important

- Extend VPN access (only valid until 30.12.19)
- Visit ithelpdesk on Sunday for VPN credentials

## TODO

- actually do some literature search!!!

- Put more comments to code

- run FAMSA on Ibex

- Run Alignment on min_score 700 to avoid that many 0 target drugs
- execute hmm_pipeline on MSA results

- how to map e values from [0,inf) with 0 as hit and inf as no hit to [0,1] where 0 no hit and 1 is hit
- how to compare both similarity measures? (-> Maxat)


- get distibution over similarities and then take lowest, median and highest and fit above mapping to it

- iSCB conference deadline is 30. January


--------------------------------------------------------------------
- move repo to project folder (reach out to Ramy)

- Build alignment over all targets
  - build hmms over alignments
  - Query against only human target data
  - Build STRING graph only on human data (Y)
  - Use human interactions data to build training/testing data
    - Utilize non-binary score? Sigmoid over score?

- get node2vec embedding with node features
(- use node2vec and add my embedding on top of that) -> destroys idea of approach

- get phenotypes for drugs

- Finally execute hmm\_pipeline !!! 

- check if VPN is working

- Run actual model (-> Without node features at first)
  - When doing train/test split actually exclude test target? Given in from other species

./famsa-1.2.5-linux -t 4 ../data/fasta_files/CIDm00000753_fasta_700_min_score.fasta ../data/alignment_targets/CIDm00000753_famsa_aligned_700_min_score.afa
./famsa-1.2.5-linux -t 4 ../data/fasta_files/CIDm00000888_fasta_700_min_score.fasta ../data/alignment_targets/CIDm00000888_famsa_aligned_700_min_score.afa

---------------------------------------------------
- read Roberts Hypothesis testing paper for title and outline

- What are GCNs exactly doing?

- What do links in PPI graph mean? (-> STRING data/documentation)
  - what does closeness in enriched PPI graph means?
- Why do i need the DDIs and the Side effect similarity

- Read PREDICT paper again

- are there papers on DTI using PPIs?

- PPI radius based on confidence? -> Ego Graph



## Ideas for presentation

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

