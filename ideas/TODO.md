## TODO
- acutally do some literature search!!!

- Put more comments to code

- Test FAMSA(-> get lib fixed)
- run Mafft 

- Run Alignment on min_score 700 to avoid that many 0 target drugs
- execute hmm_pipeline on MSA results


- calculate side effect similarity matrix (See jupyter notebook)
- get similarity matrix from bipartite graph only (Jaccard, others suggested by Maxat)

- How to stitch together MedDRA and SIDER graph in networkx (-> visit Maxat)
  - SPARQL against MEDDRA.ttl to get UMLS identifier for each MedDRA entity (preferably mapping)
  - Add my SIDER graph to MEDDRA.ttl in RDF format (watch out for correct URIs)

- how to compare alignments (use blastp for e values)
- how to map e values from [0,inf) with 0 as hit and inf as no hit to [0,1] where 0 no hit and 1 is hit
- how to compare both similarity measures? (-> Maxat)

- get distibution over similarities and then take lowest, median and highest and fit above mapping to it

