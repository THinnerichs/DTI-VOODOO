#/bin/bash

# This command line runs fastsemsim comman-line with the parameters in file example_cmdline_params.

fastsemsim --load_params example_cmdline_params.txt

# The previous command is equivalent to running the following command:
# fastsemsim -vvv --task SS -o GeneOntology --ac fly  --query_input file --query_mode pairs --query_ss_type obj  --ontology_ignore positively_regulates --ontology_ignore negatively_regulates --ontology_ignore regulates --cut 0.9 --query_file data/GO_fly_pairs_query.txt  --tss Resnik --tmix BMA --root biological_process --ignore_EC TAS --query_file_sep "   "

# The file example_cmdline_params.txt was automatically built by fastSemSim by running this line adding the parameter  --save_params example_cmdline_params.txt:
# fastsemsim -vvv --task SS -o GeneOntology --ac fly  --query_input file --query_mode pairs --query_ss_type obj  --ontology_ignore positively_regulates --ontology_ignore negatively_regulates --ontology_ignore regulates --cut 0.9 --query_file data/GO_fly_pairs_query.txt  --tss Resnik --tmix BMA --root biological_process --ignore_EC TAS --query_file_sep "   " --save_params example_cmdline_params.txt

