This is a draft readme file that will be extended as more users
download STITCH and ask questions.

CHEMICAL SET
============

Chemicals are derived from PubChem. As described in the STITCH paper,
we merge salt forms and isomers. However, since STITCH 3, isomeric
compounds can also be investigated separately. In the download files,
the following convention holds:

CIDs / CID0... - this is a stereo-specific compound, and the suffix is the 
PubChem compound id.

CIDm / CID1... - this is a "flat" compound, i.e. with merged stereo-isomers
The suffix (without the leading "1") is the PubChem compound id.

(Note that the download files contain the prefixes CIDs/CIDm, while the API
still returns CID0/CID1.)

PROTEIN SET
===========

Proteins are from STRING 10 / 10.5. The format of the identifier is e.g.
9606.ENSP000123456 where 9606 is the NCBI taxonomic identifier and
ENSP... is the source name in the database which we imported.

To map our proteins to your proteins, please use the sequences
and identifiers (aliases) found on the STRING 9 download page:

http://string9.embl.de/newstring_cgi/show_download_page.pl


SCORES
======

The scores in the flat files are the same those on the website,
multiplied by 1000 (to make it a nice integer instead of an ugly
floating-point number). 

The scores are combined in a naive Bayesian fashion, i.e.

total Score = 1 - (1 - p1) * (1 - p2) * (1 - p3)

However, we do implement some further tweaks like correcting for the
random expectation (prior) of seeing interactions, as a naive Bayesian
combination of scores would overestimate the effect of very low 
contributions from multipleesources.


DETAILED SCORES
===============

The scores for different kinds of sources are combined and listed
in a column of the detailed scores file. For example, the database
column combines interactions from KEGG, DrugBank etc. 

For further information about the included sources, please check the 
details of interactions by clicking on edges in the network, or refer
to our papers:

http://nar.oxfordjournals.org/cgi/content/full/gkp937v1
http://nar.oxfordjournals.org/cgi/content/full/gkm795v1


ACTIONS
=======

Actions are only available for a subset of interactions, therefore
you will find many interactions that don't have a corresponding 
action in this file. 

In the actions file, mode is one of:

activation
phenotype (phenotypic effects, or: predicted to have the same phenotype)
binding
pred_bind (predicted to bind)
catalysis
inhibition
reaction

"action" can be only activation and inhibition. The idea is that you
can have a source that tells you "this is binding, and activating at
the same time". This is different from separate sources for "binding"
and "activation" because the activation might be via a different route
(e.g. changes in gene expression).

Actions are listed twice (A-B, B-A), in order to distinguish cases where
one of the partners is known to be the active partner (e.g. A inhibits B:
here a_is_acting would be "t" for A-B, but not for B-A).

If a_is_acting is "f", you can't infer that B is acting: sometimes none
of them is acting. For example, when a binding constant is measured,
we don't assign an active interaction partner.