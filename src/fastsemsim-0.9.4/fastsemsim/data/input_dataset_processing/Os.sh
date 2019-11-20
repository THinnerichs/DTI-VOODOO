#!/bin/bash
# Retrieve ontologies

# wget --no-check-certificate https://cell-ontology.googlecode.com/svn/trunk/src/ontology/cl.obo
# mv cl.obo ../Os/CellOntology_full_`date "+%Y.%m.%d"`.obo

wget --no-check-certificate https://cell-ontology.googlecode.com/svn/trunk/src/ontology/cl-basic.obo
mv cl-basic.obo ../Os/CellOntology_`date "+%Y.%m.%d"`.obo

wget --no-check-certificate http://purl.obolibrary.org/obo/doid.obo
mv doid.obo ../Os/DiseaseOntology_`date "+%Y.%m.%d"`.obo

wget http://purl.obolibrary.org/obo/go.obo
mv go.obo ../Os/GeneOntology_`date "+%Y.%m.%d"`.obo

wget http://archive.geneontology.org/latest-termdb/go_daily-termdb.obo-xml.gz
mv go_daily-termdb.obo-xml.gz ../Os/GeneOntology_`date "+%Y.%m.%d"`.obo-xml.gz
