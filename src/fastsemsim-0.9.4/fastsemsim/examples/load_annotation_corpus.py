# -*- coding: iso-8859-1 -*-
'''
Copyright 2011 Marco Mina. All rights reserved.

This file is part of fastSemSim

fastSemSim is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fastSemSim is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with fastSemSim.  If not, see <http://www.gnu.org/licenses/>.
'''

from fastsemsim.Ontology import ontologies
from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim import data
import sys
import os

if __name__ == "__main__":

	# Select the type of ontology (GeneOntology, ...)
	# ontology_type = 'GeneOntology'
	# ontology_type = 'CellOntology'
	ontology_type = 'DiseaseOntology'

	# Select the relatioships to be ignored. For the GeneOntology, has_part is ignore by default, for CellOntology, lacks_plasma_membrane_part is ignored by default
	# ontology_parameters =	{}
	ontology_parameters =	{'ignore':{}}
	# ontology_parameters =	{'ignore':{'has_part':True, 'occurs_in':True, 'happens_during':True}}
	# ontology_parameters =	{'ignore':{'regulates':True, 'has_part':True, 'negatively_regulates':True, 'positively_regulates':True, 'occurs_in':True, 'happens_during':True}}

	# Select the source file type (obo or obo-xml)
	ontology_source_type = 'obo'

	# Select the ontology source file name. If None, the default ontology_type included in fastsemsim will be used
	ontology_source = None

	# Select the ac source file name. If None, the default ac included in fastsemsim for the ac_species will be used
	ac_source = None

	ac_species = 'human'
	# ac_species = 'arabidopsis'
	# ac_species = 'fly'
	# ac_species = 'mouse'
	# ac_species = 'rat'
	# ac_species = 'worm'
	# ac_species = 'zebrafish'

	ac_source_type = 'plain'
	# ac_source_type = 'gaf-2.0'

	ac_params = {}

# gaf-2.0 ac
	ac_params['filter'] = {} # filter section is useful to remove undesired annotations
	ac_params['filter']['EC'] = {} # EC filtering: select annotations depending on their EC
	# ac_params['filter']['EC']['EC'] = EC_include # select which EC accept or reject
	# ac_params['filter']['EC']['inclusive'] = True # select which EC accept or reject
	# ac_params['filter']['EC'] = {} # EC filtering: select annotations depending on their EC
	# ac_params['filter']['EC']['EC'] = EC_ignore # select which EC accept or reject
	# ac_params['filter']['EC']['inclusive'] = False # select which EC accept or reject

	ac_params['filter']['taxonomy'] = {}
	# ac_params['filter']['taxonomy']['taxonomy'] = tax_include # set properly this field to load only annotations involving proteins/genes of a specific species
	# ac_params['filter']['taxonomy']['inclusive'] = True # select which EC accept or reject
	# ac_params['filter']['taxonomy'] = {}
	# ac_params['filter']['taxonomy']['taxonomy'] = tax_ignore
	# ac_params['filter']['taxonomy']['inclusive'] = False # select which EC accept or reject
	# ac_params['simplify'] = True # after parsing and filtering, removes additional information such as taxonomy or EC. Useful if you have a huge amount of annotations and not enough memory

# Plain ac
	ac_params['multiple'] = True # Set to True if there are many associations per line (the object in the first field is associated to all the objects in the other fields within the same line)
	ac_params['term first'] = False # set to True if the first field of each row is a GO term. Set to False if the first field represents a protein/gene
	ac_params['separator'] = "\t" # select the separtor used to divide fields
	

	print "\n######################"
	print "# Loading ontology... #"
	print "######################\n"

	ontology = ontologies.load(source=ontology_source, source_type=ontology_source_type, ontology_type = ontology_type, parameters=ontology_parameters)

	print "\n######################"
	print "# Loading annotation corpus... #"
	print "######################\n"

	ac = AnnotationCorpus.AnnotationCorpus(ontology)
	builtin_dataset = data.dataset.Dataset()
	if ac_source is None:
		ac_data = builtin_dataset.get_default_annotation_corpus(ontology.name, ac_species)
		ac_source = ac_data['file']
		ac_source_type = ac_data['filetype']
	if ac_source is None:
		ac = None
	ac.parse(ac_source, ac_source_type, ac_params)
	ac.isConsistent()

	print "\n#################################"
	print "# Annotation corpus successfully loaded."
	print "#################################\n"


	print "Ontology source: " + str(ontology_source)
	print "Ontology source_type: " + str(ontology_source_type)
	print "ontology_type: " + str(ontology_type)
	print "ontology_parameters: " + str(ontology_parameters)
	print "Ontology - Number of nodes: " + str(ontology.node_number())
	print "Ontology - Number of edges: " + str(ontology.edge_number())
	print "\nOntology - Type and number of edges:\n-------------\n" + str(ontology.edges['type'].value_counts())
	print "-------------"
	print "\nOntology - Inner edge number (within the ontology):\n-------------\n" + str(ontology.edges['inner'].value_counts())
	print "-------------"
	print "\nOntology - Intra edge number (within the same namespace):\n-------------\n" + str(ontology.edges['intra'].value_counts())
	print "-------------"
	print "\nOntology - Outer edges (link to other ontologies):\n-------------\n" + str(ontology.edges.loc[ontology.edges['inner'] == False])
	print "-------------"
	print "\nOntology - Inter edges (link between different namespaces - within the same ontology):\n-------------\n" + str(ontology.edges.loc[(ontology.edges['intra'] == False) & (ontology.edges['inner'] == True)])
	print "-------------"
	#

	print "\n\n"
	print "AC source: " + str(ac_source)
	print "ac source_type: " + str(ac_source_type)
	print "ac_parameters: " + str(ac_params)
	print "AC species: " + str(ac_species)
	print "ac - Number of annotated proteins: " + str(len(ac.annotations))
	print "ac - Number of annotated terms: " + str(len(ac.reverse_annotations))
	print "-------------"
	#
#





