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
import sys
import os

if __name__ == "__main__":

	# Select the type of ontology (GeneOntology, ...)
	ontology_type = 'GeneOntology'
	# ontology_type = 'CellOntology'
	# ontology_type = 'DiseaseOntology'

	# Select the relatioships to be ignored. For the GeneOntology, has_part is ignore by default, for CellOntology, lacks_plasma_membrane_part is ignored by default
	# ignore_parameters =	{}
	ignore_parameters =	{'ignore':{}}
	# ignore_parameters =	{'ignore':{'has_part':True, 'occurs_in':True, 'happens_during':True}}
	# ignore_parameters =	{'ignore':{'regulates':True, 'has_part':True, 'negatively_regulates':True, 'positively_regulates':True, 'occurs_in':True, 'happens_during':True}}

	# Select the source file type (obo or obo-xml)
	source_type = 'obo'

	# Select the source file name. If None, the default GeneOntology included in fastsemsim will be used
	source = None

	print "\n######################"
	print "# Loading ontology... #"
	print "######################\n"

	ontology = ontologies.load(source=source, source_type=source_type, ontology_type = ontology_type, parameters=ignore_parameters)

	print "\n#################################"
	print "# Ontology successfully loaded."
	print "#################################\n"

	print "source: " + str(source)
	print "source_type: " + str(source_type)
	print "ontology_type: " + str(ontology_type)
	print "ignore_parameters: " + str(ignore_parameters)
	print "Number of nodes: " + str(ontology.node_number())
	print "Number of edges: " + str(ontology.edge_number())
	print "\nType and number of edges:\n-------------\n" + str(ontology.edges['type'].value_counts())
	print "-------------"
	print "\nInner edge number (within the ontology):\n-------------\n" + str(ontology.edges['inner'].value_counts())
	print "-------------"
	print "\nIntra edge number (within the same namespace):\n-------------\n" + str(ontology.edges['intra'].value_counts())
	print "-------------"
	print "\nOuter edges (link to other ontologies):\n-------------\n" + str(ontology.edges.loc[ontology.edges['inner'] == False])
	print "-------------"
	print "\nInter edges (link between different namespaces - within the same ontology):\n-------------\n" + str(ontology.edges.loc[(ontology.edges['intra'] == False) & (ontology.edges['inner'] == True)])
	print "-------------"
	#

#
