# -*- coding: iso-8859-1 -*-

# Copyright 2011 Marco Mina. All rights reserved.

# This file is part of fastSemSim

# fastSemSim is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# fastSemSim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with fastSemSim.  If not, see <http://www.gnu.org/licenses/>.


# print "fastsemsim/data/dataset.py"

import pandas as pd
# import sys
import os

class Dataset(object):
	'''
	This class keeps track of the dataset of ontologies and annotation corpora included in fastSemSim.
	The file data/dataset.txt is read to collect the list of embedded ontologies and annotation corpora.
	'''

	def __init__(self, descriptor=None):
		'''
		Initialize class structures. Use the descriptor parameter to specify a dataset descriptor file.
		By default, the file data/dataset.txt will be used.
		'''
		# program_dir = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
		# print "dataset.py: " + program_dir
		self.populate(descriptor)
		# print self.dataset
	#

	def populate(self, descriptor=None):
		'''
		Initialize class structures. Use the descriptor parameter to specify a dataset descriptor file.
		By default, the file data/dataset.txt will be used.
		'''
		descriptor_null = False
		if descriptor == None:
			descriptor_null = True
		if descriptor_null:
			program_dir = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
			# print program_dir
			descriptor = program_dir + '/dataset.txt'
		self.dataset = pd.read_csv(descriptor, sep="\t", comment="#", header=0).dropna(how='all')
		self.dataset.index = self.dataset['name']
		if descriptor_null:
			program_dir = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")  + "/"
			self.dataset['file'] = program_dir + self.dataset['file']
			# print self.dataset
	#

	def get_default_ontology(self, ontology_type):
		'''
		Return the default embedded ontology of the ontology_type type
		'''
		selected = self.dataset.loc[(self.dataset['type'] == 'O') & (self.dataset['ontology'] == ontology_type)]
		if selected.shape[0] == 0:
			return None
		return selected.iloc[0] # return the first (preferred) ontology
		pass
	#

	def get_dataset(self, dataset_name):
		'''
		Return the required dataset
		'''
		if not dataset_name in self.dataset.index:
			return None
		return self.dataset.loc[dataset_name] # return the selected ontology or ac
	#

	def get_ontology(self, dataset_name):
		'''
		Return the required ontology
		'''
		if not dataset_name in self.dataset.index:
			return None
		descriptor = self.dataset.loc[dataset_name] # return the selected ontology or ac
		if descriptor['type'] == 'O':
			return descriptor # return the selected ontology
		return None

	def get_annotation_corpus(self, dataset_name):
		'''
		Return the required annotation corpus
		'''
		if not dataset_name in self.dataset.index:
			return None
		descriptor = self.dataset.loc[dataset_name] # return the selected ontology or ac
		if descriptor['type'] == 'AC':
			return descriptor # return the selected ac
		return None
	#

	def get_annotation_corpus_by_species(self, ontology=None, species=None):
		'''
		Return the annotation corpus for the selected species, and compatible with the ontology specified by the ontology parameter.
		'''
		if ontology == None and species == None:
			selected = self.dataset.loc[(self.dataset['type'] == 'AC') & (self.dataset['ontology'] == ontology) & (self.dataset['species'] == species)]
		elif ontology == None:
			selected = self.dataset.loc[(self.dataset['type'] == 'AC') & (self.dataset['species'] == species)]
		elif species == None:
			selected = self.dataset.loc[(self.dataset['type'] == 'AC') & (self.dataset['ontology'] == ontology)]
		else:
			selected = self.dataset.loc[(self.dataset['type'] == 'AC') & (self.dataset['ontology'] == ontology) & (self.dataset['species'] == species)]
		return selected
	#

	def get_default_annotation_corpus(self, ontology=None, species=None):
		'''
		Return the default annotation corpus for the selected species, and compatible with the ontology specified by the ontology parameter.
		'''
		selected = self.get_annotation_corpus_by_species(ontology, species)
		if selected.shape[0] == 0:
			return None
		return selected.iloc[0] # return the first (preferred) ontology
		pass
	#
#
