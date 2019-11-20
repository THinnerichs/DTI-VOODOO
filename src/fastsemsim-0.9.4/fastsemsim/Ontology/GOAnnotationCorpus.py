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

'''
This class provides an interface to handle Annotation Corpora.

annotations: dict with protein ids as primary key. Each key is associated with a dictionary of GO Terms annotated for the protein. Detailed information, when available, are included as values within the latter dictionary.

reverse_annotations: dict with GO Terms as primary key. Each key is associated with a dict of proteins/gene products annotated with the GO term.

obj_set: set of proteins/gene products present in the annotation table, connected with the taxon id of the organism they belong to, when this information is available. This table is useful to filter out proteins from uninteresting species.

term_set: set of terms present in the annotation table.

If a GO object is passed as input data, annotation corpus is corrected removing obsolete annotations and resolving alternative ids.
This can be done later by calling sanitize method after supplying a valid GO object.

general_parameters: filtering options and parameters that apply in general
specific_parameters: parameter that should be used to load a particular file format

Each type of file carries different types of information. How to deal with that? Every operation is rerouted to the original file parser, that will take care of it. This is good since it avoids to duplicate data. 

'''

import sys
import copy

# from GeneOntology import *
from PlainAnnotationCorpus import PlainAnnotationCorpus
from GAF2AnnotationCorpus import GAF2AnnotationCorpus

INT_DEBUG = True
FILTER_PARAM = 'filter'
RESET_PARAM = 'reset'

AnnotationCorpusFormat = {'gaf-2.0':GAF2AnnotationCorpus,
													'gaf-1.0':None,
													'GOA':GAF2AnnotationCorpus,
													'plain':PlainAnnotationCorpus
													}

class GOAnnotationCorpus:
	#----------------------------------------------------------------------------------------
	int_exclude_GO_root = True

	#----------------------------------------------------------------------------------------

	def reset(self):
		self.annotations = {}
		self.reverse_annotations = {}
		self.obj_set = {}
		self.term_set = {}
		self.int_resetFields()
		self.filters = {}

	def int_resetFields(self):
		self.obj_fields = []
		self.term_fields = []
		self.annotations_fields = []
		self.reverse_annotations_fields = []
		self.obj_field2pos= {}
		self.term_field2pos= {}
		self.annotations_field2pos= {}
		self.reverse_annotations_field2pos= {}

	def __init__(self, go=None):
		self.go = go
		self.initCommonFilter()
		self.reset()

	def __deepcopy__(self, memo):
		a = AnnotationCorpus(self.go)
		a.exclude_GO_root = self.exclude_GO_root
		a.annotations = copy.deepcopy(self.annotations, memo)
		a.reverse_annotations = copy.deepcopy(self.reverse_annotations, memo)
		a.obj_set = copy.deepcopy(self.obj_set, memo)
		a.term_set= copy.deepcopy(self.term_set, memo)
		a.filter_taxonomy = copy.deepcopy(self.filter_taxonomy, memo)
		a.filter_EC = copy.deepcopy(self.filter_EC, memo)
		a.filter_EC_inclusive = self.filter_EC_inclusive
		a.filters = copy.deepcopy(self.filters, memo)
		a.obj_fields = self.obj_fields
		a.term_fields = copy.deepcopy(self.term_fields, memo)
		a.annotations_fields = copy.deepcopy(self.annotations_fields, memo)
		a.reverse_annotations_fields = copy.deepcopy(self.reverse_annotations_fields, memo)
		a.obj_field2pos= copy.deepcopy(self.obj_field2pos, memo)
		a.term_field2pos= copy.deepcopy(self.term_field2pos, memo)
		a.annotations_field2pos= copy.deepcopy(self.annotations_field2pos, memo)
		a.reverse_annotations_field2pos= copy.deepcopy(self.reverse_annotations_field2pos, memo)
		return a
		
#-----------------------------------------------------------------------------
# Load the annotation corpus from a file.
# - params is supposed to be a dict. Each key specifies a type of parameter (the associated value)
#   - currently supported parameters:
#     - FILTER_PARAM: parameters regarding the filtering of annotation corpus. See set_filters function for more details
#     - RESET_PARAM: if set to False, do not clean current annotations, but integrate the new data.
#-----------------------------------------------------------------------------
	
	def load(self, fname, ftype, params={}):
		self.parse(fname, ftype, params)

	def parse(self, fname, ftype, params={}):
		#print "AnnotationCorpus: parse"
		if params == None:
			self.reset()
		elif RESET_PARAM in params and params[RESET_PARAM]:
			self.reset()

		if not params == None and FILTER_PARAM in params:
			self.setCommonfilters(params[FILTER_PARAM])

		if ftype in AnnotationCorpusFormat:
			temp = AnnotationCorpusFormat[ftype](self, params)
			return temp.parse(fname)
		else:
			if INT_DEBUG:
				print "AnnotationCorpus.py: Format not recognized"
			raise Exception
#



#-----------------------------------------------------------------------------
# Consistency check and sanitizer, Align the annotation corpus to a Gene Ontology
# i.e. removing obsolete terms, mapping alternative ids to respective primary ones, ...
#-----------------------------------------------------------------------------
	
	def sanitize(self):
		if self.go is None:
			if INT_DEBUG:
				print("No GO specified. All the annotations will be considered valid.")
			return True
		for i in self.reverse_annotations.keys():
			if not i in self.go.alt_ids:
				#print("Term " + str(i) + " not found in GO.")
				for j in self.reverse_annotations[i]:
					del self.annotations[j][i]
				del self.reverse_annotations[i]
				continue
			if not i in self.go.nodes:
				if self.go.alt_ids[i] == i:
					#print("Term " + str(i) + " is an obsolete id.")
					for j in self.reverse_annotations[i]:
						del self.annotations[j][i]
						if len(self.annotations[j])==0:
							del self.annotations[j]
							del self.obj_set[j]
					del self.reverse_annotations[i]
					del self.term_set[i]
				else:
					#print("Term " + str(i) + " is an alternative id.")
					for j in self.reverse_annotations[i]:
						if not self.go.alt_ids[i] in self.annotations[j]:
							self.annotations[j][self.go.alt_ids[i]] = self.annotations[j][i]
						else:
							for k in self.annotations[j][i]:
								if not k in self.annotations[j][self.go.alt_ids[i]]:
									self.annotations[j][self.go.alt_ids[i]][k] = self.annotations[j][i][k]
						del self.annotations[j][i]
					del self.reverse_annotations[i]
					del self.term_set[i]
					continue
		return True

	def isConsistent(self):
		return self.int_checkConsistency()

	def int_checkConsistency(self):
		if self.go is None:
			if INT_DEBUG:
				print("No GO specified. All the annotations will be considered valid.")
			return True
		valid = True
		for i in self.reverse_annotations:
			if not i in self.go.alt_ids:
				if INT_DEBUG:
					print("Term " + str(i) + " not found in GO.")
				valid = False
				continue
			if not i in self.go.nodes:
				if self.go.alt_ids[i] == i:
					if INT_DEBUG:
						print("Term " + str(i) + " is an obsolete id.")
					#self.obsoletes[i] = {}
					valid = False
				else:
					if INT_DEBUG:
						print("Term " + str(i) + " is an alternative id.")
					valid = False
					continue
		return valid

#-----------------------------------------------------------------------------
# Filtering routines. Remove annotations that do not meet the requirements specified in *_field variables
#-----------------------------------------------------------------------------
	def setFilter(self, field, selector):
		self.filters[field] = selector

	def resetFilter(self, field):
		if field in self.filters:
			del self.filters[field]

	def isOk(self, field, value):
		if field not in self.filters:
			return True
		else:
			return self.filters[field](value)
			
	def constrain(self):
		if len(self.filters)==0:
			return

		for cf in self.filters:
			if cf in self.obj_field2pos:
				cp = self.obj_field2pos[cf]
				cff = self.filters[cf]
				temp_obj_set = {}
				for i in self.obj_set:
					if cff(self.obj_set[i]):
						temp_obj_set[i] = self.obj_set[i]
				self.obj_set = temp_obj_set
		for cf in self.filters:
			if cf in self.term_field2pos:
				cp = self.term_field2pos[cf]
				cff = self.filters[cf]
				temp_term_set = {}
				for i in self.term_set:
					if cff(self.term_set[i]):
						temp_term_set[i] = self.term_set[i]
				self.term_set = temp_term_set

		temp_annotations = {}
		temp_reverse_annotations = {}
		for i in self.obj_set:
			for j in self.annotations[i]:
				if j in self.term_set:
							if i not in temp_annotations:
								temp_annotations[i] = {}
							if j not in temp_annotations[i]:
								temp_annotations[i][j] = self.annotations[i][j]
							if j not in temp_reverse_annotations:
								temp_reverse_annotations[j] = {}
							if i not in temp_reverse_annotations[j]:
								temp_reverse_annotations[j][i] = self.reverse_annotations[j][i]
		self.annotations = temp_annotations
		self.reverse_annotations = temp_reverse_annotations

		temp_to_apply = []
		for cf in self.filters:
			if cf in self.annotations_field2pos:
				cp = self.annotations_field2pos[cf]
				cff = self.filters[cf]
				temp_to_apply.append((cff, cp))

		if len(temp_to_apply) > 0:
			temp_annotations = {}
			temp_reverse_annotations = {}
			for i in self.annotations:
				for j in self.annotations[i]:
					for k in self.annotations[i][j]: # assume k is ... what? a tuple or a key???
						temp_keep = True
						for ct in temp_to_apply:
							if not ct[0](k[ct[1]]): # version for list
								temp_keep = False
							#if not ct[0](self.annotations[i][j][k][ct[1]]): # version for dict
								#temp_keep = False
							if temp_keep:
								#del self.annotations[i][j][k]
								if i not in temp_annotations:
									temp_annotations[i] = {}
								if j not in temp_annotations[i]:
									temp_annotations[i][j] = {} # ? are you sure?
								temp_annotations[i][j][k] = self.annotations[i][j][k] # works for dict only, not for lists!
								if j not in temp_reverse_annotations:
									temp_reverse_annotations[j] = {}
								if i not in temp_annotations[j]:
									temp_reverse_annotations[j][i] = {} # ? are you sure?
								temp_reverse_annotations[j][i][k] = self.reverse_annotations[j][i][k] # works for dict only, not for lists!
			self.annotations = temp_annotations
			self.reverse_annotations = temp_reverse_annotations

			temp_obj_set = {}
			temp_term_set = {}
			for i in self.annotations:
				temp_obj_set[i] = {}
				for j in self.annotations[i]:
					if j not in temp_term_set:
						temp_term_set[j] = {}
			self.obj_set = temp_obj_set
			self.term_set = temp_term_set

#-----------------------------------------------------------------------------
# Filtering classes of common use.
#-----------------------------------------------------------------------------

	class TaxonomyFilter:
		name = 'taxonomy'
		taxonomy = {}
		inclusive = False

		def __init__(self, params):
			if 'taxonomy' in params:
				self.taxonomy = params['taxonomy']
				if type(self.taxonomy) == str or type(self.taxonomy) == unicode:
					self.taxonomy = {str(self.taxonomy):None}
			if 'inclusive' in params:
				self.inclusive = params['inclusive']

		def filter(self, taxonomy):
			#print self.taxonomy
			#print taxonomy
			if taxonomy in self.taxonomy and self.inclusive:
				return True
			elif not taxonomy in self.taxonomy and not self.inclusive:
				return True
			return False

	class ECFilter:
		name = 'EC'
		EC = {}
		inclusive = False
		def __init__(self, params):
			if 'EC' in params:
				self.EC = params['EC']
				if type(self.EC) == str or type(self.EC) == unicode:
					self.EC = {str(self.EC):None}
			if 'inclusive' in params:
				self.inclusive = params['inclusive']

		def filter(self, EC):
			if EC in self.EC and self.inclusive:
				return True
			elif not EC in self.EC and not self.inclusive:
				return True
			return False

	class GOFilter:
		name = 'GO'
		GO = None
		inclusive = True
		int_go = None
		def __init__(self, params=None):
			self.set(params)

		def set(self, params):
			if params == None:
				return True
			if 'int_current_go' in params:
				self.int_go = params['int_current_go']
			if 'GO' in params:
				self.GO = params['GO']
			if 'inclusive' in params:
				self.inclusive = params['inclusive']
			return True

		def filter(self, GO):
			if self.int_go == None:
				raise Exception
			if self.GO == None and not self.inclusive:
				return True
			elif self.GO == None and self.inclusive:
				return False
			print "Fix filter in GOFilter class."
			return True

#-----------------------------------------------------------------------------
# Simplified methods to set common filters
#-----------------------------------------------------------------------------
	def initCommonFilter(self):
		self.common_filters = {
											self.TaxonomyFilter.name:(self.TaxonomyFilter,),
											self.ECFilter.name:(self.ECFilter,),
											self.GOFilter.name:(self.GOFilter, {'int_current_go':self.go})
											}

	def setCommonfilters(self, inf):
		if type(inf) == dict:
			for i in inf:
				if i in self.common_filters:
					if len(self.common_filters[i]) > 1:
						new_filter = self.common_filters[i][0](self.common_filters[i][1])
						new_filter.set(inf[i])
					else:
						new_filter = self.common_filters[i][0](inf[i])
					self.setFilter(i, new_filter.filter)
		else:
			raise Exception

	def resetCommonfilter(self, i):
		self.resetFilter(i)
		
	#def reset_EC_filter(self):
		#self.resetFilter('EC')
#
