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
This class provides the prototype for Term semantic similarity measures (TSS)

There are two types of Term semantic similarity: a first group that can evaluate the semantic similarity between two sets of terms (groupwise - G_TSS), and a second group that can only evaluate the similarity between pairs of GO terms (pairwise - P_TSS). Each class extending TermSemSim should declare whether it is groupwise or pairwise.

TermSemSim relies on SemSimUtils to perform a lot of tasks (e.g. evaluating Term IC or common ancestors).
A SemSimUtils object can be passed to the constructor as input data. Otherwise, a new instance will be created. Using only one copy of SemSimUtils helps reducing time and space requirements and is strongly recommended.
'''

from fastsemsim.Ontology import AnnotationCorpus
from fastsemsim.Ontology import Ontology
from SemSimUtils import *
import sys
import os
import math


	#-#-#-#-#-#-#-#-#-#-#
	# class TermSemSim  #
	#-#-#-#-#-#-#-#-#-#-#

class MissingAcException(Exception):
	def __init__(self, message):
		self.message = message
	#
#

class TermSemSim(object):
	
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	format_and_check_data = True
	P_TSS = "Pairwise"
	G_TSS = "Groupwise"
	SS_type = None
	IC_based = None
	# SS_type # Type of Term Sem Sim: Can be P_TSS or G_TSS 
	# IC_based # Tells whether te Term Sem Sim is based on Information Content

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#### internal functions

	def __init__(self, ontology, ac = None, util = None, do_log = False):
		# self.SS_type = None # Type of Term Sem Sim: Can be P_TSS or G_TSS
		# self.IC_based = False # Tells whether te Term Sem Sim is based on Information Content
		if self.SS_type == None: # Must be set
			raise Exception
		if self.IC_based == None: # Must be set
			raise Exception

		self.ontology = ontology
		self.ac = ac
		self.util = util
		self.log = []
		self.do_log = do_log

		if self.IC_based and self.ac == None:
			raise MissingAcException("The selected semantic measure is based on IC and requires an annotation corpus.")
		if self.util == None:
			self.util = SemSimUtils(self.ontology, self.ac)
		if self.IC_based and self.util.IC == None:
			self.util.det_IC_table()
	#

	def _has_IC(self, term):
		# if self.util.IC == None:
			# return None
		if not term in self.util.IC:
			return False
		if self.util.IC[term] == None:
			return False
		return True
	#

	def _format_data(self, term1):
	# Format input query
	# 1) convert Terms to proper ontology format
	# 2) verify terms are valid (and have an IC)
	# 3) verify terms have the same root. (this blocks cross-ontological Term SemSim measures)
		id1 = self.ontology.id2node(term1, alt_check = False)

		if self.SS_type == self.P_TSS: # only single terms allowed
			if id1 == None:
				if self.do_log:
					reason = 'Unmapped Term'
					self.log.append(reason)
				return None
			if type(id1) is dict or type(id1) is list:
				if self.do_log:
					reason = 'Mulitple Terms passed to Pairwise SS Measure. Perhaps the term has multiple alternative ids.'
					self.log.append(reason)
				return None
			if not self.ontology.is_valid(id1):
				if self.do_log:
					reason = 'Invalid Term'
					self.log.append(reason)
				return None
			if self.IC_based and not self._has_IC(id1):
				if self.do_log:
					reason = 'No IC for Term'
					self.log.append(reason)
				return None
			return id1

		elif self.SS_type == self.G_TSS: # multiple terms allowed
			temp_id1 = []
			new_temp_id1 = []
			if type(id1) is dict:
				temp_id1 = id1.keys()
			elif type(id1) is list:
				temp_id1 = id1
			else:
				temp_id1 = [ id1 ]

			current_onto = None
			for i in temp_id1:
				if i == None:
					if self.do_log:
						reason = 'Unmapped Term'
						self.log.append(reason)
					continue
				if not self.ontology.is_valid(i):
					if self.do_log:
						reason = 'Invalid Term'
						self.log.append(reason)
					continue
				if self.IC_based and not self._has_IC(i):
					if self.do_log:
						reason = 'No IC for Term'
						self.log.append(reason)
					continue

				if current_onto is None: 
					current_onto = self.util.lineage[i]
				elif not current_onto == self.util.lineage[i]:
					if self.do_log:
						reason = 'Different namespaces'
						self.log.append(reason)
					return None

				new_temp_id1.append(i)

			if len(new_temp_id1) == 0:
				if self.do_log:
					reason = 'No Terms'
					self.log.append(reason)
				return None

			return new_temp_id1
	#

	def _SemSim(self, term1, term2):
		if self.do_log:
			reason = 'Generic _SemSim_'
			self.log.append(reason)
		return None
	#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#### public functions

	def setSanityCheck(self, en):
	# This enables or disables sanity checks and conversion on query data.
	# If speeds up computation, but if input data are not correct expect exceptions to be raised
		self.format_and_check_data = en


	def SemSim(self, term1, term2, ontology = None):
	# This is the main function that should be called to evaluate the Term Sem Sim
	# It takes care of verifying data, format them, and evaluate the Sem Sim.
	# It might be necessary to Overload this function for cross-ontological Term Sem Sim measures 
		if self.do_log:
			self.log = []
		# reason = 'Generic _SemSim_'
		# self.log.append(reason)

		if self.format_and_check_data:
			if term1 is None or term2 is None:
				if self.do_log:
					reason = 'Null input Terms'
					self.log.append(reason)
				return None
			id1 = self._format_data(term1)
			id2 = self._format_data(term2)
			#print "\""+term1+"\""
			#print "\""+term2+"\""
			#print id1
			#print id2
			if id1 is None or id2 is None or (self.SS_type == self.G_TSS and len(id1) == 0) or (self.SS_type == self.G_TSS and len(id2) == 0):
				#print(str(term1) + " or " + str(term2) + "   not valid.")
				# if self.do_log:
					# reason = 'Invalid Terms'
					# self.log.append(reason)
				return None
			if self.SS_type == self.P_TSS:
				if not self.util.lineage[id1] == self.util.lineage[id2]:
					#raise "Terms are not from the same ontology"
					if self.do_log:
						reason = 'Different namespaces between query terms'
						self.log.append(reason)
					return None
			elif self.SS_type == self.G_TSS:
				for i in id1:
					t1 = i
					break
				for i in id2:
					t2 = i
					break
				if not self.util.lineage[t1] == self.util.lineage[t2]:
					#raise "Terms are not from the same ontology"
					if self.do_log:
						reason = 'Different namespaces between query terms'
						self.log.append(reason)
					return None
		else:
			id1 = term1
			id2 = term2
		return self._SemSim(id1, id2)
	#
#
