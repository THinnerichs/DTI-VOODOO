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

from ..Ontology import AnnotationCorpus
from ..Ontology import Ontology
from SemSimUtils import *
import sys
import os
import math

# types of semantic similarity

P_SS = "Pairwise"
G_SS = "Groupwise"

"interdag"
"intradag"

	#-#-#-#-#-#-#-#-#
	# class SemSim  #
	#-#-#-#-#-#-#-#-#

class MissingAcException(Exception):
	def __init__(self, message = None):
		self.message = message
	#
#


class SemSim(object):
	"""
	This class provides the prototype for a generic Semantic Similarity measure (SS)
	There are two types of SemSim: those which evaluate the semantic similarity between two sets of terms (groupwise - G_SS), and those which can only evaluate the similarity between pairs of GO terms (pairwise - P_SS). Each class extending TermSemSim should declare whether it is groupwise or pairwise.

	SemSim relies on SemSimUtils to perform a lot of tasks (such as evaluating Term IC or common ancestors)
	a SemSimUtils object can be passed to the constructor as input data. Otherwise, a new instance will be created. Using only one copy of SemSimUtils helps reducing time and spece requirements and is strongly adviced.
	This class DOES NOT perform any check on the input query
	"""

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
	SemSim_type = None # SemSim_type # Type of Term Sem Sim: Can be P_SS or G_SS 
	IC_based = None # IC_based # Tells whether te Term Sem Sim is based on Information Content

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

	def __init__(self, ontology, ac = None, util = None, do_log = False):
		# self.SemSim_type = None # Type of Sem Sim: Can be P_SS or G_SS
		# self.IC_based = False # Tells whether te Sem Sim is based on Information Content
		if self.SemSim_type == None: # Must be set
			raise Exception
		if self.IC_based == None: # Must be set
			raise Exception

		self.ontology = ontology
		self.ac = ac
		self.util = util
		self.log = []
		self.do_log = do_log

		if self.IC_based and self.ac == None:
			raise MissingAcException()
		if self.util == None:
			self.util = SemSimUtils(self.ontology, self.ac)
		if self.IC_based and self.util.IC == None:
			self.util.det_IC_table()
	#

	def _SemSim(self, term1, term2):
		if self.do_log:
			reason = 'Generic _SemSim_'
			self.log.append(reason)
		return None
	#

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

	def SemSim(self, term1, term2, ontology = None):
	""" 
	This is the main function that should be called to evaluate the semantic similarity
	ONLY EVALUATES THE SEMANTIC SIMILARITY. NO CHECK OR CONVERSION IS PERFORMED
	It is necessary to overload this function for cross-ontological SS
	"""
		if self.do_log:
			self.log = []

		if isinstance(term1, None.__class__):
			if self.do_log:
				reason = 'Term 1 is None.'
				self.log.append(reason)
			return None
		if isinstance(term2, None.__class__):
			if self.do_log:
				reason = 'Term 2 is None.'
				self.log.append(reason)
			return None

		if self.SemSim_type == P_SS:
			if isinstance(term1, (dict, list)):
				raise Exception()
			if isinstance(term2, (dict, list)):
				raise Exception()

		if self.SemSim_type == G_SS:
			if not isinstance(term1, (list, tuple)):
				raise Exception()
			if not isinstance(term2, (list, tuple)):
				raise Exception()

		return self._SemSim(term1, term2)
	#
#
