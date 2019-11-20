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


"""
This class provides the prototype for a generic Object Set Semantic Similarity measure (PSS)
"""

# from fastSemSim.Ontology import AnnotationCorpus
# from fastSemSim.Ontology import Ontology
from SemSimUtils import SemSimUtils
# from TermSemSim import *
# from MixSemSim import *
from . import select_mix_SemSim
from . import select_term_SemSim

# import sys
# import os
import math

class ObjSetSemSim:
	# def pick_TSS(self):
	# 	if not self.TSS in SemSimMeasures:
	# 		raise "Semantic Similarity Measure not available."
	# 		return TermSemSim(self.ac, self.ontology, self.util)
	# 	else:
	# 		return SemSimMeasures[self.TSS][0](self.ac, self.ontology, self.util)

	# def pick_mixSS(self):
	# 	if not self.mixSS in MixingStrategies:
	# 		raise "Mixing Strategy not available."
	# 		return MixSemSim(self.ac, self.ontology)
	# 	else:
	# 		return MixingStrategies[self.mixSS](self.ac, self.ontology)

	def __init__(self, ontology, ac, TSS = None, MSS = None, util = None, do_log = False):
	# def __init__(self, ontology, ac, TSS = None, MSS = None, OMSS = None, util = None, do_log = False):
		
		self.ontology = ontology
		self.ac = ac
		self.do_log = do_log
		self.log = []

		self.util = util
		if self.util == None:
			self.util = SemSimUtils(self.ontology, self.ac)
			# self.util.det_IC_table()

		self.term_SS_class = select_term_SemSim(TSS)
		self.term_SS = self.term_SS_class(self.ontology, self.ac, self.util)

		if not MSS == None:
			self.mix_SS_class = select_mix_SemSim(MSS)
			self.mix_SS = self.mix_SS_class(self.ontology, self.ac, self.util)
		else:
			self.mix_SS = None

		# if not OMSS == None:
		# 	self.o_mix_SS_class = select_mix_SemSim(OMSS)
		# 	self.o_mix_SS = self.mix_SS_class(self.ontology, self.ac, self.util)
		# else:
		# 	self.o_mix_SS = None

		# self.TSS = TSS
		# self.mixSS = MSS
		# if self.TSS is None:
		# 	self.TSS = TermSemSim(self.ac, self.ontology, self.util)
		# elif type(self.TSS) is str or unicode:
		# 	self.TSS = str(self.TSS)
		# 	self.TSS = self.pick_TSS()
		# else:
		# 	raise Exception
		# if self.mixSS is None:
		# 	self.mixSS = MixSemSim(self.ac, self.ontology)
		# elif type(self.mixSS) is str or unicode:
		# 	self.mixSS = str(self.mixSS)
		# 	self.mixSS = self.pick_mixSS()
		# else:
		# 	raise Exception
	#

	def _format_data(self, obj, onto):
		if not type(obj) == list:
			obj = [obj,]
		terms = []
		for j in obj:
			if not j in self.ac.annotations:
				#print(str(obj) + " not found in Annotation Corpus.")
				if self.do_log:
					reason = 'Object not in annotation corpus'
					self.log.append(reason)
				continue
			
			for i in self.ac.annotations[j]:
				#if i in self.ontology.obsolete_ids: # not present in GO_root
					#continue
				# print i
				# print self.util.lineage[i]
				# print onto
				# print "--"
				if i in self.util.lineage and self.util.lineage[i] == onto:
					terms.append(i)
		return terms

	def _SemSim(self, term1, term2):
		if term1 is None or term2 is None:
			if self.do_log:
				reason = 'Invalid Term set'
				self.log.append(reason)
			return None
		if len(term1) == 0 or len(term2) == 0:
			if self.do_log:
				reason = 'Empty Term set'
				self.log.append(reason)
			return None
		if self.term_SS.SS_type == self.term_SS.P_TSS:
			sscore = self.mix_SS.SemSim(term1, term2, self.term_SS)
		elif self.term_SS.SS_type == self.term_SS.G_TSS:
			sscore = self.term_SS.SemSim(term1, term2)
		else:
			raise "Semantic Similarity measure not properly configured."
		return sscore
	#

	def SemSim(self, obj1, obj2, root = None):
		if self.do_log:
			self.log = []
		# print root
		if root == None:
			root = self.ontology.roots.keys()[0]
			# print root
		if not root in self.ontology.roots:
			# raise Exception(str(root) + " is not an ontology root.")
			if self.do_log:
				reason = 'Selected root ' + str(root) + ' is not in the ontology.'
				self.log.append(reason)
			return None
		# print obj1
		# print obj2
		t1 = self._format_data(obj1, root)
		t2 = self._format_data(obj2, root)
		# print t1
		# print t2
		return self._SemSim(t1, t2)
	#
#
