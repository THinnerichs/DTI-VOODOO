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
This class provides the prototype for a generic mixing strategy for pairwise Term Semantic Similarity measures
"""

import sys
import os
import math

class MixSemSim(object):

	def __init__(self, ontology, ac, util = None, do_log=False):
		self.ontology = ontology
		self.annotation_corpus = ac
		self.util = util
		self.do_log = do_log

		#if self.util == None:
			#self.util = SemSimUtils(ac, go)
			#self.ssu.det_IC_table()
	#

	def _format_data(self, term1): # simply organize input data
		if type(term1) is list or type(term1) is dict or type(term1) is set:
			return term1
		else:
			return [term1,]
	#
		
	def SemSim(self, set1, set2, TSS):
		if self.do_log:
			self.log = []
		lset1 = self._format_data(set1)
		lset2 = self._format_data(set2)
		if lset1 is None or lset2 is None or len(lset1) == 0 or len(lset2) == 0:
			if self.do_log:
				reason = 'Null input Terms'
				self.log.append(reason)
			return None
		temp_scores = []
		for i in lset1:
			for j in lset2:
				newscore = TSS.SemSim(i,j)
				if newscore == None:
					if self.do_log:
						reason = 'None score'
						self.log.append(reason)
					continue
				temp_scores.append( (i, j, newscore) )
		if len(temp_scores) == 0:
			if self.do_log:
				reason = 'No scores available'
				self.log.append(reason)
			return None
		return self._SemSim(temp_scores)
	#

	def _SemSim(self, scores):
		if self.do_log:
			reason = 'Generic mixing strategy'
			self.log.append(reason)
		raise "No mixing strategy selected."
		return None
	#
#
