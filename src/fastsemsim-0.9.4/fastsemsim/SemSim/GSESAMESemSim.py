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
G-SESAME Semantic Similarity Measure

**Reference**: 
"""

# from fastSemSim.Ontology import AnnotationCorpus
# from fastSemSim.Ontology import Ontology
# from SemSimUtils import *
from TermSemSim import * 
# import sys
# import os
import math

class GSESAMESemSim(TermSemSim) :
	SS_type = TermSemSim.P_TSS
	IC_based = False
	is_a_score = 0.8
	part_of_score = 0.6
	regulates_score = 0.6
	pos_regulates_score = regulates_score
	neg_regulates_score = regulates_score
	generic_score = 0.5

	# def __init__(self, go, ac, util = None):
		# super(GSESAMESemSim, self).__init__(go, ac, util)
		
	def score_ancestors(self, term):
		processed = {}
		queue = []
		processed[term] = 1
		queue.append(term)
		while len(queue) > 0:
			t = queue.pop()
			for tp in self.util.ontology.parents[t]:
				if tp not in processed:
					queue.append(tp)
					processed[tp] = processed[t] * self.score_edge(tp, t)
		return processed

	def score_edge(self, tp, t): # t = child, tp = parent
		# print str(tp) + " " + str(t)
		for j in self.ontology.nodes[tp]:
			# print self.ontology.edges['nodes'][j]
			if self.ontology.edges.ix[j,'child'] == t: # can replace with find. Way faster!
				# print self.ontology.edges.ix[j,'type']
				if self.ontology.edges.ix[j,'type'] == 'is_a':
					return self.is_a_score
				elif self.ontology.edges.ix[j,'type'] == 'part_of':
					return self.part_of_score
				elif self.ontology.edges.ix[j,'type'] == 'regulates':
					return self.regulates_score
				elif self.ontology.edges.ix[j,'type'] == 'positively_regulates':
					return self.pos_regulates_score
				elif self.ontology.edges.ix[j,'type'] == 'negatively_regulates':
					return self.neg_regulates_score
				else:
					return self.generic_score
		print "Error"
		raise Exception

	def _SemSim(self, term1, term2):
		ca = self.util.det_common_ancestors(term1, term2)

		s1 = self.score_ancestors(term1)
		s2 = self.score_ancestors(term2)

		num = 0
		for i in ca:
			num += s1[i] + s2[i]
		den = 0
		for i in s1:
			den += s1[i]
		for i in s2:
			den += s2[i]
		return float(num)/float(den)