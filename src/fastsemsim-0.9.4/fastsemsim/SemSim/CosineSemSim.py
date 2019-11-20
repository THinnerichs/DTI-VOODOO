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
Cosine Semantic Similarity Measure

**Reference**: 
"""


from SemSimUtils import *
from TermSemSim import * 
# import sys
# import os
import math

class CosineSemSim(TermSemSim):
	SS_type = TermSemSim.G_TSS
	IC_based = False
	extend_annotations = True

	# def __init__(self, go, ac, util = None):
		# super(CosineSemSim, self).__init__(go, ac, util)
	
	def dotprod(self, vector1, vector2):
		dotprod = 0.0
		for i in range(0,len(vector1)):
			dotprod += vector1[i]*vector2[i]
		return dotprod
		
	def _SemSim(self, term1, term2):
		if self.extend_annotations:
			anc1 = self.util.get_ancestors(term1)
			anc2 = self.util.get_ancestors(term2)
		else:
			anc1 = term1
			anc2 = term2

		self.vector1 = []
		self.vector2 = []
		for i in anc1:
			self.vector1.append(1)
			if i in anc2:
				self.vector2.append(1) 
			else:
				self.vector2.append(0)
		for i in anc2:
			if i in anc1:
				pass
			else:
				self.vector2.append(1)
				self.vector1.append(0)

		num = self.dotprod(self.vector1, self.vector2)
		den1 = math.sqrt(self.dotprod(self.vector1, self.vector1))
		den2 = math.sqrt(self.dotprod(self.vector2, self.vector2))
		return float(num)/float(den1*den2)
	#
#
