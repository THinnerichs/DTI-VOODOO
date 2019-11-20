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
Czekanowski and Dice Semantic Similarity Measure

**Reference**: 
"""

from SemSimUtils import *
from TermSemSim import * 
# import sys
# import os
import math

class CzekanowskiDiceSemSim(TermSemSim):
	SS_type = TermSemSim.G_TSS
	IC_based = False
	extend_annotations = True

	# def __init__(self, go, ac, util = None):
		# super(CzekanowskiDiceSemSim, self).__init__(go, ac, util)
		
	def _SemSim(self, term1, term2):
		if self.extend_annotations:
			anc1 = self.util.get_ancestors(term1)
			anc2 = self.util.get_ancestors(term2)
		else:
			anc1 = term1
			anc2 = term2
		inters = self.util.intersection(anc1, anc2)
		union = self.util.int_merge_sets(anc1, anc2)
		allanc = self.util.int_merge_sets(self.util.difference(anc1, anc2) , self.util.difference(anc2, anc1))
		return float(len(allanc))/float(len(union) + len(inters))