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
Dice Semantic Similarity Measure

**Reference**: 
"""

from SemSimUtils import *
from TermSemSim import * 
# import sys
# import os
import math

class DiceSemSim(TermSemSim):
	SS_type = TermSemSim.G_TSS
	IC_based = False
	extend_annotations = True

	# def __init__(self, go, ac, util = None):
		# super(DiceSemSim, self).__init__(go, ac, util)
		
	def _SemSim(self, term1, term2):
		if self.extend_annotations:
			inters = self.util.det_common_ancestors(term1, term2)
			anc1 = self.util.get_ancestors(term1)
			anc2 = self.util.get_ancestors(term2)
			allanc = len(anc1) + len(anc2)
		else:
			inters = self.util.intersection(term1, term2)
			allanc = len(term1) + len(term2)
		return float(len(inters))/float(allanc)