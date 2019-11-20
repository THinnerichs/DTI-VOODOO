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
Jiang and Conrath Semantic Similarity Measure

**Reference**: 
"""

# from fastSemSim.Ontology import AnnotationCorpus
# from fastSemSim.Ontology import Ontology
# from SemSimUtils import *
from TermSemSim import * 
# import sys
# import os
import math

class JiangConrathSemSim(TermSemSim) :
	SS_type = TermSemSim.P_TSS
	IC_based = True

	# def __init__(self, go, ac, util = None):
		# super(JiangConrathSemSim, self).__init__(go, ac, util)
		
	def _SemSim(self, term1, term2):
		termid = self.util.det_MICA(term1, term2)
		sim = 1/(-2.0 * self.util.IC[termid] + self.util.IC[term1] + self.util.IC[term2] + 1)
		return sim