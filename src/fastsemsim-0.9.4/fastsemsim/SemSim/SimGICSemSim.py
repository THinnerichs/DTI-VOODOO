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
SimGIC Semantic Similarity Measure

**Reference**: 
"""

# from fastSemSim.Ontology import AnnotationCorpus
# from fastSemSim.Ontology import Ontology
# from SemSimUtils import *
from TermSemSim import * 
import sys
import os
import math

class SimGICSemSim(TermSemSim) :
	SS_type = TermSemSim.G_TSS
	IC_based = True

	# def __init__(self, go, ac, util = None):
		# super(SimGICSemSim, self).__init__(go, ac, util)
		
	def _SemSim(self, term1, term2):
		#print term1
		#print term2
		inters = self.util.det_common_ancestors(term1, term2)
		union = self.util.det_ancestors_union(term1, term2)
		#print inters
		#print union
		intIC = 0
		for i in inters:
			intIC += self.util.IC[i]
		uniIC = 0
		for i in union:
			uniIC += self.util.IC[i]
		if uniIC == 0 and intIC == 0:
			return 0
		return intIC/uniIC
	#
#
