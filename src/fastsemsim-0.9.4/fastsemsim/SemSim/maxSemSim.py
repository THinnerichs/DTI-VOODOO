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
Max mixing strategy
"""

from MixSemSim import *
import sys
import os
import math

class maxSemSim(MixSemSim):

	def _SemSim(self, scores):
		if len(scores) == 0:
			return 0
		finale = 0
		for i in scores:
			#if i[2] == -1:
				#print("Errore in avgSemSim")
			if i[2] >= finale:
				finale = i[2]
		return finale
