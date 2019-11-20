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
Best Match Average (BMA) mixing strategy for pairwise term Protein Semantic Similarity measures
"""

from MixSemSim import *
# import sys
# import os
import math

class BMASemSim(MixSemSim):
	fair = True
	
	def _SemSim(self, scores):
		if len(scores) == 0:
			return 0
		rebuild1 = {}
		rebuild2 = {}
		for i in scores:
			#print(i
			#print(scores[i]
			#if i[2] == -1:
				#print("Errore in BMASemSim")
				#sys.exit()
			if i[0] not in rebuild1:
				rebuild1[i[0]] = {}
			if i[1] not in rebuild1[i[0]]:
				rebuild1[i[0]][i[1]] = i[2]
			if i[1] not in rebuild2:
				rebuild2[i[1]] = {}
			if i[0] not in rebuild2[i[1]]:
				rebuild2[i[1]][i[0]] = i[2]
		rescore1 = []
		rescore2 = []
		for i in rebuild1:
			temp_score = 0
			for j in rebuild1[i]:
				if rebuild1[i][j] >= temp_score:
					temp_score = rebuild1[i][j]
			rescore1.append(temp_score)
		for i in rebuild2:
			temp_score = 0
			for j in rebuild2[i]:
				if rebuild2[i][j] >= temp_score:
					temp_score = rebuild2[i][j]
			rescore2.append(temp_score)
		avg1 = 0
		for i in rescore1:
			avg1 += i
		avg2 = 0
		for i in rescore2:
			avg2 += i
		unfair_score = (avg1 + avg2) / (float(len(rescore1) + len(rescore2)))
		fair_score = (avg1/float(len(rescore1)) + avg2/float(len(rescore2))) / 2
		if self.fair:
			return fair_score
		return unfair_score
