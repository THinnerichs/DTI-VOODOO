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


'''
This class provides some routines to calculate basic properties used by different SS measures.
In particular this class provides code for evaluating:

- term ICs
- term frequency within an annotation corpus
- term's ancestors
- term's offspring
- terms's children
- terms's parents
- MICA/DCA/LCA
- term's distance
'''

# TODO:
# - term depth
# - terms common ancestors (between two terms or two sets of terms)
#             I could use pairs module to to this!


# from fastSemSim.Ontology import AnnotationCorpus
# from fastSemSim.Ontology import Ontology
# import sys
# import os
import math
import numpy as np

class SemSimUtils(object):

# variables

	# go = None
	# ac = None
	# ancestors = None
	# offspring = None
	# IC = None
	# freq = None
	# GO_root = None
	# p = None

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#internal functions

	def __init__(self, ontology, ac=None):
		self.ontology = ontology
		if self.ontology == None:
			raise Exception
		# self.ontology._s1_to_s2() # already done when loading the go. no need here?

		self.ac = ac
		
		self.ancestors = None
		self.offspring = None
		self.IC = None
		self.freq = None
		self.p = None
		self.lineage = None

		self.int_det_offspring_table()
		self.int_det_ancestors_table()
		self.int_det_lineage()

		#self.det_freq_table()
		#self.det_ICs_table()

	def int_det_offspring_table(self):
		self.offspring = {}
		# conta  = 0
		temp_intra = set(np.where(self.ontology.edges['intra'])[0])
		for i in self.ontology.nodes:
			# conta += 1
			# print str(conta) + "on " + str(len(self.ontology.nodes)) 
			self.offspring[i] = self.int_det_offspring(i, temp_intra)

	def int_det_ancestors_table(self):
		self.ancestors = {}
		temp_intra = set(np.where(self.ontology.edges['intra'])[0])
		for i in self.ontology.nodes:
			self.ancestors[i] = self.int_det_ancestors(i, temp_intra)

	def int_det_offspring(self, goid, temp_intra):
		if goid not in self.ontology.children:
			return set()
		anc = set()
		# anc.add(goid)
		processed = {}
		queue = [goid]
		# for i in self.ontology.children[goid]:
			# if self.ontology.children[goid][i] in self.ontology.edges['inter'] and self.ontology.edges['inter'][self.ontology.children[goid][i]]:
				# continue
			# queue.append(i)
		while len(queue) > 0:
			# print queue
			t = queue.pop()
			anc.add(t)
			for tp in self.ontology.children[t]:
				edid = self.ontology.children[t][tp]
				if not edid in temp_intra:
				# if self.ontology.children[t][tp] in self.ontology.edges['inter'] and self.ontology.edges['inter'][self.ontology.children[t][tp]]:
					continue
				if tp not in processed:
					queue.append(tp)
			processed[t] = None
		return anc

	def int_det_ancestors(self, goid, temp_intra):
		if goid not in self.ontology.parents:
			return set()
		anc = set()
		# anc.add(goid)
		processed = {}
		queue = [goid]
		# for i in self.ontology.parents[goid]:
			# queue.append(i)
		#print(queue
		while len(queue) > 0:
			t = queue.pop()
			anc.add(t)
			#print(parent_going[t]
			for tp in self.ontology.parents[t]:
				edid = self.ontology.parents[t][tp]
				if not edid in temp_intra:
				# if self.ontology.parents[t][tp] in self.ontology.edges['inter'] and self.ontology.edges['inter'][self.ontology.parents[t][tp]]:
					continue
				if tp not in processed:
					queue.append(tp)
			processed[t] = 0
		return anc

	def int_det_lineage(self):
		self.lineage = {}
		for j in self.ontology.roots:
			temp = self.offspring[j]
			for i in temp:
				if i in self.lineage:
					# print i
					# print self.lineage[i]
					raise Exception
				self.lineage[i] = j
		return self.lineage

	def int_det_freq(self,term_id):
		freq = 0
		children_set = self.offspring[term_id]
		for j in children_set:
			if j in self.ac.reverse_annotations:
				freq += len(self.ac.reverse_annotations[j])
		return freq

	def int_det_freq_table(self):
		self.freq = {}
		for i in self.ontology.nodes:
			self.freq[i] = self.int_det_freq(i)

	def int_det_p_table(self):
		self.p = {}
		for i in self.ontology.nodes:
			self.p[i] = self.int_det_p(i)

	def int_det_p(self,term_id):
		if self.freq == None:
			self.int_det_freq_table()
		if not term_id in self.freq:
			return None
		if self.freq[term_id] == float(0):
			return float(0)
		rootf = self.freq[self.lineage[term_id]]
		temp_p = float(self.freq[term_id])/float(rootf)
		return temp_p

	def int_det_IC(self, term_id):
		pr = self.int_det_p(term_id)
		if pr == None:
			return None
		if pr == float(0):
			return None
		pr = -math.log(pr)
		if pr == float(-0):
			pr = -pr 
		return pr

	def int_det_IC_table(self):
		self.IC = {}
		conta = 0
		for i in self.ontology.nodes:
			conta+= 1
			#print conta
			#print len(self.ontology.nodes_edges)
			temp_IC = self.int_det_IC(i)
			#if not temp_IC == None:
			self.IC[i] = temp_IC

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# public functions

# ASSUMPTION: Terms are passed as integers or lists of integers, not as strings!

	def det_IC_table(self):
		self.int_det_IC_table()
		return self.IC

	def det_IC(self,term):
		id = self.ontology.name2id(term)
		if self.IC == None:
			self.det_IC_table()
		if id in self.IC:
			return self.IC[id]
		return None

	def det_MICA(self, term1, term2):
		gene1anc = self.ancestors[term1]
		gene2anc = self.ancestors[term2]
		maxIC = 0
		maxterm = -1
		for i in gene1anc:
			if i in gene2anc:
				if self.IC[i] >= maxIC:
					maxIC = self.IC[i]
					maxterm = i
		return maxterm

	def int_merge_sets(self, set1, set2):
		ca = {}
		for i in set1:
			ca[i] = None
		for i in set2:
			ca[i] = None
		return ca

	def intersection(self, set1, set2):
		ca = {}
		for i in set1:
			if i in set2:
				ca[i] = None
		return ca

	def difference(self, set1, set2):
		ca = {}
		for i in set1:
			if not i in set2:
				ca[i] = None
		return ca

	def det_common_ancestors(self, term1, term2):
		if type(term1) is int or type(term1) is str:
			gene1anc = self.ancestors[term1]
		else:
			gene1anc = {}
			for i in term1:
				gene1anc = self.int_merge_sets(gene1anc, self.ancestors[i])
		if type(term2) is int or type(term2) is str:
			gene2anc = self.ancestors[term2]
		else:
			gene2anc = {}
			for i in term2:
				gene2anc = self.int_merge_sets(gene2anc, self.ancestors[i])
		#gene1anc = self.ancestors[term1]
		#gene2anc = self.ancestors[term2]
		ca = {}
		for i in gene1anc:
			if i in gene2anc:
				ca[i] = None
		return ca

	def get_ancestors(self, term1):
		if type(term1) is int or type(term1) is str:
			if term1 not in self.ancestors:
				return {}
			gene1anc = self.ancestors[term1]
		else:
			gene1anc = {}
			for i in term1:
				if i not in self.ancestors:
					continue
				gene1anc = self.int_merge_sets(gene1anc, self.ancestors[i])
		ca = {}
		for i in gene1anc:
			ca[i] = None
		return ca
		
	def det_ancestors_union(self, term1, term2):
		if type(term1) is int or type(term1) is str:
			if term1 not in self.ancestors:
				return {}
			gene1anc = self.ancestors[term1]
		else:
			gene1anc = {}
			for i in term1:
				if i not in self.ancestors:
					continue
				gene1anc = self.int_merge_sets(gene1anc, self.ancestors[i])
		#print(gene1anc
		if type(term2) is int or type(term2) is str:
			if term2 not in self.ancestors:
				return {}
			gene2anc = self.ancestors[term2]
		else:
			gene2anc = {}
			for i in term2:
				if i not in self.ancestors:
					continue
				gene2anc = self.int_merge_sets(gene2anc, self.ancestors[i])
		#print(gene2anc
		ca = {}
		for i in gene1anc:
			ca[i] = None
		for i in gene2anc:
			ca[i] = None
		return ca
	#

	def _has_IC(self, term):
		if self.IC == None:
			return False
		if not term in self.IC:
			return False
		if self.IC[term] == None:
			return False
		return True
	#
#
