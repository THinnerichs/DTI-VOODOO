# -*- coding: iso-8859-1 -*-

# Copyright 2011-2013 Marco Mina. All rights reserved.

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

# @mail marco.mina.85@gmail.com
# @version 2.0
# @desc GeneOntology class: extend Ontology to handle the GeneOntology

"""
Gene Ontology class
"""

# Out-of-class variables:
# 	types of relationships: is_a, part_of, regulates, ...
# Out-of-class functions:
# 	go_name2id
# 	go_id2name

# import types
# import os
# from xml.sax import make_parser
# from xml.sax.handler import ContentHandler
# import gzip
import Ontology

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# constants and macro
# assume the GO ids are in the standard format "GO:" + 7 digit number
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# IS_A = 0
# PART_OF = 1
# REGULATES = 2
# POS_REG = 3
# NEG_REG = 4
# HAS_PART = 5

# def go_name2id(code):
# 	return int(code[3:])

# def go_id2name(code):
# 	# assumption: GO terms are 3 + 7 characters long.
# 	return "GO:" + '0'*(7 - len(str(code))) + str(code)
	
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# GeneOntology class
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

class GeneOntology(Ontology.Ontology):
	
	# BP_root_str = "GO:0008150"
	# MF_root_str = "GO:0003674"
	# CC_root_str = "GO:0005575"

	# BP_root = go_name2id(BP_root_str)
	# MF_root = go_name2id(MF_root_str)
	# CC_root = go_name2id(CC_root_str)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# public functions and variables that should be used 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

	@staticmethod
	def _node2id(code):
		return "GO:" + '0'*(7 - len(str(code))) + str(code)
	#

	@staticmethod
	def _id2node(code, strict=True):
		if code.startswith('GO:'):
			try:
				return int(code[3:])
			except Exception:
				return None
		if strict:
			return None
		return Ontology.Ontology._id2node(code, strict)
	#


	# def name2id(self, codes, alt_check = True):
	# 	nid = None
	# 	if codes == None:
	# 		return nid
	# 	if type(codes) is str:
	# 		# nid = go_name2id(codes)
	# 		nid = GeneOntology._name2id(codes,strict=True)
	# 		nid = self.name2id(nid, alt_check)
	# 	elif type(codes) is int:
	# 		nid = codes
	# 		if alt_check:
	# 			if nid in self.alt_id:
	# 				nid = self.alt_id[nid]
	# 	elif type(codes) is dict or type(codes) is list:
	# 		nid = []
	# 		for i in codes:
	# 			if type(i) is str:
	# 				# tnid = go_name2id(i)
	# 				tnid = GeneOntology._name2id(i,strict=True)
	# 			else:
	# 				tnid = i
	# 			if alt_check:
	# 				if tnid in self.alt_id:
	# 					tnid = self.alt_id[tnid]
	# 			nid.append(tnid)
	# 	return nid

	# def id2name(self, codes, alt_check = False):
	# 	if alt_check:
	# 		print "id2name - alt_check not yet implemented."
	# 	sid = None
	# 	if codes == None:
	# 		return sid
	# 	if type(codes) is int:
	# 		sid = GeneOntology._id2name(codes)
	# 	elif type(codes) is str:
	# 		sid = codes
	# 	elif type(codes) is dict or type(codes) is list:
	# 		sid= []
	# 		for i in codes:
	# 			if type(i) is int:
	# 				tnid = GeneOntology._id2name(i)
	# 			else:
	# 				tnid = i
	# 			sid.append(tnid)
	# 	return sid
	# #

	def __init__(self, terms, edges, parameters):
		# terms['namespace'] = {} # impose this if current Ontology is faulty
		if parameters == None:
			parameters = {}
		if not 'ignore' in parameters:
			parameters['ignore'] = {'has_part':True, 'happens_during':True, 'occurs_in':True}
		Ontology.Ontology.__init__(self, terms = terms, edges = edges, parameters = parameters)
		self.name = 'GeneOntology'
	#
#
