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


"""
@mail marco.mina.85@gmail.com
@version 2.0
@desc DiseaseOntology class handles DiseaseOntology

Disease Ontology class
"""
# import types
# import os
# from xml.sax import make_parser
# from xml.sax.handler import ContentHandler
# import gzip
import Ontology

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# constants and macro
# assume the DO ids are in the standard format "DOID:" + 7 digit number
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

# IS_A = 0
# PART_OF = 1
# REGULATES = 2
# POS_REG = 3
# NEG_REG = 4
# HAS_PART = 5

# def do_id2node(code):
# 	return int(code[5:])

# def do_node2id(code):
# 	# assumption: GO terms are 5 + 7 characters long.
# 	return "DOID:" + '0'*(7 - len(str(code))) + str(code)
	
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# GeneOntology class
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

class DiseaseOntology(Ontology.Ontology):

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# public functions and variables that should be used 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

	alt_ids = None
	obsolete_ids = None

	@staticmethod
	def _node2id(code):
		return "DOID:" + '0'*(7 - len(str(code))) + str(code)
	#

	@staticmethod
	def _id2node(code, strict=True):
		if code.startswith('DOID:'):
			try:
				return int(code[5:])
			except Exception:
				return None
		if strict:
			return None
		return Ontology.Ontology._id2node(code)
	#

	# DO_root_str = "DOID:0000004"
	# DO_root = _id2node(DO_root_str)

	# def name2id(self, codes, alt_check = True):
	# 	nid = None
	# 	if codes == None:
	# 		return None
	# 	if type(codes) is str:
	# 		nid = DiseaseOntology._id2node(codes)
	# 		nid = self.name2id(nid, alt_check)
	# 	elif type(codes) is int:
	# 		nid = codes
	# 		if alt_check:
	# 			if nid in self.alt_ids:
	# 				nid = self.alt_ids[nid]
	# 	elif type(codes) is dict or type(codes) is list:
	# 		nid = []
	# 		for i in codes:
	# 			if type(i) is str:
	# 				tnid = DiseaseOntology._id2node(i)
	# 			else:
	# 				tnid = i
	# 			if alt_check:
	# 				if tnid in self.alt_ids:
	# 					tnid = self.alt_ids[tnid]
	# 			nid.append(tnid)
	# 	return nid

	# def id2name(self, codes, alt_check = False):
	# 	if codes == None:
	# 		return None
	# 	if alt_check:
	# 		print "id2name - alt_check not yet implemented."
	# 	sid = None
	# 	if type(codes) is int:
	# 		sid = DiseaseOntology._node2id(codes)
	# 	elif type(codes) is str:
	# 		sid = codes
	# 	elif type(codes) is dict or type(codes) is list:
	# 		sid= []
	# 		for i in codes:
	# 			if type(i) is int:
	# 				tnid = DiseaseOntology._node2id(i)
	# 			else:
	# 				tnid = i
	# 			sid.append(tnid)
	# 	return sid
	# #

	def __init__(self, terms, edges, parameters):
		Ontology.Ontology.__init__(self, terms = terms, edges = edges, parameters = parameters)
		self.name = 'DiseaseOntology'
	#
#
