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
@desc CellOntology class handles CellOntology

Cell Ontology class
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

# def co_name2id(code):
# 	return int(code[3:])

# def co_id2name(code):
# 	# assumption: CO terms are 3 + 7 characters long.
# 	return "CL:" + '0'*(7 - len(str(code))) + str(code)
	
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# GeneOntology class
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

class CellOntology(Ontology.Ontology):
	
	# CO_root_str = "CL:0000000"
	# CO_root = CellOntology._name2id(CO_root_str)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# public functions and variables that should be used 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

	alt_ids = None
	obsolete_ids = None

	@staticmethod
	def _node2id(code):
		return "CL:" + '0'*(7 - len(str(code))) + str(code)
	#

	@staticmethod
	def _id2node(code, strict=True):
		if code.startswith('CL:'):
			try:
				return int(code[3:])
			except Exception:
				return None
		if strict:
			return None
		print code
		return Ontology.Ontology._id2node(code, strict)
	#

	# def name2id(self, codes, alt_check = True):
	# 	nid = None
	# 	if codes == None:
	# 		return None
	# 	if type(codes) is str:
	# 		nid = CellOntology._name2id(codes)
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
	# 				tnid = CellOntology._name2id(i)
	# 			else:
	# 				tnid = i
	# 			if alt_check:
	# 				if tnid in self.alt_ids:
	# 					tnid = self.alt_ids[tnid]
	# 			nid.append(tnid)
	# 	return nid

	# def id2name(self, codes, alt_check = False):
	# 	if alt_check:
	# 		print "id2name - alt_check not yet implemented."
	# 	sid = None
	# 	if codes == None:
	# 		return None
	# 	if type(codes) is int:
	# 		sid = CellOntology._id2name(codes)
	# 	elif type(codes) is str:
	# 		sid = codes
	# 	elif type(codes) is dict or type(codes) is list:
	# 		sid= []
	# 		for i in codes:
	# 			if type(i) is int:
	# 				tnid = CellOntology._id2name(i)
	# 			else:
	# 				tnid = i
	# 			sid.append(tnid)
	# 	return sid
	# #

	def __init__(self, terms, edges, parameters):
		terms['namespace'] = {} # impose this if current Ontology is faulty
		if parameters == None:
			parameters = {}
		if not 'ignore' in parameters:
			parameters['ignore'] = {'lacks_plasma_membrane_part':True} #, 'develops_from':True}
		Ontology.Ontology.__init__(self, terms = terms, edges = edges, parameters = parameters)
		self.name = 'CellOntology'
	#
#
