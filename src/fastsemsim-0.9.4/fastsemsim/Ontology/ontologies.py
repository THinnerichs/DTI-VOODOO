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
# @desc Parser module: load/store ontologies


"""
Set of functions to parse and handle ontologies.
"""

import types
import os
from xml.sax import make_parser
from xml.sax.handler import ContentHandler
import gzip
import Ontology
import DiseaseOntology
import CellOntology
import GeneOntology
import FFOntology
from fastsemsim import data

'''
Struct ontologies.
Contains a list of all available ontologies.
It is built as a dictionary. ontology names are used as keys. Each entry is a tuple with the following structure:
(Class Name, )
'''
ontologies = {
	'GeneOntology' : (GeneOntology.GeneOntology, ),
	'DiseaseOntology' : (DiseaseOntology.DiseaseOntology, ),
	'CellOntology' : (CellOntology.CellOntology, ),
	'FFOntology' : (FFOntology.FFOntology, ),
	'GenericOntology' : (Ontology.Ontology, ),
	'Ontology' : (Ontology.Ontology, )
}

ontology_source_types = ['obo-xml', 'obo']

def parse(source = None, source_type = 'obo', ontology_type = 'GeneOntology', parameters={}):
	return load(source, source_type, ontology_type, parameters)
#

def load(source = None, source_type = 'obo', ontology_type = 'GeneOntology', parameters={}):
	ontology = None
	# namespace = None

	if source == None:
		# program_dir = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
		# print "ontologies.py: " + program_dir
		builtin_dataset = data.dataset.Dataset()
		selected_source = builtin_dataset.get_default_ontology(ontology_type)
		# print selected_source
		if selected_source is None:
			return None
		source = selected_source['file']
		source_type = selected_source['filetype']

	# generate source file handle
	if type(source) == unicode:
		source = str(source)
	if type(source) == str:
		fn,fe = os.path.splitext(source)
		if fe == '.gz':
			source_handle = gzip.open(source, 'rb')
		else:
			source_handle = open(source, 'rU')
	else: # assume that the passed object is a file stream
		source_handle = source

	# select proper input parser
	if 'ontology_type' in parameters:
		ontology_type = parameters['ontology_type']
	if ontology_type in ontologies:
		ontology_class = ontologies[ontology_type][0]
	else:
		raise Exception
	if 'source_type' in parameters:
		source_type = parameters['source_type']

	# parse data
	if source_type == 'obo-xml':
		parser = make_parser()
		handler = OboXmlParser(ontology_class, parameters)
		parser.setContentHandler(handler)
		# print "A"
		parser.parse(source_handle)
		# print "B"
	elif source_type == 'obo':
		handler = OboParser(ontology_class, parameters)
		handler.parse(source_handle)
		# namespace = handler.namespace
	else:
		# print "GeneOntology load: Unknown file format: " + str(parameters['type'])
		raise Exception

	if type(source) == str: # if original source was a handle, close input file
		source_handle.close()

	# # postprocess data, if required
	# if 'ignore' in parameters and 'inter' in parameters['ignore'] and parameters['ignore']['inter']:
	# 	for i in range(0, len(handler.edges)):
	# 		(u,v,z) = handler.edges[i] 
	# 		if not namespace == None:
	# 			vn = None
	# 			un = None
	# 			if u in namespace:
	# 				un = namespace[u]
	# 			if v in namespace:
	# 				vn = namespace[v]
	# 			if not un == vn:
	# 				handler.edges[i] = None	
	
	# build ontology
	# print(handler.terms.keys())
	# print(handler.edges)
	# return(handler.terms, handler.edges)

	ontology = ontology_class(handler.terms, handler.edges, parameters)
	# print "LOADED"
	# print len(handler.terms['id'])
	# print len(ontology.nodes)
	# for i in handler.terms:
		# if not i in ontology.nodes:
			# print i
	# for i in ontology.nodes:
		# if not i in handler.terms['id']:
			# print i


	return ontology
#


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# OboParser: Class to parse obo files
# Parse input file
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

class OboParser:

	term_tag = '[Term]'
	typedef_tag = '[Typedef]'
	comment_tag = '!'
	delimiter_tag = ":"
	relationship_delimiter_tag = " "

	id_tag = "id"
	name_tag = "name"
	def_tag = "def"
	alt_id_tag = "alt_id"
	replaced_by_tag = "replaced_by"
	consider_tag = "consider"

	namespace_tag = "namespace"
	is_obsolete_tag = "is_obsolete"
	is_a_tag = "is_a"
	relationship_tag = "relationship"

	def __init__(self, ontology_class, parameters = {}):
		self.ontology_class = ontology_class

		# self.ignore_part_of = False
		# self.ignore_regulates = True ### NOTE regulates ignored by default!
		# self.ignore_has_part = True ### NOTE has part ignored by default!
		# self.ignore_is_a = False
		
		# if 'ignore' in parameters and type(parameters['ignore']) == dict:
		# 	ignore = parameters['ignore']
		# 	if 'part_of' in ignore:
		# 		self.ignore_part_of = bool(ignore['part_of'])
		# 	if 'regulates' in ignore:
		# 		self.ignore_regulates = bool(ignore['regulates'])
		# 	if 'has_part' in ignore:
		# 		self.ignore_has_part = bool(ignore['has_part'])
		# 	if 'is_a' in ignore:
		# 		self.ignore_is_a = bool(ignore['is_a'])
	
		# Initialize structures

		self.tags = {}
		self.tags[self.id_tag] = True
		self.tags[self.namespace_tag] = True
		self.tags[self.alt_id_tag] = True
		self.tags[self.name_tag] = True
		self.tags[self.def_tag] = True
		self.tags[self.replaced_by_tag] = True
		self.tags[self.is_obsolete_tag] = True
		self.tags[self.consider_tag] = False
		self.tags[self.is_a_tag] = True
		self.tags[self.relationship_tag] = True
		
		self.relationship_tags = {}
		self.relationship_tags[self.is_a_tag] = True

		self.terms = {} # ids of terms contained in the ontology
		self.terms[self.namespace_tag] = {}
		self.terms[self.id_tag] = {}
		self.terms[self.alt_id_tag] = {}
		self.terms[self.is_obsolete_tag] = {}
		self.terms[self.replaced_by_tag] = {}
		self.terms[self.name_tag] = {}
		self.terms[self.def_tag] = {}

		self.edges = [] # set of relationships between terms

	#


	'''
	Find next block of informations
	'''
	def find_next_term(self):
	    # read each line until it has a certain start, and then puts the start tag back
	    while True:
	        # pos = self.handle.tell()
	        line = self.handle.readline()
	        if not line:
	            break
	        if line.startswith(self.term_tag):
	            # self.handle.seek(pos)
	            return True
	    return False
	#

	'''
	Split tag and value line; strip trailing and ending spaces and return characters
	'''
	def split_tag(self, st):
		st = st.split(self.comment_tag, 1)[0].strip() # remove comments 
		st = st.split(self.delimiter_tag, 1)
		# print st
		for i in range(0,len(st)):
			st[i] = st[i].strip()

		return(st)
	#

	'''
	Split relationship tag; strip trailing and ending spaces and return characters
	'''
	def split_rel_tag(self, st):
		st = st.split(self.relationship_tag, 1)[0].strip() # remove comments 
		st = st.split(self.relationship_delimiter_tag, 1)
		# print st
		for i in range(0,len(st)):
			st[i] = st[i].strip()

		return(st)
	#
	
	'''
	Check if lines have a key:value structure
	'''
	def has_tag_structure(self, st):
		if st == None:
			return False
		if st.startswith(self.comment_tag):
			return False
		st2 = st.split(self.delimiter_tag, 1)
		if(len(st2)) < 2:
			return False
		return True
	#

	'''
	Parse the file
	'''
	def parse(self, _handle):
		self.handle = _handle
		while True:
			if not self.find_next_term(): # go to next [Term]
				break

			# collect lines from current block. Stop when find another block
			lines = []
			while 1:
				pos = self.handle.tell()
				line = self.handle.readline()
				if not line:
					break
				if line.startswith(self.typedef_tag) or line.startswith(self.term_tag):
					self.handle.seek(pos)
					break
				lines.append(line)

			# process info in captured block
			got_term = False
			got_rel = False

			temp_data = {}
			for i in self.terms.keys():
				temp_data[i] = None
			#

			temp_rel = []

			for line in lines:
				if not self.has_tag_structure(line):
					continue
				key, value = self.split_tag(line)

				if not key in self.tags:
					# print "New tag found: " + str(key)
					self.tags[key] = True
					self.terms[key] = {}
					temp_data[key] = None
				#
				if not self.tags[key]:
					# print "Ignoring tag: " + str(key)
					continue
				#

				if key == self.id_tag:
					if got_term:
						raise Exception
					term_id = value
					temp_data[self.id_tag] = value
					got_term = True
				#

				elif key == self.alt_id_tag:
					# term_alt_ids.append(value)
					if temp_data[self.alt_id_tag] == None:
						temp_data[self.alt_id_tag] = []
					temp_data[self.alt_id_tag].append(value)
				#

				elif key == self.replaced_by_tag:
					if temp_data[self.replaced_by_tag] == None:
						temp_data[self.replaced_by_tag] = []
					temp_data[self.replaced_by_tag].append(value)
					# term_replaced_ids.append(value)
				#			

				elif key == self.namespace_tag:
					temp_data[self.namespace_tag] = value

				elif key == self.is_obsolete_tag and value=="true":
					if not got_term:
						raise Exception
					if got_rel:
						# print "Obsolete with edges!"
						# raise Exception
						pass
					temp_data[self.is_obsolete_tag] = True
				#

				elif key == self.is_a_tag:
					# if not got_term:
						# raise Exception
					if temp_data[self.is_obsolete_tag]:
						# print "Obsolete with edges!"
						# raise Exception
						pass
					isa = value
					temp_rel.append( (isa, self.is_a_tag ) )
					got_rel = True
				#

				elif key == self.relationship_tag:
					# if not got_term:
						# pass
					if temp_data[self.is_obsolete_tag]:
						# print "Obsolete with edges!"
						# raise Exception
						pass
					cline = self.split_rel_tag(value)
					ctype = cline[0]
					cto = cline[1]
					if not ctype in self.relationship_tags:
						# print "Found new relationship: " + ctype
						self.relationship_tags[ctype] = True
					if self.relationship_tags[ctype]:
						temp_rel.append( (cto, ctype ) )
						got_rel = True
				#

				else: # it already passed the True/False check
					if temp_data[key]==None:
						temp_data[key] = []
					temp_data[key].append(value) # treat as unique value.

			# commit Term info
			if not got_term:
				print "Missing term - skipping block"
				# raise Exception
				continue
			if term_id in self.terms[self.id_tag]: # add term
				print "Duplicated term: " + str(term_id)
			for k in temp_data:
				if not temp_data[k] == None:
					self.terms[k][term_id] = temp_data[k]
			for i in temp_rel:
				self.edges.append( [ term_id, i[0], i[1] ] )

		#
	#
#


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# OboXmlParser: Class to parse GO obo-xml files
# Given an obo-xml file builds a GeneOntology object parsing it
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

class OboXmlParser(ContentHandler):
	
	term_tag = 'term'
	comment_tag = '!'
	delimiter_tag = ":"
	relationship_delimiter_tag = " "

	id_tag = "id"
	name_tag = "name"
	def_tag = "def"
	alt_id_tag = "alt_id"
	replaced_by_tag = "replaced_by"
	consider_tag = "consider"

	namespace_tag = "namespace"
	is_obsolete_tag = "is_obsolete"
	is_a_tag = "is_a"
	part_of_tag = "part_of"
	relationship_tag = "relationship"
	relationship_type_tag = "type"
	relationship_to_tag = "to"

	def __init__(self, ontology_class, parameters = {}):
		self.ontology_class = ontology_class # type of ontology to load

		self.tags = {}
		self.tags[self.id_tag] = True
		self.tags[self.namespace_tag] = True
		self.tags[self.alt_id_tag] = True
		self.tags[self.name_tag] = True
		self.tags[self.def_tag] = True
		self.tags[self.replaced_by_tag] = True
		self.tags[self.is_obsolete_tag] = True
		self.tags[self.consider_tag] = False
		self.tags[self.is_a_tag] = True
		self.tags[self.relationship_tag] = True
		
		self.relationship_tags = {}
		self.relationship_tags[self.is_a_tag] = True

		self.terms = {} # ids of terms contained in the ontology
		self.terms[self.namespace_tag] = {}
		self.terms[self.id_tag] = {}
		self.terms[self.alt_id_tag] = {}
		self.terms[self.is_obsolete_tag] = {}
		self.terms[self.replaced_by_tag] = {}
		self.terms[self.name_tag] = {}
		self.terms[self.def_tag] = {}

		self.edges = [] # set of relationships between terms

#

		self.isId, self.isIsA, self.isPartOf, self.isaltId, self.isRelationship = 0,0,0,0,0
		self.isObsolete, self.isReplacedBy, self.isConsider = 0,0,0
		self.isName, self.isDef = 0,0
		self.isRelationshipTo, self.isRelationshipType = 0,0
		self.isNamespace = 0
		self.inTerm = 0
		self.got_term = False
		self.got_rel = False
		self.temp_data = {}
		for i in self.terms.keys():
			self.temp_data[i] = None
		#
		self.temp_rel = []


	def startElement(self, name, attrs):
		if name == self.term_tag:
			self.inTerm = 1
		elif self.inTerm == 1:
			if name == self.id_tag:
				self.isId = 1
				self.id = ''
			elif name == self.is_a_tag:
				self.isIsA = 1
				self.isa = ''
			elif name == self.part_of_tag:
				self.isPartOf = 1
				self.partof = ''
			elif name == self.alt_id_tag:
				self.isaltId = 1
				self.curaltid = ''
			elif name == self.relationship_tag:
				self.isRelationship = 1
			elif name == self.relationship_type_tag:
				if self.isRelationship:
					self.isRelationshipType = 1
					self.parent_type = ''
			elif name == self.relationship_to_tag:
				if self.isRelationship:
					self.isRelationshipTo = 1
					self.parent = ''
			elif name == self.is_obsolete_tag:
				self.isObsolete = 1
			elif name == self.replaced_by_tag:
				self.isReplacedBy = 1
				self.currepid = ''
			elif name == self.consider_tag:
				self.isConsider = 1
			elif name == self.namespace_tag:
				self.isNamespace = 1
				self.namespace = ''
			elif name == self.name_tag:
				self.isName = 1
				self.name = ''
			elif name == self.def_tag:
				self.isDef = 1
				self.defi = ''
			else:
				# print "Unknown tag: " + name 
				pass
#


	def endElement(self, name):
		if self.inTerm == 1:

			if name == self.term_tag:
				if not self.got_term:
					raise Exception
				self.terms[self.id_tag][self.id] = {}

				for k in self.temp_data:
					if not self.temp_data[k] == None:
						self.terms[k][self.id] = self.temp_data[k]
				for i in self.temp_rel:
					self.edges.append( [ self.id, i[0], i[1] ] )
				
				self.isId, self.isIsA, self.isPartOf, self.isaltId, self.isRelationship = 0,0,0,0,0
				self.isObsolete, self.isReplacedBy, self.isConsider = 0,0,0
				self.isRelationshipTo, self.isRelationshipType = 0,0
				self.inTerm = 0
				self.got_term = False
				self.got_rel = False
				self.temp_data = {}
				for i in self.terms.keys():
					self.temp_data[i] = None
				#
				self.temp_rel = []
			#

			elif name == self.id_tag:
				self.isId = 0
				self.got_term = True
				self.temp_data[self.id_tag] = self.id
			#
			elif name == self.namespace_tag:
				self.isNamespace = 0
				self.temp_data[self.namespace_tag] = self.namespace
			#
			elif name == self.name_tag:
				self.isName = 0
				self.temp_data[self.name_tag] = self.name
			#
			elif name == self.def_tag:
				self.isDef = 0
				self.temp_data[self.def_tag] = self.defi
			#
			elif name == self.alt_id_tag:
				self.isaltId = 0
				if not self.got_term:
					raise Exception
				if self.temp_data[self.alt_id_tag] == None:
					self.temp_data[self.alt_id_tag] = []
				self.temp_data[self.alt_id_tag].append(self.curaltid)
			#
			elif name == self.is_obsolete_tag:
				self.isObsolete = 0
				self.temp_data[self.is_obsolete_tag] = True
			#
			elif name == self.replaced_by_tag:
				self.isReplacedBy = 0
				if not self.got_term:
					raise Exception
				if self.temp_data[self.replaced_by_tag] == None:
					self.temp_data[self.replaced_by_tag] = []
				self.temp_data[self.replaced_by_tag].append(self.currepid)
			#
			elif name == self.consider_tag:
				self.isConsider = 0
			#
			elif name == self.is_a_tag:
				self.isIsA = 0
				# if self.curobsolete:
					# raise Exception
				isaid = self.isa
				self.got_rel = True
				self.temp_rel.append( (isaid, self.is_a_tag ) )
			#
			elif name == self.part_of_tag:
				self.isPartOf = 0
				# if self.curobsolete:
					# raise Exception
				partofid = self.partof
				self.got_rel = True
				self.temp_rel.append( (partofid, self.part_of_tag ) )
			#
			elif name == self.relationship_tag:
				# if self.curobsolete:
					# raise Exception
				self.isRelationship = 0
			#
			elif name == self.relationship_type_tag:
				self.isRelationshipType = 0
			#
			elif name == self.relationship_to_tag:
				if self.isRelationshipTo:
					self.isRelationshipTo = 0
					if not self.got_term:
						raise Exception
					# if self.curobsolete:
						# raise Exception
					cto = self.parent
					self.got_rel = True
					self.temp_rel.append( (cto, self.parent_type ) )
			#
	#

	def characters(self, ch):
		if self.isId == 1:
			self.id += ch
		if self.isNamespace == 1:
			self.namespace += ch
		elif self.isIsA == 1:
			self.isa += ch
		elif self.isPartOf == 1:
			self.partof += ch
		elif self.isaltId == 1:
			self.curaltid += ch
		elif self.isRelationshipTo == 1:
			self.parent += ch
		elif self.isRelationshipType == 1:
			self.parent_type += ch
		elif self.isObsolete == 1:
			self.curobsolete = bool(ch)
		elif self.isReplacedBy == 1:
			self.currepid += ch
		elif self.isName == 1:
			self.name += ch
		elif self.isDef == 1:
			self.defi += ch
		elif self.isConsider == 1:
			pass
	#
#
