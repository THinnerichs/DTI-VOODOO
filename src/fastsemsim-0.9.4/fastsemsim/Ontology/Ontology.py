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
# @desc Ontology is the basic class representing ontologies representable with multirooted connected DAGs.


'''
Supported ontologies are those representable as multirooted DAGs. It is not required DAGs to be disconnected, but 'inter-DAG' edges are required to be specified. Class Ontology provides a function is_consistent that checks whether this contraints is satisfied. Inconsistent DAGs are NOT currently usable.

Different datastructures can be used to represent ontologies. The section Variables lists a set of different alternatives. Currently Ontology is tuned for using a parent-children and a node-edge representation.

Superclasses can extend the basic datastructure with additional layers of information. 
'''


# 	_use_parent_children_ = True
# 	_use_node_edges_ = True

# 	# available_data_structs = ('dicts', '')
# 	roots = {} # ids of roots 
# 	# data_structs = []

# # data struct 1: nodes and edges
# 	nodes = {} # each node is a list of edge ids
# 	edges = {} # dict. Each key is an attrib of the edge

# 	# variant 1(1a)
# 	# edges['inter'] = [] # list. Whether the edge is inter-category or not
# 	# variant 1(1b)
# 	edges['inter'] = {} # list. Whether the edge is inter-category or not

# 	# variant 1(2a)
# 	# edges['parent'] = [] # list. Parent node id
# 	# edges['child'] = [] # list. Child node id
# 	# variant 1(2b)
# 	edges['nodes'] = [] # list. (parent, child) node ids

# 	# additional 1(3a)
# 	# edges['attrib'] = [] # list. Additional attributes
	
# # data struct 2: nodes and edges, edges handled as table using numpy
# 	# import numpy as np
# 	# nodes = {} # each node is a list of edge ids
# 	# edges = np.zeros(shape=(0,3)) # numpy array. column 1: parent, column 2: child, column 3: inter-category, column 4: attribs
# 	# FIX PROBLEM OF ONLY INT IN MATRIX

# # data struct 3: parents and children
# 	#variant 3(1a)
# 	parents = {} # dict. key: node id, value: dict of parents {}
# 	children = {} # dict. key: node id, value: dict of children {}
# 	#variant 3(1b)
# 	# parents = [] # list. pos: node id, value: dict of parents {}
# 	# children = [] # list. pos: node id, value: dict of children {}

# # ? node names are int and sequential? We assume they are NOT. It implies variant 3(1b) requires a mapping structure as well. It implies a potential overhead?

# # data struct 3: ??
# 	# alt_id = None
# 	# obsolete_ids = None



import pandas as pd
import numpy as np

	# id_tag = "id"
	# name_tag = "name"
	# def_tag = "def"
	# alt_id_tag = "alt_id"
	# replaced_by_tag = "replaced_by"
	# consider_tag = "consider"

	# namespace_tag = "namespace"
	# is_obsolete_tag = "is_obsolete"
	# is_a_tag = "is_a"
	# relationship_tag = "relationship"

class Ontology(object):
	'''
	Base class for the representation of an ontology. It currently supports any multi-rooted DAG (Directed Acyclic Graph).
	'''
	debug = False
	gen_error = False

	@staticmethod
	def _node2id(code, strict=True):
		''' Converts an ontology term node id into its ontology form
		'''
		# return "generic:" + '0'*(7 - len(str(code))) + str(code)
		return str(code) # watch out for Unicode strings
	#

	@staticmethod
	def _id2node(code, strict=True):
		return str(code) # watch out for Unicode strings
	#

	def id2node(self, codes, alt_check = True):
			nid = None
			if codes == None:
				return nid
			if isinstance(codes,str):
				# nid = go_name2id(codes)
				nid = self._id2node(codes,strict=True)
				nid = self.id2node(nid, alt_check)
			elif isinstance(codes,int):
				nid = codes
				if alt_check:
					if nid in self.alt_id:
						nid = self.alt_id[nid]
			elif isinstance(codes,(dict,list,tuple)):
				nid = []
				for i in codes:
					if type(i) is str:
						# tnid = go_name2id(i)
						tnid = self._id2node(i,strict=True)
					else:
						tnid = i
					if alt_check:
						if tnid in self.alt_id:
							tnid = self.alt_id[tnid]
					nid.append(tnid)
			return nid
	#

	def node2id(self, codes, alt_check = False):
		if alt_check:
			if self.debug:
				print "id2name - alt_check not yet implemented."
		sid = None
		if codes == None:
			return sid
		if type(codes) is dict or type(codes) is list:
			sid= []
			for i in codes:
				# if type(i) is int:
				tnid = self._node2id(i)
				# else:
					# tnid = i
				sid.append(tnid)
		# if type(codes) is int:
			# sid = Ontology._node2id(codes)
		else:
		 # type(codes) is str:
			# sid = codes
			sid = self._node2id(codes)
		return sid
	#

	def name2node(self, codes):
		match = []
		for i in self.node_attributes:
			if 'name' in self.node_attributes[i]:
				for j in self.node_attributes[i]['name']:
					if j == codes:
						match.append(i)
						return i
		return  None
	#


	def node_number(self):
		return len(self.nodes)

	def edge_number(self):
		return(self.edges.shape[0])

	def __init__(self, terms, edges, parameters=None):
		'''Initialization.

		:Parameters:
			terms : dict
				Data of the ontology terms
			edges : list
				List of ontological relationships
			parameters: dict [Default = None]
				parameters affecting the construction of the ontology
		'''
		
		self.name = 'Ontology'
		self.nodes = {}
		self.alt_id = {}
		self.obsolete = {}
		self.namespace = {}
		self.node_attributes = {} # ['subset', 'comment', 'xref', 'synonym', 'intersection_of'] # Ignored attributes
		self.parents = {}
		self.children = {}

		self.ignore = None
		if not parameters == None and 'ignore' in parameters:
			self.ignore = parameters['ignore']

		self._add_nodes(terms)
		self._add_edges(edges)
		self.det_roots()

		# Check consistency of alt_names
		if not self.alt_id == None:
			for i in self.alt_id:
				if i in self.nodes:
					# if self.gen_error:
					print "Warning: Inconsistent redefinition of valid term " + str(self._node2id(i)) + " as an alternative of " + str(self._node2id(self.alt_id[i]))
					# raise Exception # it means there are inconsistencies!
					# if self.debug:
					# raise Exception # it means there are inconsistencies!
		#
	#


	def _add_node(self,n):
		if n not in self.nodes:
			self.nodes[n] = []
		else:
			raise(Exception, 'Node ' + str(n) + ' already in the ontology')
	#

	def _add_nodes(self, terms): #input: disctionary of terms to be added.
		for i in terms['id']:
			if i == None:
				continue
			iid = self._id2node(i, strict = True)
			if iid == None: # not in ontology. Do not save data for now
				continue
			
			# process alt_id and replaced_by
			if i in terms['alt_id']:
				for j in terms['alt_id'][i]:
					jid = self._id2node(j, strict = True)
					if jid == None: 
						pass
					else:
						if jid in self.nodes:
							print "Warning: ignoring inconsistent alt id: " + j
							# raise Exception
						else:
							if jid not in self.alt_id:
								self.alt_id[jid] = []
							self.alt_id[jid].append(iid)
			#
			if i in terms['replaced_by']:
				if iid in self.nodes:
					# print 
					raise (Exception, "Inconsistent replaced_by id: " + i)
				for j in terms['replaced_by'][i]:
					jid = self._id2node(j, strict = True)
					if jid == None: 
						pass
					else:
						if iid not in self.alt_id:
							self.alt_id[iid] = []
						self.alt_id[iid].append(jid)
			#

			if i in terms['is_obsolete']: # do not save data for now
				# idata = {}
				if iid not in self.obsolete:
					self.obsolete[iid] = terms['is_obsolete'][i]
			else:
				if iid not in self.nodes:
					self.nodes[iid] = []
				else:
					print "Duplicated term: " + i
					raise Exception

				if i in terms['namespace']:
					self.namespace[iid] = terms['namespace'][i]

				idata = {}
				for j in { 'name', 'def'}:
					if i in terms[j]:
						idata[j] = terms[j][i]
				self.node_attributes[iid] = idata
			#
		#
	#

	def _add_edges(self, edges): #input: list of edges to be added. Each item of the list is a list [child, parent, type]
		if not self.ignore == None:
			newedges = []
			for i in range(0,len(edges)):
				# print edges[i][2]
				if not edges[i][2] in self.ignore:
					newedges.append(edges[i])
			edges = newedges
		#
		for i in range(0,len(edges)):
			childid = self._id2node(edges[i][0], strict = True)
			parentid = self._id2node(edges[i][1], strict = True)
			inner = True
			intra = True
			if childid == None :
				inner = False
				childid = edges[i][0]
			if parentid == None :
				inner = False
				parentid = edges[i][1]
			if inner and (childid in self.namespace) and (parentid in self.namespace) and (not self.namespace[parentid] == self.namespace[childid]):
				intra = False
			elif not inner:
				intra = False
			edges[i][0] = childid
			edges[i][1] = parentid
			edges[i].append(inner)
			edges[i].append(intra)

			if inner:
				if parentid not in self.nodes:
					if parentid in self.obsolete:
						print "Skipping relationship " + str(edges[i]) + ": parent node is obsolete: " + str(parentid)
						continue
					else:
						print "Adding ghost node " + str(parentid)
						self._add_node(parentid)
				if childid not in self.nodes:
					if childid in self.obsolete:
						print "Skipping relationship " + str(edges[i]) + ": child node is obsolete: " + str(childid)
						continue
					else:
						print "Adding ghost node " + str(childid)
						self._add_node(childid)
				self.nodes[parentid].append(i)
				if parentid != childid:
					self.nodes[childid].append(i)

		self.edges = pd.DataFrame(edges, columns=['child','parent', 'type', 'inner', 'intra'])
	#

	def _s1_to_s2(self): # given data struct 1, generates data struct 2
		self.parents = {}
		self.children = {}
		for j in self.nodes:
			self.parents[j] = {}
			self.children[j] = {}
		for count, row in self.edges.iterrows():
			if row['child'] in self.nodes and row['parent'] in self.nodes:
				self.parents[row['child']][row['parent']] = count
				self.children[row['parent']][row['child']] = count
	#

	def det_roots(self): # determines roots (nodes without parents)
		# using data struct 1
		if len(self.parents) == 0:
			self._s1_to_s2()
		self.roots = {}
		for i in self.parents:
			if len(self.parents[i]) == 0:
				self.roots[i] = None
	#

	def is_consistent(self): # verifies the Ontology is consistent
		self.det_roots()
		consistent = True
		inc = {}
		temp = {}
		tempq = []
		for i in self.roots:
			tempq.append(i)
			temp[i] = i
		while len(tempq) > 0:
			current = tempq.pop()
			current_cat = temp[current]
			current_children = self.children[current]
			# print current_children
			for i in current_children:
				if current_children[i] in self.edges['inter'] and self.edges['inter'][current_children[i]]:
					# print('CASO INTER')
					continue
				if not i in temp:
					temp[i] = current_cat
					# print str(i) + ' gets ' + str(current_cat)
				else:
					if not temp[i] == current_cat:
						consistent = False
						# print "happens"
						# print str(i) + "-" + str(current)
						break
					pass
				tempq.append(i)
		# print len(temp)
		# print self.node_number()
		if len(temp) < self.node_number():
			consistent = False
		return consistent
	#

	# verify if a term is valid (not obsolete and inside the ontology). Single terms are required. Does not work with lists/dicts.
	def is_valid(self, term):
		if term in self.nodes:
			return True
		return False
#
