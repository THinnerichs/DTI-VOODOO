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
Plain annotation corpus files parsing utility.
	Plain format 1: object (eg. gene) ID - Term ID
	Plain format 1: Term ID - object (eg. gene) ID
'''

import sys
import os
# import GeneOntology

MANYASSPERROW = 'multiple'
TERMFIRST = 'term first'
SEPARATOR = 'separator'
COMMENT = 'comment'

class PlainAnnotationCorpus(object):
	
#------------------------------------------------------
# Inizialization routines
#------------------------------------------------------

	def __init__(self, ac, parameters=None):
		self.ac = ac
		if self.ac == None:
			raise Exception
		self.parameters = parameters
		self.int_interpretParameters()

	def int_interpretParameters(self):
		self.int_obj_first = True
		self.int_one_association_per_row = True
		self.int_separator = '\t'
		self.int_comment = '#'

		if self.parameters == None:
			return
		if len(self.parameters) > 0:
			if MANYASSPERROW in self.parameters:
				if self.parameters[MANYASSPERROW]:
					self.int_one_association_per_row = False
				else:
					self.int_one_association_per_row = True
			if TERMFIRST in self.parameters:
				if self.parameters[TERMFIRST]:
					self.int_obj_first = False
				else:
					self.int_obj_first = True
			if SEPARATOR in self.parameters:
				self.int_separator = self.parameters[SEPARATOR]
			if COMMENT in self.parameters:
				self.int_comment = self.parameters[COMMENT]

#------------------------------------------------------
# Parsing routine
#------------------------------------------------------
	def setFields(self):
		# there are no additional fields in plain annotation files,
		pass

	def isOk(self):

		temp_term = self.temp_term
		if self.ac._exclude_roots:
			if temp_term in self.ac.go.roots:
				return False
		if not temp_term in self.ac.go.nodes:
			#print(str(self.temp_term) + " not found in GO.")
			return False
		# if not temp_term in self.ac.go.nodes:
			#print(str(self.temp_term) + " is obsolete.")
			#self.obso[self.temp_term] = None
			# return False
		return True

	def parse(self, fname):
		#self.obso = {}
		self.setFields()
		
		if type(fname) == unicode:
			fname = str(fname)
		if type(fname) is str:
			fn,fe = os.path.splitext(fname)
			if fe == '.gz':
				stream = gzip.open(fname, 'rb')
			else:
				stream = open(fname, 'r')
		else:
			stream = fname
		#filenum = rowcount(fname);
		#if filenum > self.SHOW_PROCESS_THRESHOLD:
			#print(fname + " has " + str(filenum) + " lines."
			#self.SHOW_PROCESS = True;
		#else:
		#self.SHOW_PROCESS = False;
		lines_counter = 0
		#ignored = 0
		for line in stream:
			line = line.rstrip('\n')
			line = line.rstrip('\r')
			if line[0] == self.int_comment:
				continue
			line = line.split(self.int_separator)
			if len(line) == 0:
				continue
			if len(line) < 2:
				#print("Strange line: " + str(line))
				continue
			temp_to_add = []
			if self.int_one_association_per_row:
				if self.int_obj_first:
					obj_id = line[0]
					term = line[1]
				else:
					obj_id = line[1]
					term = line[0]
				temp_to_add.append((obj_id, term))
			else:
				obj_id = line[0]
				for i in range(1,len(line)):
					term = line[i]
					temp_to_add.append((obj_id, term))

			for i in temp_to_add:
				obj_id = i[0]
				term = i[1]
				self.temp_term = term

				if obj_id not in self.ac.obj_set:
					self.ac.obj_set[obj_id] = {}

				self.temp_term_ar = self.ac.go.id2node(term, alt_check = True)
				if not type(self.temp_term_ar) == list:
					self.temp_term_ar = [self.temp_term_ar]

				for k in self.temp_term_ar:
					self.temp_term = k
					term = k
					# print self.temp_term
					if not self.isOk():
						continue

					if term not in self.ac.term_set:
						self.ac.term_set[term] = {}
					#### Build up annotations set
					if obj_id not in self.ac.annotations:
						self.ac.annotations[obj_id] = {}
					if term not in self.ac.annotations[obj_id]:
						self.ac.annotations[obj_id][term] = {}
					#### Build up reverse annotations set
					if term not in self.ac.reverse_annotations:
						self.ac.reverse_annotations[term] = {}
					if obj_id not in self.ac.reverse_annotations[term]:
						self.ac.reverse_annotations[term][obj_id] = {}
			lines_counter += 1
			#if self.ac.SHOW_PROCESS and (lines_counter%(filenum/20)==0):
				#print("Lines processed: " + str(lines_counter) + " on " + str(filenum) + " (" + str(int(100*float(lines_counter)/float(filenum))) + "%)")
		#print "Ignored: " + str(ignored)
		#print "Obso: " + str(len(self.obso))
		if type(fname) is str:
			stream.close()
		return True
