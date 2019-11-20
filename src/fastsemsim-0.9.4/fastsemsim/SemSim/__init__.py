#!/usr/bin/env python
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
This module contains the classes for the evaluation of the Semantic Similarity.
Please refer to the single classes for details on the implemented measures.
'''

__author__="Marco Mina"
__email__="marco.mina.85@gmail.com"

# __all__ = ["SemSimUtils", "ObjSemSim", "SemSimMeasures", "TermSemSim", "MixSemSim", "ObjSetSemSim"]

from .ResnikSemSim import *
from .LinSemSim import *
from .JiangConrathSemSim import *
from .SimGICSemSim import *
from .SimUISemSim import *
from .SimICSemSim import *
from .SimRelSemSim import *
from .avgSemSim import *
from .maxSemSim import *
from .BMASemSim import *
from .DiceSemSim import *
from .SimTOSemSim import *
from .SimNTOSemSim import *
from .JaccardSemSim import *
from .CzekanowskiDiceSemSim import *
from .CosineSemSim import *
from .GSESAMESemSim import *
from .SimICNDSemSim import *
from .SimICNPSemSim import *
# from .ObjSemSim import *
# from .TermSemSim import *

'''
Struct term_SemSim_measures.
Contains a list of all available term SS measures.
It is built as a dictionary. SS measure names are used as keys. Each entry is a tuple with the following structure:
(Class Name, is Pairwise flag, )
'''
term_SemSim_measures = {
# present in version 0.6
	'Resnik' : (ResnikSemSim, True),
	'SimGIC': (SimGICSemSim, False),
	'SimUI': (SimUISemSim, False),
	'Lin' :(LinSemSim, True),
	'Jiang-Conrath' :(JiangConrathSemSim, True),
	'SimRel' :(SimRelSemSim, True),
	'SimIC' :(SimICSemSim, True),
	'Dice' :(DiceSemSim, False),
	'TO' :(SimTOSemSim, False),
	'NTO' :(SimNTOSemSim, False),
	'Jaccard' :(JaccardSemSim, False),
	'Czekanowski-Dice' :(CzekanowskiDiceSemSim, False),
	'Cosine' :(CosineSemSim, False),
	'G-SESAME' :(GSESAMESemSim, True),
	'SimICND' :(ICNDSemSim, True),
	'SimICNP' :(ICNPSemSim, True),
}


'''
Struct mix_strategies.
Contains a list of all available mixing strategies
It is built as a dictionary. Mixing strategy names are used as keys. Each entry is a tuple with the following structure:
(class pointer, )
'''
mix_strategies = {
	'max':(maxSemSim, ),
	'BMA':(BMASemSim, ),
	'avg':(avgSemSim, ),
}

# '''
# Struct obj_SemSim_measures.
# Contains a list of all available obj sem sim measures
# It is built as a dictionary. Mixing strategy names are used as keys. Each entry is a tuple with the following structure:
# (class pointer, )
# '''
# obj_SemSim_measures = {
# 	'obj':(ObjSemSim)
# }

	#-#-#-#-#-#-#-#-#-#-#-#-#-#
	# select Term Sem Sim     #
	#-#-#-#-#-#-#-#-#-#-#-#-#-#

	# the function selectTermSemSim helps retrieving the proper class implementing a given Term Sem Sim.
	# It takes in input the name (str) of the Term Sem Similarity
	# It returns the class to be used. Just call the class constructor to instantiante an object.

def select_term_SemSim(tss_name):
	if not tss_name in term_SemSim_measures:
		raise "Semantic Similarity measure not available."
		return None
	else:
		return term_SemSim_measures[tss_name][0]
#

def select_mix_SemSim(mix_name):
	if not mix_name in mix_strategies:
		raise "Semantic Similarity measure not available."
		return None
	else:
		return mix_strategies[mix_name][0]
#

# def select_obj_SemSim(oss_name='obj'):
# 	if not oss_name in obj_SemSim_measures:
# 		raise "Semantic Similarity measure not available."
# 		return None
# 	else:
# 		return obj_SemSim_measures[oss_name][0]
# #

# TermSemSim

# ss = ObjSetSemSim(ontology, ac, s, mix, ssutil, do_log = False)
# ss = ObjSemSim(ontology, ac, s, mix, ssutil, do_log = False) 
# ss = SetSemSim(ontology, ac, s, mix, ssutil, do_log = False)
# tss_class = SemSim.select_term_SemSim(s)
# ss = tss_class(ontology, ac, ssutil, do_log=False)

