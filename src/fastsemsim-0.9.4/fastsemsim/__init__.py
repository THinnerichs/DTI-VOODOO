#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
'''
Copyright 2011 Marco Mina. All rights reserved.

This file is part of fastSemSim

fastSemSim is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fastSemSim is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with fastSemSim.  If not, see <http://www.gnu.org/licenses/>.
'''

__author__="Marco Mina"
__email__="marco.mina.85@gmail.com"

import os
vh = open(os.path.join(os.path.dirname(os.path.abspath(__file__)).replace("\\", "/"),'version'),'r')
lines = vh.readlines()
vh.close()
__version__ = lines[0].rstrip('\n').rstrip('\r')

from . import Ontology