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

import os

try:
    from setuptools import setup # improved package
except ImportError:
    from distutils.core import setup # standard package

## Package path
pkg_path = os.path.dirname(__file__)

## Get package description
README = os.path.join(pkg_path, 'README')
lines = open(README).readlines()
description = ''.join(lines[:3])
long_description = ''.join(lines[:4])

## Get package Version
vh = open('fastsemsim/version','r')
lines = vh.readlines()
version = lines[0].rstrip('\n').rstrip('\r')
vh.close()

setup(
    name='fastsemsim',
    version=version,
    url='https://sites.google.com/site/fastsemsim/',
    description=description,
    long_description=long_description,
    keywords='semantic similarity, Ontology, Gene Ontology',
    author='Marco Mina',
    author_email='marco.mina.85@gmail.com',
    license='GNU GPL version 3',
    download_url = 'https://sourceforge.net/projects/fastsemsim/files/',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
    ],
    package_dir={'fastsemsim':'fastsemsim'},
    package_data={'fastsemsim.data':['dataset*','ACs/*','Os/*'], 'fastsemsim':['version'], 'fastsemsim.examples':['data/*','*.sh', '*.txt']},
    packages=['fastsemsim', 'fastsemsim.Ontology', 'fastsemsim.SemSim', 'fastsemsim.data', 'fastsemsim.examples'],
    install_requires=[
          'pandas',
      ],
    requires=[
          'pandas',
      ],
    scripts=[
        # 'fastsemsim/bin/fastsemsim.sh',
        # 'fastsemsim/bin/fastsemsim',
    ],
    entry_points = {
        'console_scripts': ['fastsemsim = fastsemsim.fastsemsim_cmdline:start']
    },
    include_package_data=True,
)
