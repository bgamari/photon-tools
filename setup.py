#!/usr/bin/python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
        Extension('photon_tools.bin_photons', ['photon_tools/bin_photons.pyx'], include_dirs=['.']),
        Extension('photon_tools.timetag_parse', ['photon_tools/timetag_parse.pyx'], include_dirs=['.']),
        Extension('photon_tools.filter_photons', ['photon_tools/filter_photons.pyx'], include_dirs=['.']),
]

setup(name = 'photon-tools',
      author = 'Ben Gamari',
      author_email = 'bgamari@physics.umass.edu',
      url = 'http://goldnerlab.physics.umass.edu/',
      description = 'Tools for manipulating photon data from single-molecule experiments',
      version = '1.0',
      packages = ['photon_tools'],
      scripts = ['bin_photons', 'fcs-fit', 'plot-fret', 'alex.py', 'plot-bins'],
      license = 'GPLv3',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules,
)

