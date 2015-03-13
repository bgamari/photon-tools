#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

setup(name = 'photon-tools',
      author = 'Ben Gamari',
      author_email = 'bgamari@physics.umass.edu',
      url = 'http://goldnerlab.physics.umass.edu/',
      description = 'Tools for manipulating photon data from single-molecule experiments',
      version = '1.0',
      packages = ['photon_tools', 'photon_tools.io', 'photon_tools.correlate'],
      scripts = ['bin_photons', 'fcs-fit', 'fcs-corr', 'plot-fret', 'plot-bins', 'lifetime-deconvolve',
                 'trim-stamps', 'anisotropy', 'fcs-mem', 'summarize-timestamps'],
      license = 'GPLv3',
      install_requires=[
          'numpy',
          'scipy',
          'squmfit',
      ],
      ext_modules = cythonize([
          Extension('photon_tools',['photon_tools/bin_photons.pyx',
                                    'photon_tools/filter_photons.pyx',
                                    'photon_tools/io/timetag_parse.pyx'],
                                    include_dirs=[np.get_include()])
          ],
          include_path=['photon_tools'],
      )
)
