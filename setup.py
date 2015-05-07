#!/usr/bin/env python

from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy as np

from Cython.Distutils import build_ext

setup(name = 'photon-tools',
      author = 'Ben Gamari',
      author_email = 'bgamari@physics.umass.edu',
      url = 'http://goldnerlab.physics.umass.edu/',
      description = 'Tools for manipulating photon data from single-molecule experiments',
      version = '1.0',
      packages = ['photon_tools', 'photon_tools.io', 'photon_tools.correlate', 'photon_tools.utils'],
      scripts = ['bin_photons', 'fcs-fit', 'fcs-corr', 'plot-fret', 'plot-bins', 'lifetime-deconvolve',
                 'trim-stamps', 'anisotropy', 'fcs-mem', 'summarize-timestamps', 'imbalance', 'dls-mem'],
      license = 'GPLv3',
      install_requires=[
          'numpy',
          'scipy',
          'squmfit',
      ],
      cmdclass = { 'build_ext': build_ext },
      ext_modules = [
          Extension('photon_tools.bin_photons', ['photon_tools/bin_photons.pyx'],
                    include_dirs=['photon_tools']),
          Extension('photon_tools.filter_photons', ['photon_tools/filter_photons.pyx'],
                    include_dirs=['photon_tools']),
          Extension('photon_tools.io.timetag_parse', ['photon_tools/io/timetag_parse.pyx'],
                    include_dirs=[np.get_include(), 'photon_tools']),
      ],
)
