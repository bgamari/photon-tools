from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
        Extension('photon_tools.bin_photons', ['bin_photons.pyx'])
]

setup(name = 'photon-tools',
      author = 'Ben Gamari',
      author_email = 'bgamari@physics.umass.edu',
      url = 'http://goldnerlab.physics.umass.edu/',
      description = 'Tools for manipulating photon data from single-molecule experiments',
      version = '1.0',
      package_dir = {'photon_tools': '.'},
      packages = ['photon_tools'],
      scripts = ['bin_photons', 'fcs-fit'],
      license = 'GPLv3',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules,
)

