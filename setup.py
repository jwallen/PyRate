#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   PyRate - Python tools for computing chemical reaction rates
#
#   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
#                         Yury V. Suleimanov (ysuleyma@mit.edu)
#                         William H. Green (whgreen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a 
#   copy of this software and associated documentation files (the "Software"), 
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
#   and/or sell copies of the Software, and to permit persons to whom the 
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
#   DEALINGS IN THE SOFTWARE. 
#
################################################################################

try:
    from distutils.core import setup
    from distutils.extension import Extension
    from distutils.cmd import Command
except ImportError:
    print('The distutils package is required to install PyRate.')
    quit()

try:
    from Cython.Distutils import build_ext
    import Cython.Compiler.Options
    Cython.Compiler.Options.annotate = True
except ImportError:
    print('The cython package is required to install PyRate.')
    quit()

try:
    import numpy
except ImportError:
    print('The numpy package is required to install PyRate.')
    quit()

################################################################################

class test(Command):
    """
    Run the unit test suite. All files in the `tests` directory and its
    subdirectories that match the pattern ``*test.py`` will have their unit
    test suites executed.
    """

    description = "Run the unit test suite"
    user_options = []

    def initialize_options(self):
        self.verbosity = 2

    def finalize_options(self):
        try:
            self.verbosity = int(self.verbosity)
        except ValueError:
            raise ValueError('Verbosity must be an integer.')

    def run(self):
        import sys
        if sys.version.startswith('2.6') or sys.version.startswith('3.1'):
            import unittest2 as unittest
        else:
            import unittest
        suite = unittest.TestLoader().discover('tests', pattern='*test.py', top_level_dir='.')
        unittest.TextTestRunner(verbosity=self.verbosity).run(suite)

################################################################################

# The Cython extension modules to build
ext_modules = [
    Extension('pyrate.constants', ['pyrate/constants.py']),
    Extension('pyrate.kinetics.model', ['pyrate/kinetics/model.pyx']),
    Extension('pyrate.statmech.schrodinger', ['pyrate/statmech/schrodinger.pyx']),
    Extension('pyrate.statmech.translation', ['pyrate/statmech/translation.pyx']),
    Extension('pyrate.statmech.rotation', ['pyrate/statmech/rotation.pyx']),
    Extension('pyrate.statmech.vibration', ['pyrate/statmech/vibration.pyx']),
    Extension('pyrate.statmech.torsion', ['pyrate/statmech/torsion.pyx']),
    Extension('pyrate.thermo.converter', ['pyrate/thermo/converter.pyx']),
    Extension('pyrate.thermo.model', ['pyrate/thermo/model.pyx']),
    Extension('pyrate.species', ['pyrate/species.pyx']),
    Extension('pyrate.reaction', ['pyrate/reaction.pyx']),
]
for module in ext_modules:
    module.pyrex_directives = {'embedsignature': True}
    module.pyrex_c_in_temp = True

setup(
    name = 'PyRate',
    version = '0.1.0',
    description = 'Python tools for computing chemical reaction rates',
    author = 'Joshua W. Allen, Yury V. Suleimanov, William H. Green',
    author_email = 'pyrate_dev@mit.edu',
    url = 'http://github.com/GreenGroup/PyRate',
    packages = ['pyrate'],
    cmdclass = {'build_ext': build_ext, 'test': test},
    ext_modules = ext_modules,
    include_dirs = ['.', numpy.get_include()],
    requires = ['cython (>=0.15)', 'numpy (>=1.5.0)', 'scipy (>=0.9.0)', 'quantities (>=0.9.0)'],
    provides = ['pyrate'],
)

################################################################################

from numpy.distutils.core import setup, Extension

# The Fortran extension modules to build using f2py
ext_modules = [
    Extension('pyrate.rpmd._main', ['pyrate/rpmd/_main.pyf', 'pyrate/rpmd/_main.f90', 'pyrate/rpmd/_surface.f90'], libraries = ['fftw3']),
    Extension('pyrate.rpmd._surface', ['pyrate/rpmd/_surface.f90']),
]

setup(
    ext_modules = ext_modules,
)
