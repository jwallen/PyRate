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

"""
This script contains unit tests of the :mod:`pyrate.statmech.system` module.
"""

import unittest
import math
import numpy
import quantities as pq

from pyrate.statmech import *
import pyrate.constants as constants

################################################################################

class TestStatMech(unittest.TestCase):
    """
    Contains unit tests of the StatMech class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.ethylene = StatMech(
            modes = [
                IdealGasTranslation(mass=(28.03,"g/mol")),
                NonlinearRotor(inertia=([3.41526,16.6498,20.065],"amu*angstrom^2"), symmetry=4),
                HarmonicOscillator(frequencies=([828.397,970.652,977.223,1052.93,1233.55,1367.56,1465.09,1672.25,3098.46,3111.7,3165.79,3193.54], "cm^-1")),
            ],
            E0 = (0.0,"kJ/mol"),
            spinMultiplicity = 1,
        )
        self.oxygen = StatMech(
            modes = [
                IdealGasTranslation(mass=(31.99,"g/mol")),
                LinearRotor(inertia=(11.6056,"amu*angstrom^2"), symmetry=2),
                HarmonicOscillator(frequencies=([1621.54],"cm^-1")),
            ],
            E0 = (0.0,"kJ/mol"),
            spinMultiplicity = 3,
        )

    def test_getPartitionFunction_ethylene(self):
        """
        Test the StatMech.getPartitionFunction() method for ethylene.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([4.05311e+09, 4.19728e+10, 2.82309e+12, 7.51135e+13, 1.16538e+15])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.ethylene.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))

    def test_getHeatCapacity_ethylene(self):
        """
        Test the StatMech.getHeatCapacity() method for ethylene.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([5.11186, 7.40447, 11.1659, 13.1221, 14.1617]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.ethylene.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
    
    def test_getEnthalpy_ethylene(self):
        """
        Test the StatMech.getEnthalpy() method for ethylene.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([4.23129, 5.04826, 7.27337, 8.93167, 10.1223]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.ethylene.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))
    
    def test_getEntropy_ethylene(self):
        """
        Test the StatMech.getEntropy() method for ethylene.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([27.3540, 30.5085, 36.9422, 41.8817, 45.8142]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.ethylene.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))
    
    def test_getSumOfStates_ethylene(self):
        """
        Test the StatMech.getSumOfStates() method for ethylene.
        """
        Elist = numpy.arange(0, 5000*0.01196, 2*0.01196)
        sumStates = self.ethylene.getSumOfStates(Elist)
        densStates = self.ethylene.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n+1]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n+1]), sumStates[n]))
            
    def test_getDensityOfStates_ethylene(self):
        """
        Test the StatMech.getDensityOfStates() method for ethylene.
        """
        Elist = numpy.arange(0, 5000*0.01196, 2*0.01196)
        densStates = self.ethylene.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.ethylene.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 1, '{0} != {1} within 1 place'.format(Qexp, Qact))

    def test_getPartitionFunction_oxygen(self):
        """
        Test the StatMech.getPartitionFunction() method for oxygen.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([1.55584e+09, 9.38339e+09, 1.16459e+11, 5.51016e+11, 1.72794e+12])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.oxygen.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))

    def test_getHeatCapacity_oxygen(self):
        """
        Test the StatMech.getHeatCapacity() method for oxygen.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([3.52538, 3.70877, 4.14751, 4.32063, 4.39392]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.oxygen.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
    
    def test_getEnthalpy_oxygen(self):
        """
        Test the StatMech.getEnthalpy() method for oxygen.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([3.50326, 3.54432, 3.75062, 3.91623, 4.02765]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.oxygen.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))
    
    def test_getEntropy_oxygen(self):
        """
        Test the StatMech.getEntropy() method for oxygen.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([24.5699, 26.4079, 29.1328, 30.8526, 32.1070]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.oxygen.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))
    
    def test_getSumOfStates_oxygen(self):
        """
        Test the StatMech.getSumOfStates() method for oxygen.
        """
        Elist = numpy.arange(0, 5000*0.01196, 2*0.01196)
        sumStates = self.oxygen.getSumOfStates(Elist)
        densStates = self.oxygen.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n+1]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n+1]), sumStates[n]))
            
    def test_getDensityOfStates_oxygen(self):
        """
        Test the StatMech.getDensityOfStates() method for oxygen.
        """
        Elist = numpy.arange(0, 5000*0.01196, 2*0.01196)
        densStates = self.oxygen.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.oxygen.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 1, '{0} != {1} within 1 place'.format(Qexp, Qact))
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
