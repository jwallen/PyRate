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
This script contains unit tests of the :mod:`pyrate.statmech.translation` module.
"""

import unittest
import math
import numpy
import quantities as pq

from pyrate.statmech.translation import *
import pyrate.constants as constants

################################################################################

class TestIdealGasTranslation(unittest.TestCase):
    """
    Contains unit tests of the IdealGasTranslation class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.mass = 32.0 * pq.u
        self.quantum = False
        self.active = False
        self.mode = IdealGasTranslation(
            mass = self.mass, 
            quantum = self.quantum, 
            active = self.active, 
        )
        
    def test_getPartitionFunction_classical(self):
        """
        Test the IdealGasTranslation.getPartitionFunction() method for a
        classical translator.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([7.22597e+06, 2.59130e+07, 1.46586e+08, 4.03944e+08, 8.29217e+08])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getHeatCapacity_classical(self):
        """
        Test the IdealGasTranslation.getHeatCapacity() method using a classical
        translator.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
    
    def test_getEnthalpy_classical(self):
        """
        Test the IdealGasTranslation.getEnthalpy() method using a classical
        translator.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1.5, 1.5, 1.5, 1.5, 1.5]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))
    
    def test_getEntropy_classical(self):
        """
        Test the IdealGasTranslation.getEntropy() method using a classical
        translator.
        """
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([18.2932, 19.5703, 21.3031, 22.3168, 23.0360]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))
    
    def test_getSumOfStates_classical(self):
        """
        Test the IdealGasTranslation.getSumOfStates() method using a classical
        translator.
        """
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getDensityOfStates_classical(self):
        """
        Test the IdealGasTranslation.getDensityOfStates() method using a 
        classical translator.
        """
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 1, '{0} != {1} within 1 place'.format(Qexp, Qact))

    def test_repr(self):
        """
        Test that an IdealGasTranslation object can be successfully
        reconstructed from its repr() output with no loss of information.
        """
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.mass, mode.mass, 6)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)
        
    def test_pickle(self):
        """
        Test that a IdealGasTranslation object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertAlmostEqual(self.mode.mass, mode.mass, 6)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
