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
This script contains unit tests of the :mod:`pyrate.statmech.rotation` module.
"""

import unittest
import math
import numpy
import quantities as pq

from pyrate.statmech.rotation import *
import pyrate.constants as constants

################################################################################

class TestLinearRotor(unittest.TestCase):
    """
    Contains unit tests of the LinearRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 11.75 * pq.u * pq.angstrom * pq.angstrom
        self.symmetry = 2
        self.quantum = False
        self.mode = LinearRotor(
            inertia = self.inertia, 
            symmetry = self.symmetry, 
            quantum = self.quantum,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the LinearRotor.rotationalConstant property.
        """
        Bexp = 1.434692 * pq.wavenumber
        Bact = self.mode.rotationalConstant
        self.assertAlmostEqual(Bexp / Bact, 1.0, 4, '{0} != {1}'.format(Bexp, Bact))
        
    def test_setRotationalConstant(self):
        """
        Test setting the LinearRotor.rotationalConstant property.
        """
        self.mode.rotationalConstant *= 2
        self.assertAlmostEqual(self.mode.inertia / (self.inertia / 2), 1.0, 4, '{0} != {1}'.format(self.mode.inertia, self.inertia * 2))
        
    def test_getLevelEnergy(self):
        """
        Test the LinearRotor.getLevelEnergy() method.
        """
        B = float(self.mode.rotationalConstant) * constants.h * constants.c * 100.
        for J in range(0, 100):
            Eexp = B * J * (J + 1)
            Eact = self.mode.getLevelEnergy(J)
            if J == 0:
                self.assertEqual(Eact, 0, '{0} != {1}'.format(Eact, 0.0))
            else:
                self.assertAlmostEqual(Eexp / Eact, 1.0, 4, '{0} != {1}'.format(Eexp, Eact))
    
    def test_getLevelDegeneracy(self):
        """
        Test the LinearRotor.getLevelDegeneracy() method.
        """
        for J in range(0, 100):
            gexp = 2 * J + 1
            gact = self.mode.getLevelDegeneracy(J)
            self.assertEqual(gexp, gact, '{0} != {1}'.format(gact, gexp))
    
    def test_getPartitionFunction_classical(self):
        """
        Test the LinearRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([72.6691, 121.115, 242.230, 363.346, 484.461])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getPartitionFunction_quantum(self):
        """
        Test the LinearRotor.getPartitionFunction() method for a quantum
        rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([72.8360, 121.282, 242.391, 363.512, 484.627])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getHeatCapacity_classical(self):
        """
        Test the LinearRotor.getHeatCapacity() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1, 1, 1, 1, 1]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
    
    def test_getHeatCapacity_quantum(self):
        """
        Test the LinearRotor.getHeatCapacity() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1, 1, 1, 1, 1]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
       
    def test_getEnthalpy_classical(self):
        """
        Test the LinearRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = constants.R * Tlist / 1000.
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))
    
    def test_getEnthalpy_quantum(self):
        """
        Test the LinearRotor.getEnthalpy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.997705, 0.998624, 0.999312, 0.999541, 0.999656]) * constants.R * Tlist / 1000.
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))

    def test_getEntropy_classical(self):
        """
        Test the LinearRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([5.28592, 5.79674, 6.48989, 6.89535, 7.18304]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))
    
    def test_getEntropy_quantum(self):
        """
        Test the LinearRotor.getEntropy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([5.28592, 5.79674, 6.48989, 6.89535, 7.18304]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))

    def test_getSumOfStates_classical(self):
        """
        Test the LinearRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 2000*0.01196, 1.0*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        sumStates = self.mode.getSumOfStates(Elist)
        for n in range(1, len(Elist)):
            self.assertAlmostEqual(numpy.sum(densStates[0:n]) / sumStates[n], 1.0, 3)

    def test_getSumOfStates_quantum(self):
        """
        Test the LinearRotor.getSumOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 4000.*0.01196, 2.0*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        sumStates = self.mode.getSumOfStates(Elist)
        for n in range(1, len(Elist)):
            self.assertAlmostEqual(numpy.sum(densStates[0:n+1]) / sumStates[n], 1.0, 3)

    def test_getDensityOfStates_classical(self):
        """
        Test the LinearRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,400,500])
        Elist = numpy.arange(0, 4000.*0.01196, 1.0*0.01196)
        for T in Tlist:
            densStates = self.mode.getDensityOfStates(Elist)
            Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
            Qexp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 2, '{0} != {1} within 2 figures'.format(Qexp, Qact))

    def test_getDensityOfStates_quantum(self):
        """
        Test the LinearRotor.getDensityOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,400,500])
        Elist = numpy.arange(0, 4000.*0.01196, 2.0*0.01196)
        for T in Tlist:
            densStates = self.mode.getDensityOfStates(Elist)
            Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
            Qexp = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 2, '{0} != {1} within 2 figures'.format(Qexp, Qact))

    def test_repr(self):
        """
        Test that a LinearRotor object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia, mode.inertia, 6)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)
        
    def test_pickle(self):
        """
        Test that a LinearRotor object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertAlmostEqual(self.mode.inertia, mode.inertia, 6)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
