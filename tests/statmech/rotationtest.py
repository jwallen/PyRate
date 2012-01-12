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

class TestNonlinearRotor(unittest.TestCase):
    """
    Contains unit tests of the NonlinearRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = numpy.array([3.415, 16.65, 20.07]) * pq.amu * pq.angstrom * pq.angstrom
        self.symmetry = 4
        self.quantum = False
        self.mode = NonlinearRotor(
            inertia = self.inertia, 
            symmetry = self.symmetry, 
            quantum = self.quantum,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the NonlinearRotor.rotationalConstant property.
        """
        Bexp = numpy.array([4.93635, 1.0125, 0.839942]) * pq.wavenumber
        Bact = self.mode.rotationalConstant
        for B0, B in zip(Bexp, Bact):
            self.assertAlmostEqual(B0 / B, 1.0, 4, '{0} != {1}'.format(B0, B))
        
    def test_setRotationalConstant(self):
        """
        Test setting the NonlinearRotor.rotationalConstant property.
        """
        self.mode.rotationalConstant *= 2
        for I, I0 in zip(self.mode.inertia, self.inertia):
            self.assertAlmostEqual(I / (I0 / 2), 1.0, 4, '{0} != {1}'.format(I, I0 * 2))
        
    def test_getPartitionFunction_classical(self):
        """
        Test the NonlinearRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([651.162, 1401.08, 3962.84, 7280.21, 11208.6])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getHeatCapacity_classical(self):
        """
        Test the NonlinearRotor.getHeatCapacity() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.5, 1.5, 1.5, 1.5, 1.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
    
    def test_getEnthalpy_classical(self):
        """
        Test the NonlinearRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = 1.5 * constants.R * Tlist / 1000.
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))
    
    def test_getEntropy_classical(self):
        """
        Test the NonlinearRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([7.97876, 8.74500, 9.78472, 10.3929, 10.8244]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))
    
    def test_getSumOfStates_classical(self):
        """
        Test the NonlinearRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 1000*0.01196, 1*0.01196)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getDensityOfStates_classical(self):
        """
        Test the NonlinearRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 1000*0.01196, 1*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 3, '{0} != {1} within 3 places'.format(Qexp, Qact))

    def test_repr(self):
        """
        Test that a NonlinearRotor object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('mode = {0!r}'.format(self.mode))
        self.assertEqual(self.mode.inertia.shape, mode.inertia.shape)
        for I0, I in zip(self.mode.inertia, mode.inertia):
            self.assertAlmostEqual(I0, I, 6)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)
        
    def test_pickle(self):
        """
        Test that a NonlinearRotor object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertEqual(self.mode.inertia.shape, mode.inertia.shape)
        for I0, I in zip(self.mode.inertia, mode.inertia):
            self.assertAlmostEqual(I0, I, 6)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)

################################################################################

class TestKRotor(unittest.TestCase):
    """
    Contains unit tests of the KRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 11.75 * pq.u * pq.angstrom * pq.angstrom
        self.symmetry = 2
        self.quantum = False
        self.mode = KRotor(
            inertia = self.inertia,
            symmetry = self.symmetry, 
            quantum = self.quantum,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the KRotor.rotationalConstant property.
        """
        Bexp = 1.434692 * pq.wavenumber
        Bact = self.mode.rotationalConstant
        self.assertAlmostEqual(Bexp / Bact, 1.0, 4, '{0} != {1}'.format(Bexp, Bact))
        
    def test_setRotationalConstant(self):
        """
        Test setting the KRotor.rotationalConstant property.
        """
        self.mode.rotationalConstant *= 2
        self.assertAlmostEqual(self.mode.inertia / (self.inertia / 2), 1.0, 4, '{0} != {1}'.format(self.mode.inertia, self.inertia * 2))
        
    def test_getLevelEnergy(self):
        """
        Test the KRotor.getLevelEnergy() method.
        """
        B = float(self.mode.rotationalConstant) * constants.h * constants.c * 100.
        for J in range(0, 100):
            Eexp = float(B * J * J)
            Eact = float(self.mode.getLevelEnergy(J))
            if J == 0:
                self.assertEqual(Eact, 0, '{0} != {1}'.format(Eact, 0.0))
            else:
                self.assertAlmostEqual(Eexp / Eact, 1.0, 4, '{0} != {1}'.format(Eexp, Eact))
    
    def test_getLevelDegeneracy(self):
        """
        Test the KRotor.getLevelDegeneracy() method.
        """
        for J in range(0, 100):
            gexp = 1 if J == 0 else 2
            gact = self.mode.getLevelDegeneracy(J)
            self.assertEqual(gexp, gact, '{0} != {1}'.format(gact, gexp))
    
    def test_getPartitionFunction_classical(self):
        """
        Test the KRotor.getPartitionFunction() method for a classical
        rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([10.6839, 13.7929, 19.5060, 23.8899, 27.5857])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getPartitionFunction_quantum(self):
        """
        Test the KRotor.getPartitionFunction() method for a quantum
        rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([10.6839, 13.7929, 19.5060, 23.8899, 27.5857])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getHeatCapacity_classical(self):
        """
        Test the KRotor.getHeatCapacity() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
    
    def test_getHeatCapacity_quantum(self):
        """
        Test the KRotor.getHeatCapacity() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
       
    def test_getEnthalpy_classical(self):
        """
        Test the KRotor.getEnthalpy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))
    
    def test_getEnthalpy_quantum(self):
        """
        Test the KRotor.getEnthalpy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))

    def test_getEntropy_classical(self):
        """
        Test the KRotor.getEntropy() method using a classical rotor.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([2.86874, 3.12415, 3.47072, 3.67346, 3.81730]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))
    
    def test_getEntropy_quantum(self):
        """
        Test the KRotor.getEntropy() method using a quantum rotor.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([2.86874, 3.12415, 3.47072, 3.67346, 3.81730]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))

    def test_getSumOfStates_classical(self):
        """
        Test the KRotor.getSumOfStates() method using a classical rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 1000*0.01196, 1*0.01196)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.75 < numpy.sum(densStates[0:n+1]) / sumStates[n] < 1.3333, '{0} != {1}'.format(numpy.sum(densStates[0:n+1]), sumStates[n]))

    def test_getSumOfStates_quantum(self):
        """
        Test the KRotor.getSumOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 1000*0.01196, 1*0.01196)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n+1]) / sumStates[n] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n+1]), sumStates[n]))

    def test_getDensityOfStates_classical(self):
        """
        Test the KRotor.getDensityOfStates() method using a classical
        rotor.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 1000*0.01196, 1*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 500
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 0, '{0} != {1} within 0 places'.format(Qexp, Qact))

    def test_getDensityOfStates_quantum(self):
        """
        Test the KRotor.getDensityOfStates() method using a quantum rotor.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 1000*0.01196, 1*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 500
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 1, '{0} != {1} within 1 place'.format(Qexp, Qact))

    def test_repr(self):
        """
        Test that a KRotor object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia, mode.inertia, 6)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)
        
    def test_pickle(self):
        """
        Test that a KRotor object can be successfully pickled and unpickled
        with no loss of information.
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
