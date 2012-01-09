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
This script contains unit tests of the :mod:`pyrate.statmech.torsion` module.
"""

import unittest
import math
import numpy
import quantities as pq

from pyrate.statmech.torsion import *
import pyrate.constants as constants

################################################################################

class TestHinderedRotor(unittest.TestCase):
    """
    Contains unit tests of the HinderedRotor class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.inertia = 1.56764 * pq.amu * pq.angstrom * pq.angstrom
        self.symmetry = 3
        self.barrier = pq.Quantity(11.373,"kJ/mol") 
        self.fourier = numpy.array([ [4.58375, 0.841648, -5702.71, 6.02657, 4.7446], [0.726951, -0.677255, 0.207032, 0.553307, -0.503303] ]) * pq.J / pq.mol
        self.quantum = True
        self.active = True
        self.mode = HinderedRotor(
            inertia = self.inertia, 
            symmetry = self.symmetry, 
            barrier = self.barrier,
            fourier = self.fourier,
            quantum = self.quantum,
            active = self.active,
        )
        
    def test_getRotationalConstant(self):
        """
        Test getting the HinderedRotor.rotationalConstant property.
        """
        Bexp = 10.7535 * pq.wavenumber
        Bact = self.mode.rotationalConstant
        self.assertAlmostEqual(Bexp / Bact, 1.0, 4, '{0} != {1}'.format(Bexp, Bact))
        
    def test_setRotationalConstant(self):
        """
        Test setting the HinderedRotor.rotationalConstant property.
        """
        self.mode.rotationalConstant *= 2
        self.assertAlmostEqual(self.mode.inertia / (self.inertia / 2), 1.0, 4, '{0} != {1}'.format(self.mode.inertia, self.inertia * 2))
    
    def test_getPotential_cosine(self):
        """
        Test the HinderedRotor.getPotential() method for a cosine potential.
        """
        self.mode.fourier = None
        phi = numpy.arange(0.0, 2 * constants.pi + 0.0001, constants.pi / 24.)
        V = numpy.zeros_like(phi)
        for i in range(phi.shape[0]):
            V[i] = self.mode.getPotential(phi[i])
    
    def test_getPotential_fourier(self):
        """
        Test the HinderedRotor.getPotential() method for a Fourier series
        potential.
        """
        phi = numpy.arange(0.0, 2 * constants.pi + 0.0001, constants.pi / 24.)
        V = numpy.zeros_like(phi)
        for i in range(phi.shape[0]):
            V[i] = self.mode.getPotential(phi[i])

    def test_getPartitionFunction_classical_cosine(self):
        """
        Test the HinderedRotor.getPartitionFunction() method for a cosine
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([0.741953, 1.30465, 2.68553, 3.88146, 4.91235])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getPartitionFunction_classical_fourier(self):
        """
        Test the HinderedRotor.getPartitionFunction() method for a Fourier
        series potential in the classical limit.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([0.745526, 1.30751, 2.68722, 3.88258, 4.91315])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getPartitionFunction_quantum_cosine(self):
        """
        Test the HinderedRotor.getPartitionFunction() method for a cosine
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([2.60153, 3.35860, 4.74983, 5.81734, 6.71730])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getPartitionFunction_quantum_fourier(self):
        """
        Test the HinderedRotor.getPartitionFunction() method for a Fourier
        series potential in the quantum limit.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Qexplist = numpy.array([2.60153, 3.35860, 4.74983, 5.81734, 6.71730])
        for T, Qexp in zip(Tlist, Qexplist):
            Qact = self.mode.getPartitionFunction(T)
            self.assertAlmostEqual(Qexp / Qact, 1.0, 4, '{0} != {1}'.format(Qexp, Qact))
            
    def test_getHeatCapacity_classical_cosine(self):
        """
        Test the HinderedRotor.getHeatCapacity() method using a cosine
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.01741, 0.951141, 0.681919, 0.589263, 0.552028]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
    
    def test_getHeatCapacity_classical_fourier(self):
        """
        Test the HinderedRotor.getHeatCapacity() method using a Fourier series
        potential in the classical limit.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([1.17682, 1.01369, 0.698588, 0.596797, 0.556293]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
       
    def test_getHeatCapacity_quantum_cosine(self):
        """
        Test the HinderedRotor.getHeatCapacity() method using a cosine
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
    
    def test_getHeatCapacity_quantum_fourier(self):
        """
        Test the HinderedRotor.getHeatCapacity() method using a Fourier series
        potential in the quantum limit.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Cvexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * constants.R
        for T, Cvexp in zip(Tlist, Cvexplist):
            Cvact = self.mode.getHeatCapacity(T)
            self.assertAlmostEqual(Cvexp / Cvact, 1.0, 4, '{0} != {1}'.format(Cvexp, Cvact))
       
    def test_getEnthalpy_classical_cosine(self):
        """
        Test the HinderedRotor.getEnthalpy() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1.09556, 1.09949, 0.962738, 0.854617, 0.784333]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))
    
    def test_getEnthalpy_classical_fourier(self):
        """
        Test the HinderedRotor.getEnthalpy() method using a Fourier series 
        potential in the classical limit.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([1.08882, 1.09584, 0.961543, 0.854054, 0.784009]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 3, '{0} != {1}'.format(Hexp, Hact))

    def test_getEnthalpy_quantum_cosine(self):
        """
        Test the HinderedRotor.getEnthalpy() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 4, '{0} != {1}'.format(Hexp, Hact))
    
    def test_getEnthalpy_quantum_fourier(self):
        """
        Test the HinderedRotor.getEnthalpy() method using a Fourier series 
        potential in the quantum limit.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Hexplist = numpy.array([0.5, 0.5, 0.5, 0.5, 0.5]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.mode.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 3, '{0} != {1}'.format(Hexp, Hact))

    def test_getEntropy_classical_cosine(self):
        """
        Test the HinderedRotor.getEntropy() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.798184, 1.36653, 1.95158, 2.21168, 2.37687]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))
    
    def test_getEntropy_classical_fourier(self):
        """
        Test the HinderedRotor.getEntropy() method using a Fourier series 
        potential in the classical limit.
        """
        self.mode.quantum = False
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.796242, 1.36506, 1.95101, 2.21141, 2.37671]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))

    def test_getEntropy_quantum_cosine(self):
        """
        Test the HinderedRotor.getEntropy() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.956098, 1.21152, 1.55811, 1.76084, 1.90469]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))
    
    def test_getEntropy_quantum_fourier(self):
        """
        Test the HinderedRotor.getEntropy() method using a Fourier series 
        potential in the quantum limit.
        """
        self.mode.quantum = True
        Tlist = numpy.array([300,500,1000,1500,2000])
        Sexplist = numpy.array([0.956098, 1.21152, 1.55811, 1.76084, 1.90469]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.mode.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))

    def test_getSumOfStates_classical_cosine(self):
        """
        Test the HinderedRotor.getSumOfStates() method using a cosine potential
        in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getSumOfStates_classical_fourier(self):
        """
        Test the HinderedRotor.getSumOfStates() method using a Fourier series
        potential in the classical limit.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        try:
            sumStates = self.mode.getSumOfStates(Elist)
            self.fail('NotImplementedError not raised by HinderedRotor.getSumOfStates()')
        except NotImplementedError:
            pass
        
    def test_getSumOfStates_quantum_cosine(self):
        """
        Test the HinderedRotor.getSumOfStates() method using a cosine potential
        in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getSumOfStates_quantum_fourier(self):
        """
        Test the HinderedRotor.getSumOfStates() method using a Fourier series
        potential in the quantum limit.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        sumStates = self.mode.getSumOfStates(Elist)
        densStates = self.mode.getDensityOfStates(Elist)
        for n in range(10, len(Elist)):
            self.assertTrue(0.8 < numpy.sum(densStates[0:n]) / sumStates[n-1] < 1.25, '{0} != {1}'.format(numpy.sum(densStates[0:n]), sumStates[n]))

    def test_getDensityOfStates_classical_cosine(self):
        """
        Test the HinderedRotor.getDensityOfStates() method using a classical
        potential in the classical limit.
        """
        self.mode.quantum = False
        self.mode.fourier = None
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 1, '{0} != {1} within 1 place'.format(Qexp, Qact))

    def test_getDensityOfStates_classical_fourier(self):
        """
        Test the HinderedRotor.getDensityOfStates() method using a Fourier 
        series potential in the classical limit.
        """
        self.mode.quantum = False
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        try:
            densStates = self.mode.getDensityOfStates(Elist)
            self.fail('NotImplementedError not raised by HinderedRotor.getDensityOfStates()')
        except NotImplementedError:
            pass
        
    def test_getDensityOfStates_quantum_cosine(self):
        """
        Test the HinderedRotor.getDensityOfStates() method using a classical
        potential in the quantum limit.
        """
        self.mode.quantum = True
        self.mode.fourier = None
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 1, '{0} != {1} within 1 place'.format(Qexp, Qact))

    def test_getDensityOfStates_quantum_fourier(self):
        """
        Test the HinderedRotor.getDensityOfStates() method using a Fourier 
        series potential in the quantum limit.
        """
        self.mode.quantum = True
        Elist = numpy.arange(0, 10000*0.01196, 1*0.01196)
        densStates = self.mode.getDensityOfStates(Elist)
        T = 100
        Qact = numpy.sum(densStates * numpy.exp(-Elist * 1000. / constants.R / T))
        Qexp = self.mode.getPartitionFunction(T)
        self.assertAlmostEqual(Qexp / Qact, 1.0, 1, '{0} != {1} within 1 place'.format(Qexp, Qact))

    def test_repr(self):
        """
        Test that a HinderedRotor object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('mode = {0!r}'.format(self.mode))
        self.assertAlmostEqual(self.mode.inertia, mode.inertia, 6)
        self.assertAlmostEqual(self.mode.barrier, mode.barrier, 6)
        self.assertEqual(self.mode.fourier.shape, mode.fourier.shape)
        for A0, A in zip(self.mode.fourier.flat, mode.fourier.flat):
            self.assertAlmostEqual(A0 / A, 1.0, 6)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)
        
    def test_pickle(self):
        """
        Test that a HinderedRotor object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        mode = cPickle.loads(cPickle.dumps(self.mode))
        self.assertAlmostEqual(self.mode.inertia, mode.inertia, 6)
        self.assertAlmostEqual(self.mode.barrier, mode.barrier, 6)
        self.assertEqual(self.mode.fourier.shape, mode.fourier.shape)
        for A0, A in zip(self.mode.fourier.flat, mode.fourier.flat):
            self.assertAlmostEqual(A0 / A, 1.0, 6)
        self.assertEqual(self.mode.symmetry, mode.symmetry)
        self.assertEqual(self.mode.quantum, mode.quantum)
        self.assertEqual(self.mode.active, mode.active)
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
