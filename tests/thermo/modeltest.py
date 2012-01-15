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
This script contains unit tests of the :mod:`pyrate.thermo.model` module.
"""

import unittest
import math
import numpy
import quantities as pq

from pyrate.thermo.model import *
import pyrate.constants as constants

################################################################################

class TestWilhoit(unittest.TestCase):
    """
    Contains unit tests of the :class:`Wilhoit` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.Cp0 = pq.Quantity(33.2579, "J/(mol*K)")
        self.CpInf = pq.Quantity(178.761, "J/(mol*K)")
        self.a0 = 0.0977518
        self.a1 = -16.3067
        self.a2 = 26.2524
        self.a3 = -12.6785
        self.B = pq.Quantity(1068.68, "K")
        self.H0 = pq.Quantity(-782.292, "kJ/mol")
        self.S0 = pq.Quantity(-984.932, "J/(mol*K)")
        self.Tmin = 300 * pq.K
        self.Tmax = 3000 * pq.K
        self.comment = 'C2H6'
        self.wilhoit = Wilhoit(
            Cp0 = self.Cp0,
            CpInf = self.CpInf,
            a0 = self.a0,
            a1 = self.a1,
            a2 = self.a2,
            a3 = self.a3,
            B = self.B,
            H0 = self.H0,
            S0 = self.S0,
            Tmin = self.Tmin,
            Tmax = self.Tmax,
            comment = self.comment,
        )
    
    def test_Cp0(self):
        """
        Test that the Wilhoit Cp0 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.Cp0 / self.Cp0, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.Cp0, self.Cp0))
    
    def test_CpInf(self):
        """
        Test that the Wilhoit CpInf property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.CpInf / self.CpInf, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.CpInf, self.CpInf))
    
    def test_a0(self):
        """
        Test that the Wilhoit a0 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.a0 / self.a0, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.a0, self.a0))
    
    def test_a1(self):
        """
        Test that the Wilhoit a1 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.a1 / self.a1, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.a1, self.a1))
    
    def test_a2(self):
        """
        Test that the Wilhoit a2 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.a2 / self.a2, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.a2, self.a2))
    
    def test_a3(self):
        """
        Test that the Wilhoit a3 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.a3 / self.a3, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.a3, self.a3))
    
    def test_B(self):
        """
        Test that the Wilhoit B property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.B / self.B, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.B, self.B))
    
    def test_H0(self):
        """
        Test that the Wilhoit H0 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.H0 / self.H0, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.H0, self.H0))
    
    def test_S0(self):
        """
        Test that the Wilhoit S0 property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.S0 / self.S0, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.S0, self.S0))
    
    def test_Tmin(self):
        """
        Test that the Wilhoit Tmin property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.Tmin / self.Tmin, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.Tmin, self.Tmin))
    
    def test_Tmax(self):
        """
        Test that the Wilhoit Tmax property was properly set.
        """
        self.assertAlmostEqual(self.wilhoit.Tmax / self.Tmax, 1.0, 6, '{0} != {1} within 6 places'.format(self.wilhoit.Tmax, self.Tmax))
    
    def test_Comment(self):
        """
        Test that the Wilhoit comment property was properly set.
        """
        self.assertEqual(self.wilhoit.comment, self.comment)
    
    def test_isTemperatureValid(self):
        """
        Test the Wilhoit.isTemperatureValid() method.
        """
        Tdata = [200,400,600,800,1000,1200,1400,1600,1800,2000] * pq.K
        validdata = numpy.array([False,True,True,True,True,True,True,True,True,True], numpy.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.wilhoit.isTemperatureValid(T)
            self.assertEqual(valid0, valid)
        
    def test_getHeatCapacity(self):
        """
        Test the Wilhoit.getHeatCapacity() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        Cpexplist = numpy.array([5.12003, 7.80327, 10.5528, 12.8323, 14.6013, 15.9511, 16.9842, 17.7837, 18.4114, 18.9117]) * constants.R
        for T, Cpexp in zip(Tlist, Cpexplist):
            Cpact = self.wilhoit.getHeatCapacity(T)
            self.assertAlmostEqual(Cpexp / Cpact, 1.0, 3, '{0} != {1} within 3 places'.format(Cpexp, Cpact))
       
    def test_getEnthalpy(self):
        """
        Test the Wilhoit.getEnthalpy() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        Hexplist = numpy.array([-51.9303, -22.7609, -12.1050, -6.14444, -2.16433, 0.747500, 2.99646, 4.79698, 6.27618, 7.51564]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.wilhoit.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 3, '{0} != {1}'.format(Hexp, Hact))
                
    def test_getEntropy(self):
        """
        Test the Wilhoit.getEntropy() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        Sexplist = numpy.array([25.3095, 29.6445, 33.3398, 36.7006, 39.7629, 42.5499, 45.0898, 47.4122, 49.5445, 51.5112]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.wilhoit.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))

    def test_getFreeEnergy(self):
        """
        Test the Wilhoit.getFreeEnergy() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        for T in Tlist:
            Gexp = self.wilhoit.getEnthalpy(T) - 0.001 * T * self.wilhoit.getEntropy(T)
            Gact = self.wilhoit.getFreeEnergy(T)
            self.assertAlmostEqual(Gexp / Gact, 1.0, 4, '{0} != {1}'.format(Gexp, Gact))
    
    def test_pickle(self):
        """
        Test that a Wilhoit object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        wilhoit = cPickle.loads(cPickle.dumps(self.wilhoit))
        self.assertEqual(self.wilhoit.Cp0, wilhoit.Cp0)
        self.assertEqual(self.wilhoit.CpInf, wilhoit.CpInf)
        self.assertEqual(self.wilhoit.a0, wilhoit.a0)
        self.assertEqual(self.wilhoit.a1, wilhoit.a1)
        self.assertEqual(self.wilhoit.a2, wilhoit.a2)
        self.assertEqual(self.wilhoit.a3, wilhoit.a3)
        self.assertEqual(self.wilhoit.B, wilhoit.B)
        self.assertEqual(self.wilhoit.H0, wilhoit.H0)
        self.assertEqual(self.wilhoit.S0, wilhoit.S0)
        self.assertEqual(self.wilhoit.Tmin, wilhoit.Tmin)
        self.assertEqual(self.wilhoit.Tmax, wilhoit.Tmax)
        self.assertEqual(self.wilhoit.comment, wilhoit.comment)
    
    def test_repr(self):
        """
        Test that a Wilhoit object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        exec('wilhoit = {0!r}'.format(self.wilhoit))
        self.assertEqual(self.wilhoit.Cp0, wilhoit.Cp0)
        self.assertEqual(self.wilhoit.CpInf, wilhoit.CpInf)
        self.assertEqual(self.wilhoit.a0, wilhoit.a0)
        self.assertEqual(self.wilhoit.a1, wilhoit.a1)
        self.assertEqual(self.wilhoit.a2, wilhoit.a2)
        self.assertEqual(self.wilhoit.a3, wilhoit.a3)
        self.assertEqual(self.wilhoit.B, wilhoit.B)
        self.assertEqual(self.wilhoit.H0, wilhoit.H0)
        self.assertEqual(self.wilhoit.S0, wilhoit.S0)
        self.assertEqual(self.wilhoit.Tmin, wilhoit.Tmin)
        self.assertEqual(self.wilhoit.Tmax, wilhoit.Tmax)
        self.assertEqual(self.wilhoit.comment, wilhoit.comment)

################################################################################

class TestNASA(unittest.TestCase):
    """
    Contains unit tests of the NASA class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.coeffs = [-0.307954,0.0245269,-1.2413e-05,3.07724e-09,-3.01467e-13,-10693,22.628]
        self.Tmin = 650.73 * pq.K
        self.Tmax = 3000.0 * pq.K
        self.comment = "C2H6"
        self.nasa = NASA(
            coeffs = self.coeffs,
            Tmin = self.Tmin, 
            Tmax = self.Tmax, 
            comment = self.comment,
        )
    
    def test_c0(self):
        """
        Test that the NASA c0 property was properly set.
        """
        self.assertAlmostEqual(self.nasa.coeffs[0] / self.coeffs[0], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.coeffs[0], self.coeffs[0]))
        self.assertAlmostEqual(self.nasa.c0 / self.coeffs[0], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.c0, self.coeffs[0]))
    
    def test_c1(self):
        """
        Test that the NASA c1 property was properly set.
        """
        self.assertAlmostEqual(self.nasa.coeffs[1] / self.coeffs[1], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.coeffs[1], self.coeffs[1]))
        self.assertAlmostEqual(self.nasa.c1 / self.coeffs[1], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.c1, self.coeffs[1]))
    
    def test_c2(self):
        """
        Test that the NASA c2 property was properly set.
        """
        self.assertAlmostEqual(self.nasa.coeffs[2] / self.coeffs[2], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.coeffs[2], self.coeffs[2]))
        self.assertAlmostEqual(self.nasa.c2 / self.coeffs[2], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.c2, self.coeffs[2]))
    
    def test_c3(self):
        """
        Test that the NASA c3 property was properly set.
        """
        self.assertAlmostEqual(self.nasa.coeffs[3] / self.coeffs[3], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.coeffs[3], self.coeffs[3]))
        self.assertAlmostEqual(self.nasa.c3 / self.coeffs[3], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.c3, self.coeffs[3]))
    
    def test_c4(self):
        """
        Test that the NASA c4 property was properly set.
        """
        self.assertAlmostEqual(self.nasa.coeffs[4] / self.coeffs[4], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.coeffs[4], self.coeffs[4]))
        self.assertAlmostEqual(self.nasa.c4 / self.coeffs[4], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.c4, self.coeffs[4]))
    
    def test_c5(self):
        """
        Test that the NASA c5 property was properly set.
        """
        self.assertAlmostEqual(self.nasa.coeffs[5] / self.coeffs[5], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.coeffs[5], self.coeffs[5]))
        self.assertAlmostEqual(self.nasa.c5 / self.coeffs[5], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.c5, self.coeffs[5]))
    
    def test_c6(self):
        """
        Test that the NASA c6 property was properly set.
        """
        self.assertAlmostEqual(self.nasa.coeffs[6] / self.coeffs[6], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.coeffs[6], self.coeffs[6]))
        self.assertAlmostEqual(self.nasa.c6 / self.coeffs[6], 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.c6, self.coeffs[6]))
    
    def test_cm1(self):
        """
        Test that the NASA cm1 property was properly set.
        """
        self.assertEqual(self.nasa.cm1, 0, '{0} != {1}'.format(self.nasa.cm1, 0))
    
    def test_cm2(self):
        """
        Test that the NASA cm2 property was properly set.
        """
        self.assertEqual(self.nasa.cm2, 0, '{0} != {1}'.format(self.nasa.cm2, 0))
    
    def test_Tmin(self):
        """
        Test that the NASA Tmin property was properly set.
        """
        self.assertAlmostEqual(self.nasa.Tmin / self.Tmin, 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.Tmin, self.Tmin))
    
    def test_Tmax(self):
        """
        Test that the NASA Tmax property was properly set.
        """
        self.assertAlmostEqual(self.nasa.Tmax / self.Tmax, 1.0, 6, '{0} != {1} within 6 places'.format(self.nasa.Tmax, self.Tmax))
    
    def test_comment(self):
        """
        Test that the NASA comment property was properly set.
        """
        self.assertEqual(self.nasa.comment, self.comment)
    
    def test_isTemperatureValid(self):
        """
        Test the NASA.isTemperatureValid() method.
        """
        Tdata = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        validdata = numpy.array([False,False,False,True,True,True,True,True,True,True], numpy.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.nasa.isTemperatureValid(T)
            self.assertEqual(valid0, valid)
        
    def test_getHeatCapacity(self):
        """
        Test the NASA.getHeatCapacity() method.
        """
        Tlist = numpy.array([800,1000,1200,1400,1600,1800,2000])
        Cpexplist = numpy.array([12.8213, 14.5817, 15.9420, 16.9861, 17.78645, 18.4041, 18.8883]) * constants.R
        for T, Cpexp in zip(Tlist, Cpexplist):
            Cpact = self.nasa.getHeatCapacity(T)
            self.assertAlmostEqual(Cpexp / Cpact, 1.0, 4, '{0} != {1}'.format(Cpexp, Cpact))
        
    def test_getEnthalpy(self):
        """
        Test the NASA.getEnthalpy() method.
        """
        Tlist = numpy.array([800,1000,1200,1400,1600,1800,2000])
        Hexplist = numpy.array([-6.14236, -2.16615, 0.743456, 2.99256, 4.79397, 6.27334, 7.51156]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.nasa.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 3, '{0} != {1}'.format(Hexp, Hact))
        
    def test_getEntropy(self):
        """
        Test the NASA.getEntropy() method.
        """
        Tlist = numpy.array([800,1000,1200,1400,1600,1800,2000])
        Sexplist = numpy.array([36.7131, 39.7715, 42.5557, 45.0952, 47.4179, 49.5501, 51.5152]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.nasa.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))

    def test_getFreeEnergy(self):
        """
        Test the NASA.getFreeEnergy() method.
        """
        Tlist = numpy.array([800,1000,1200,1400,1600,1800,2000])
        for T in Tlist:
            Gexp = self.nasa.getEnthalpy(T) - 0.001 * T * self.nasa.getEntropy(T)
            Gact = self.nasa.getFreeEnergy(T)
            self.assertAlmostEqual(Gexp / Gact, 1.0, 4, '{0} != {1}'.format(Gexp, Gact))
    
    def test_pickle(self):
        """
        Test that a NASA object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        nasa = cPickle.loads(cPickle.dumps(self.nasa))
        self.assertAlmostEqual(self.nasa.cm2, nasa.cm2, 4)
        self.assertAlmostEqual(self.nasa.cm1, nasa.cm1, 4)
        self.assertAlmostEqual(self.nasa.c0 / nasa.c0, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c1 / nasa.c1, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c2 / nasa.c2, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c3 / nasa.c3, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c4 / nasa.c4, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c5 / nasa.c5, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c6 / nasa.c6, 1.0, 4)
        self.assertEqual(self.nasa.Tmin, nasa.Tmin)
        self.assertEqual(self.nasa.Tmax, nasa.Tmax)
        self.assertEqual(self.nasa.comment, nasa.comment)

    def test_repr(self):
        """
        Test that a NASA object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        exec('nasa = {0!r}'.format(self.nasa))
        self.assertAlmostEqual(self.nasa.cm2, nasa.cm2, 4)
        self.assertAlmostEqual(self.nasa.cm1, nasa.cm1, 4)
        self.assertAlmostEqual(self.nasa.c0 / nasa.c0, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c1 / nasa.c1, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c2 / nasa.c2, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c3 / nasa.c3, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c4 / nasa.c4, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c5 / nasa.c5, 1.0, 4)
        self.assertAlmostEqual(self.nasa.c6 / nasa.c6, 1.0, 4)
        self.assertEqual(self.nasa.Tmin, nasa.Tmin)
        self.assertEqual(self.nasa.Tmax, nasa.Tmax)
        self.assertEqual(self.nasa.comment, nasa.comment)

################################################################################

class TestMultiNASA(unittest.TestCase):
    """
    Contains unit tests of the MultiNASA class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.nasa0 = NASA(Tmin=300 * pq.K, Tmax=650.73 * pq.K, coeffs=[4.03055,-0.00214171,4.90611e-05,-5.99027e-08,2.38945e-11,-11257.6,3.5613], comment="""Low temperature range polynomial""")
        self.nasa1 = NASA(Tmin=650.73 * pq.K, Tmax=3000 * pq.K, coeffs=[-0.307954,0.0245269,-1.2413e-05,3.07724e-09,-3.01467e-13,-10693,22.628], comment="""High temperature range polynomial""")
        self.multiNASA = MultiNASA(
            polynomials=[self.nasa0, self.nasa1],
            Tmin = 300.0*pq.K, 
            Tmax = 3000.0*pq.K, 
            comment = "C2H6",
        )
    
    def test_isTemperatureValid(self):
        """
        Test the MultiNASA.isTemperatureValid() method.
        """
        Tdata = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        validdata = numpy.array([False,True,True,True,True,True,True,True,True,True], numpy.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.multiNASA.isTemperatureValid(T)
            self.assertEqual(valid0, valid)
        
    def test_getHeatCapacity(self):
        """
        Test the MultiNASA.getHeatCapacity() method.
        """
        Tlist = numpy.array([400,600,800,1000,1200,1400,1600,1800,2000])
        Cpexplist = numpy.array([7.80157, 10.5653, 12.8213, 14.5817, 15.9420, 16.9861, 17.78645, 18.4041, 18.8883]) * constants.R
        for T, Cpexp in zip(Tlist, Cpexplist):
            Cpact = self.multiNASA.getHeatCapacity(T)
            self.assertAlmostEqual(Cpexp / Cpact, 1.0, 4, '{0} != {1}'.format(Cpexp, Cpact))
        
    def test_getEnthalpy(self):
        """
        Test the MultiNASA.getEnthalpy() method.
        """
        Tlist = numpy.array([400,600,800,1000,1200,1400,1600,1800,2000])
        Hexplist = numpy.array([-22.7613, -12.1027, -6.14236, -2.16615, 0.743456, 2.99256, 4.79397, 6.27334, 7.51156]) * 0.001 * constants.R * Tlist
        for T, Hexp in zip(Tlist, Hexplist):
            Hact = self.multiNASA.getEnthalpy(T)
            self.assertAlmostEqual(Hexp / Hact, 1.0, 3, '{0} != {1}'.format(Hexp, Hact))
        
    def test_getEntropy(self):
        """
        Test the MultiNASA.getEntropy() method.
        """
        Tlist = numpy.array([400,600,800,1000,1200,1400,1600,1800,2000])
        Sexplist = numpy.array([29.6534, 33.3516, 36.7131, 39.7715, 42.5557, 45.0952, 47.4179, 49.5501, 51.5152]) * constants.R
        for T, Sexp in zip(Tlist, Sexplist):
            Sact = self.multiNASA.getEntropy(T)
            self.assertAlmostEqual(Sexp / Sact, 1.0, 4, '{0} != {1}'.format(Sexp, Sact))

    def test_getFreeEnergy(self):
        """
        Test the MultiNASA.getFreeEnergy() method.
        """
        Tlist = numpy.array([400,600,800,1000,1200,1400,1600,1800,2000])
        for T in Tlist:
            Gexp = self.multiNASA.getEnthalpy(T) - 0.001 * T * self.multiNASA.getEntropy(T)
            Gact = self.multiNASA.getFreeEnergy(T)
            self.assertAlmostEqual(Gexp / Gact, 1.0, 4, '{0} != {1}'.format(Gexp, Gact))
    
    def test_pickle(self):
        """
        Test that a MultiNASA object can be successfully pickled and unpickled
        with no loss of information.
        """
        import cPickle
        multiNASA = cPickle.loads(cPickle.dumps(self.multiNASA))
        self.assertEqual(len(self.multiNASA.polynomials), len(multiNASA.polynomials))
        for poly0, poly in zip(self.multiNASA.polynomials, multiNASA.polynomials):
            self.assertEqual(poly0.cm2, poly.cm2)
            self.assertEqual(poly0.cm1, poly.cm1)
            self.assertEqual(poly0.c0, poly.c0)
            self.assertEqual(poly0.c1, poly.c1)
            self.assertEqual(poly0.c2, poly.c2)
            self.assertEqual(poly0.c3, poly.c3)
            self.assertEqual(poly0.c4, poly.c4)
            self.assertEqual(poly0.c5, poly.c5)
            self.assertEqual(poly0.c6, poly.c6)
            self.assertEqual(poly0.Tmin, poly.Tmin)
            self.assertEqual(poly0.Tmax, poly.Tmax)
            self.assertEqual(poly0.comment, poly.comment)
        self.assertEqual(self.multiNASA.Tmin, multiNASA.Tmin)
        self.assertEqual(self.multiNASA.Tmax, multiNASA.Tmax)
        self.assertEqual(self.multiNASA.comment, multiNASA.comment)

    def test_repr(self):
        """
        Test that a MultiNASA object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        exec('multiNASA = {0!r}'.format(self.multiNASA))
        self.assertEqual(len(self.multiNASA.polynomials), len(multiNASA.polynomials))
        for poly0, poly in zip(self.multiNASA.polynomials, multiNASA.polynomials):
            self.assertEqual(poly0.cm2, poly.cm2)
            self.assertEqual(poly0.cm1, poly.cm1)
            self.assertEqual(poly0.c0, poly.c0)
            self.assertEqual(poly0.c1, poly.c1)
            self.assertEqual(poly0.c2, poly.c2)
            self.assertEqual(poly0.c3, poly.c3)
            self.assertEqual(poly0.c4, poly.c4)
            self.assertEqual(poly0.c5, poly.c5)
            self.assertEqual(poly0.c6, poly.c6)
            self.assertEqual(poly0.Tmin, poly.Tmin)
            self.assertEqual(poly0.Tmax, poly.Tmax)
            self.assertEqual(poly0.comment, poly.comment)
        self.assertEqual(self.multiNASA.Tmin, multiNASA.Tmin)
        self.assertEqual(self.multiNASA.Tmax, multiNASA.Tmax)
        self.assertEqual(self.multiNASA.comment, multiNASA.comment)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
