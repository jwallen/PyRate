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

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
