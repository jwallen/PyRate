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
This script contains unit tests of the :mod:`pyrate.kinetics.model` module.
"""

import unittest
import math
import numpy
import quantities as pq

from pyrate.kinetics.model import *

################################################################################

class TestArrhenius(unittest.TestCase):
    """
    Contains unit tests of the :class:`Arrhenius` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.A = pq.Quantity(1.0e12,"cm^3/(mol*s)")
        self.n = 0.5
        self.Ea = pq.Quantity(10.0, "kcal/mol")
        self.T0 = 1 * pq.K
        self.Tmin = 300 * pq.K
        self.Tmax = 3000 * pq.K
        self.order = 2
        self.comment = 'C2H6'
        self.arrhenius = Arrhenius(
            A = self.A,
            n = self.n,
            Ea = self.Ea,
            T0 = self.T0,
            Tmin = self.Tmin,
            Tmax = self.Tmax,
            order = self.order,
            comment = self.comment,
        )
    
    def test_A(self):
        """
        Test that the Arrhenius A property was properly set.
        """
        self.assertAlmostEqual(self.arrhenius.A / self.A, 1.0, 6, '{0} != {1} within 6 places'.format(self.arrhenius.A, self.A))
        
    def test_n(self):
        """
        Test that the Arrhenius n property was properly set.
        """
        self.assertAlmostEqual(self.arrhenius.n / self.n, 1.0, 6, '{0} != {1} within 6 places'.format(self.arrhenius.n, self.n))
        
    def test_Ea(self):
        """
        Test that the Arrhenius Ea property was properly set.
        """
        self.assertAlmostEqual(self.arrhenius.Ea / self.Ea, 1.0, 6, '{0} != {1} within 6 places'.format(self.arrhenius.A, self.A))
        
    def test_T0(self):
        """
        Test that the Arrhenius T0 property was properly set.
        """
        self.assertAlmostEqual(self.arrhenius.T0 / self.T0, 1.0, 6, '{0} != {1} within 6 places'.format(self.arrhenius.T0, self.T0))
        
    def test_Tmin(self):
        """
        Test that the Arrhenius Tmin property was properly set.
        """
        self.assertAlmostEqual(self.arrhenius.Tmin / self.Tmin, 1.0, 6, '{0} != {1} within 6 places'.format(self.arrhenius.Tmin, self.Tmin))
        
    def test_Tmax(self):
        """
        Test that the Arrhenius Tmax property was properly set.
        """
        self.assertAlmostEqual(self.arrhenius.Tmax / self.Tmax, 1.0, 6, '{0} != {1} within 6 places'.format(self.arrhenius.Tmax, self.Tmax))
        
    def test_order(self):
        """
        Test that the Arrhenius order property was properly set.
        """
        self.assertAlmostEqual(self.arrhenius.order / self.order, 1.0, 6, '{0} != {1} within 6 places'.format(self.arrhenius.order, self.order))
        
    def test_comment(self):
        """
        Test that the Arrhenius comment property was properly set.
        """
        self.assertEqual(self.arrhenius.comment, self.comment)
        
    def test_isTemperatureValid(self):
        """
        Test the Arrhenius.isTemperatureValid() method.
        """
        Tdata = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        validdata = numpy.array([False,True,True,True,True,True,True,True,True,True], numpy.bool)
        for T, valid in zip(Tdata, validdata):
            valid0 = self.arrhenius.isTemperatureValid(T)
            self.assertEqual(valid0, valid)
                
    def test_getRateCoefficient(self):
        """
        Test the Arrhenius.getRateCoefficient() method.
        """
        Tlist = numpy.array([200,400,600,800,1000,1200,1400,1600,1800,2000])
        kexplist = numpy.array([1.6721e2, 6.8770e7, 5.5803e9, 5.2448e10, 2.0632e11, 5.2285e11, 1.0281e12, 1.7225e12, 2.5912e12, 3.6123e12])
        for T, kexp in zip(Tlist, kexplist):
            kact = self.arrhenius.getRateCoefficient(T)
            self.assertAlmostEqual(kexp / kact, 1.0, 3, '{0} != {1} within 3 places'.format(kexp, kact))

    def test_changeT0(self):
        """
        Test the Arrhenius.changeT0() method.
        """
        Tlist = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500])
        k0list = numpy.array([self.arrhenius.getRateCoefficient(T) for T in Tlist])
        self.arrhenius.changeT0(300)
        self.assertEqual(float(self.arrhenius.T0), 300)
        for T, k0 in zip(Tlist, k0list):
            self.assertAlmostEqual(k0 / self.arrhenius.getRateCoefficient(T), 1, 6)
        
    def test_fitToData(self):
        """
        Test the Arrhenius.fitToData() method.
        """
        Tdata = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500])
        kdata = numpy.array([self.arrhenius.getRateCoefficient(T) for T in Tdata])
        arrhenius = Arrhenius().fitToData(Tdata, kdata, kunits="cm^3/(mol*s)")
        self.assertEqual(float(self.arrhenius.T0), 1)
        for T, k in zip(Tdata, kdata):
            self.assertAlmostEqual(k / arrhenius.getRateCoefficient(T), 1, 6)
        self.assertAlmostEqual(arrhenius.A / self.arrhenius.A, 1, 4)
        self.assertAlmostEqual(arrhenius.n / self.arrhenius.n, 1, 4)
        self.assertAlmostEqual(arrhenius.Ea / self.arrhenius.Ea, 1, 4)
        self.assertAlmostEqual(arrhenius.T0 / self.arrhenius.T0, 1, 4)

    def test_pickle(self):
        """
        Test that an Arrhenius object can be successfully pickled and unpickled
        with no loss of information.
        """
        import cPickle
        arrhenius = cPickle.loads(cPickle.dumps(self.arrhenius))
        self.assertEqual(self.arrhenius.A, arrhenius.A)
        self.assertEqual(self.arrhenius.n, arrhenius.n)
        self.assertEqual(self.arrhenius.Ea, arrhenius.Ea)
        self.assertEqual(self.arrhenius.T0, arrhenius.T0)
        self.assertEqual(self.arrhenius.Tmin, arrhenius.Tmin)
        self.assertEqual(self.arrhenius.Tmax, arrhenius.Tmax)
        self.assertEqual(self.arrhenius.comment, arrhenius.comment)
    
    def test_repr(self):
        """
        Test that an Arrhenius object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        exec('arrhenius = {0!r}'.format(self.arrhenius))
        self.assertEqual(self.arrhenius.A, arrhenius.A)
        self.assertEqual(self.arrhenius.n, arrhenius.n)
        self.assertEqual(self.arrhenius.Ea, arrhenius.Ea)
        self.assertEqual(self.arrhenius.T0, arrhenius.T0)
        self.assertEqual(self.arrhenius.Tmin, arrhenius.Tmin)
        self.assertEqual(self.arrhenius.Tmax, arrhenius.Tmax)
        self.assertEqual(self.arrhenius.comment, arrhenius.comment)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
