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

class TestPDepArrhenius(unittest.TestCase):
    """
    Contains unit tests of the :class:`PDepArrhenius` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.arrhenius0 = Arrhenius(
            A = (1.0e6,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (10,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        self.arrhenius1 = Arrhenius(
            A = (1.0e12,"cm^3/(mol*s)"), 
            n = 1.0, 
            Ea = (20,"kJ/mol"), 
            T0 = (300,"K"), 
            Tmin = (300.0,"K"), 
            Tmax = (2000.0,"K"), 
            comment = """This data is completely made up""",
        )
        self.pressures = pq.Quantity([0.1, 10.0],"bar")
        self.arrhenius = [self.arrhenius0, self.arrhenius1]
        self.Tmin = pq.Quantity(300.0,"K")
        self.Tmax = pq.Quantity(2000.0,"K")
        self.Pmin = pq.Quantity(0.1,"bar")
        self.Pmax = pq.Quantity(10.0,"bar")
        self.order = 2
        self.comment = """This data is completely made up"""
        self.kinetics = PDepArrhenius(
            pressures = self.pressures,
            arrhenius = self.arrhenius,
            Tmin = self.Tmin, 
            Tmax = self.Tmax, 
            Pmin = self.Pmin, 
            Pmax = self.Pmax,
            order = self.order,
            comment = self.comment,
        )

    def test_pressures(self):
        """
        Test that the PDepArrhenius pressures property was properly set.
        """
        self.assertEqual(len(self.kinetics.pressures), 2)
        for i in range(2):
            self.assertEqual(self.kinetics.pressures[i], self.pressures[i])
        
    def test_arrhenius(self):
        """
        Test that the PDepArrhenius arrhenius property was properly set.
        """
        self.assertEqual(len(self.kinetics.arrhenius), 2)
        for i in range(2):
            self.assertEqual(self.kinetics.arrhenius[i].A, self.arrhenius[i].A)
            self.assertEqual(self.kinetics.arrhenius[i].n, self.arrhenius[i].n)
            self.assertEqual(self.kinetics.arrhenius[i].T0, self.arrhenius[i].T0)
            self.assertEqual(self.kinetics.arrhenius[i].Ea, self.arrhenius[i].Ea)
            self.assertEqual(self.kinetics.arrhenius[i].Tmin, self.arrhenius[i].Tmin)
            self.assertEqual(self.kinetics.arrhenius[i].Tmax, self.arrhenius[i].Tmax)
            self.assertEqual(self.kinetics.arrhenius[i].order, self.arrhenius[i].order)
            self.assertEqual(self.kinetics.arrhenius[i].comment, self.arrhenius[i].comment)
        
    def test_Tmin(self):
        """
        Test that the PDepArrhenius Tmin property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Tmin / self.Tmin, 1.0, 6, '{0} != {1} within 6 places'.format(self.kinetics.Tmin, self.Tmin))
        
    def test_Tmax(self):
        """
        Test that the PDepArrhenius Tmax property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Tmax / self.Tmax, 1.0, 6, '{0} != {1} within 6 places'.format(self.kinetics.Tmax, self.Tmax))

    def test_Pmin(self):
        """
        Test that the PDepArrhenius Pmin property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Pmin / self.Pmin, 1.0, 6, '{0} != {1} within 6 places'.format(self.kinetics.Pmin, self.Pmin))
        
    def test_Pmax(self):
        """
        Test that the PDepArrhenius Pmax property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.Pmax / self.Pmax, 1.0, 6, '{0} != {1} within 6 places'.format(self.kinetics.Pmax, self.Pmax))

    def test_order(self):
        """
        Test that the PDepArrhenius order property was properly set.
        """
        self.assertAlmostEqual(self.kinetics.order / self.order, 1.0, 6, '{0} != {1} within 6 places'.format(self.kinetics.order, self.order))
        
    def test_comment(self):
        """
        Test that the PDepArrhenius comment property was properly set.
        """
        self.assertEqual(self.kinetics.comment, self.comment)

    def test_isPressureDependent(self):
        """
        Test the PDepArrhenius.isPressureDependent() method.
        """
        self.assertTrue(self.kinetics.isPressureDependent())
    
    def test_getRateCoefficient(self):
        """
        Test the PDepArrhenius.getRateCoefficient() method.
        """
        print 
        P = 1e-1
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            k0 = self.kinetics.getRateCoefficient(T, P)
            k1 = self.arrhenius0.getRateCoefficient(T)
            self.assertAlmostEqual(k0 / k1, 1, 6, '{0} != {1} within 6 places'.format(k0, k1))
        P = 1e1
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            k0 = self.kinetics.getRateCoefficient(T, P)
            k1 = self.arrhenius1.getRateCoefficient(T)
            self.assertAlmostEqual(k0 / k1, 1, 6, '{0} != {1} within 6 places'.format(k0, k1))
        P = 1e0
        for T in [300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500]:
            k0 = self.kinetics.getRateCoefficient(T, P)
            k1 = math.sqrt(self.arrhenius0.getRateCoefficient(T) * self.arrhenius1.getRateCoefficient(T))
            self.assertAlmostEqual(k0 / k1, 1, 6, '{0} != {1} within 6 places'.format(k0, k1))
        
    def test_fitToData(self):
        """
        Test the PDepArrhenius.fitToData() method.
        """
        Tdata = numpy.array([300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500], numpy.float)
        Pdata = numpy.array([1e-1,3e-1,1e0,3e0,1e1], numpy.float)
        kdata = numpy.zeros([len(Tdata),len(Pdata)], numpy.float)
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                kdata[t,p] = self.kinetics.getRateCoefficient(Tdata[t], Pdata[p])
        kinetics = PDepArrhenius().fitToData(Tdata, Pdata, kdata, kunits='cm^3/(mol*s)')
        for t in range(len(Tdata)):
            for p in range(len(Pdata)):
                self.assertAlmostEqual(kinetics.getRateCoefficient(Tdata[t], Pdata[p]) / kdata[t,p], 1, 4)
        
    def test_pickle(self):
        """
        Test that a PDepArrhenius object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        kinetics = cPickle.loads(cPickle.dumps(self.kinetics))
        Narrh = 2
        self.assertEqual(len(self.kinetics.pressures), Narrh)
        self.assertEqual(len(kinetics.pressures), Narrh)
        self.assertEqual(len(self.kinetics.arrhenius), Narrh)
        self.assertEqual(len(kinetics.arrhenius), Narrh)
        for i in range(Narrh):
            self.assertEqual(self.kinetics.pressures[i], kinetics.pressures[i])
            self.assertEqual(self.kinetics.arrhenius[i].A, kinetics.arrhenius[i].A)
            self.assertEqual(self.kinetics.arrhenius[i].n, kinetics.arrhenius[i].n)
            self.assertEqual(self.kinetics.arrhenius[i].T0, kinetics.arrhenius[i].T0)
            self.assertEqual(self.kinetics.arrhenius[i].Ea, kinetics.arrhenius[i].Ea)
        self.assertEqual(self.kinetics.Tmin, kinetics.Tmin)
        self.assertEqual(self.kinetics.Tmax, kinetics.Tmax)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
        
    def test_repr(self):
        """
        Test that a PDepArrhenius object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('kinetics = {0!r}'.format(self.kinetics))
        Narrh = 2
        self.assertEqual(len(self.kinetics.pressures), Narrh)
        self.assertEqual(len(kinetics.pressures), Narrh)
        self.assertEqual(len(self.kinetics.arrhenius), Narrh)
        self.assertEqual(len(kinetics.arrhenius), Narrh)
        for i in range(Narrh):
            self.assertEqual(self.kinetics.pressures[i], kinetics.pressures[i])
            self.assertEqual(self.kinetics.arrhenius[i].A, kinetics.arrhenius[i].A)
            self.assertEqual(self.kinetics.arrhenius[i].n, kinetics.arrhenius[i].n)
            self.assertEqual(self.kinetics.arrhenius[i].T0, kinetics.arrhenius[i].T0)
            self.assertEqual(self.kinetics.arrhenius[i].Ea, kinetics.arrhenius[i].Ea)
        self.assertEqual(self.kinetics.Tmin, kinetics.Tmin)
        self.assertEqual(self.kinetics.Tmax, kinetics.Tmax)
        self.assertEqual(self.kinetics.Pmin, kinetics.Pmin)
        self.assertEqual(self.kinetics.Pmax, kinetics.Pmax)
        self.assertEqual(self.kinetics.comment, kinetics.comment)
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
