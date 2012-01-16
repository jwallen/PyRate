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
This script contains unit tests of the :mod:`pyrate.thermo.converter` module.
"""

import unittest
import math
import numpy
import quantities as pq

from pyrate.thermo.model import *
from pyrate.thermo.converter import *
import pyrate.constants as constants

################################################################################

class TestConverter(unittest.TestCase):
    """
    Contains unit tests of the thermodynamics model conversion functions.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.wilhoit = Wilhoit(
            Cp0 = pq.Quantity(33.2579, "J/(mol*K)"),
            CpInf = pq.Quantity(178.761, "J/(mol*K)"),
            a0 = 0.0977518,
            a1 = -16.3067,
            a2 = 26.2524,
            a3 = -12.6785,
            B = pq.Quantity(1068.68, "K"),
            H0 = pq.Quantity(-782.292, "kJ/mol"), 
            S0 = pq.Quantity(-984.932, "J/(mol*K)"),
            Tmin = 10 * pq.K,
            Tmax = 3000 * pq.K,
            comment = 'C2H6',
        )
        self.multiNASA = MultiNASA(
            polynomials = [
                NASA(Tmin=10 * pq.K, Tmax=650.73 * pq.K, coeffs=[4.03055,-0.00214171,4.90611e-05,-5.99027e-08,2.38945e-11,-11257.6,3.5613], comment="""Low temperature range polynomial"""),
                NASA(Tmin=650.73 * pq.K, Tmax=3000 * pq.K, coeffs=[-0.307954,0.0245269,-1.2413e-05,3.07724e-09,-3.01467e-13,-10693,22.628], comment="""High temperature range polynomial"""),
            ],
            Tmin = 10 * pq.K,
            Tmax = 3000 * pq.K,
            comment = 'C2H6',
        )
    
    def test_convert_Wilhoit_to_MultiNASA(self):
        """
        Test the conversion of a Wilhoit model to a MultiNASA model.
        """
        wilhoit = self.wilhoit
        multiNASA = convertWilhoitToMultiNASA(wilhoit, Tmin=10, Tmax=3000, Tint=1000)
        Tlist = numpy.arange(10, 3000, 10)
        for T in Tlist:
            Cp_wilhoit = wilhoit.getHeatCapacity(T)
            Cp_nasa = multiNASA.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_nasa / Cp_wilhoit, 1.0, 2, '{0} != {1} within 2 places'.format(Cp_nasa, Cp_wilhoit))
            H_wilhoit = wilhoit.getEnthalpy(T)
            H_nasa = multiNASA.getEnthalpy(T)
            self.assertAlmostEqual(H_nasa / H_wilhoit, 1.0, 1, '{0} != {1} within 1 place'.format(H_nasa, H_wilhoit))
            S_wilhoit = wilhoit.getEntropy(T)
            S_nasa = multiNASA.getEntropy(T)
            self.assertAlmostEqual(S_nasa / S_wilhoit, 1.0, 2, '{0} != {1} within 2 places'.format(S_nasa, S_wilhoit))

    def test_convert_MultiNASA_to_Wilhoit(self):
        """
        Test the conversion of a MultiNASA model to a Wilhoit model.
        """
        multiNASA = self.multiNASA
        wilhoit = convertMultiNASAToWilhoit(multiNASA, linear=False, Nfreq=17, Nrotors=1)
        Tlist = numpy.arange(10, 3000, 10)
        for T in Tlist:
            Cp_wilhoit = wilhoit.getHeatCapacity(T)
            Cp_nasa = multiNASA.getHeatCapacity(T)
            self.assertAlmostEqual(Cp_nasa / Cp_wilhoit, 1.0, 2, '{0} != {1} within 2 places'.format(Cp_nasa, Cp_wilhoit))
            H_wilhoit = wilhoit.getEnthalpy(T)
            H_nasa = multiNASA.getEnthalpy(T)
            self.assertAlmostEqual(H_nasa / H_wilhoit, 1.0, 1, '{0} != {1} within 1 place'.format(H_nasa, H_wilhoit))
            S_wilhoit = wilhoit.getEntropy(T)
            S_nasa = multiNASA.getEntropy(T)
            self.assertAlmostEqual(S_nasa / S_wilhoit, 1.0, 2, '{0} != {1} within 2 places'.format(S_nasa, S_wilhoit))
        
################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
