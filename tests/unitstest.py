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
This script contains unit tests of the :mod:`pyrate.units` module.
"""

import unittest
import math
import quantities as pq

from pyrate.units import *

################################################################################

class TestConvertTemperature(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertTemperature` method.
    """
    
    def test_degF_to_degC(self):
        """
        Test the conversion of a temperature from deg F to deg C.
        """
        self.assertAlmostEqual(convertTemperature(68 * pq.degF, pq.degC), 20 * pq.degC, 3)

    def test_degF_to_degR(self):
        """
        Test the conversion of a temperature from deg F to deg R.
        """
        self.assertAlmostEqual(convertTemperature(68 * pq.degF, pq.degR), 527.67 * pq.degR, 3)

    def test_degF_to_K(self):
        """
        Test the conversion of a temperature from deg F to K.
        """
        self.assertAlmostEqual(convertTemperature(68 * pq.degF, pq.K), 293.15 * pq.K, 3)

    def test_degC_to_degF(self):
        """
        Test the conversion of a temperature from deg C to deg F.
        """
        self.assertAlmostEqual(convertTemperature(20 * pq.degC, pq.degF), 68 * pq.degF, 3)

    def test_degC_to_degR(self):
        """
        Test the conversion of a temperature from deg C to deg R.
        """
        self.assertAlmostEqual(convertTemperature(20 * pq.degC, pq.degR), 527.67 * pq.degR, 3)

    def test_degC_to_K(self):
        """
        Test the conversion of a temperature from deg C to K.
        """
        self.assertAlmostEqual(convertTemperature(20 * pq.degC, pq.K), 293.15 * pq.K, 3)

    def test_degR_to_degF(self):
        """
        Test the conversion of a temperature from deg R to deg F.
        """
        self.assertAlmostEqual(convertTemperature(527.67 * pq.degR, pq.degF), 68 * pq.degF, 3)

    def test_degR_to_degC(self):
        """
        Test the conversion of a temperature from deg R to deg C.
        """
        self.assertAlmostEqual(convertTemperature(527.67 * pq.degR, pq.degC), 20 * pq.degC, 3)

    def test_degR_to_K(self):
        """
        Test the conversion of a temperature from deg R to K.
        """
        self.assertAlmostEqual(convertTemperature(527.67 * pq.degR, pq.K), 293.15 * pq.K, 3)

    def test_K_to_degF(self):
        """
        Test the conversion of a temperature from K to deg F.
        """
        self.assertAlmostEqual(convertTemperature(293.15 * pq.K, pq.degF), 68 * pq.degF, 3)

    def test_K_to_degC(self):
        """
        Test the conversion of a temperature from K to deg C.
        """
        self.assertAlmostEqual(convertTemperature(293.15 * pq.K, pq.degC), 20 * pq.degC, 3)

    def test_K_to_degR(self):
        """
        Test the conversion of a temperature from K to deg K.
        """
        self.assertAlmostEqual(convertTemperature(293.15 * pq.K, pq.degR), 527.67 * pq.degR, 3)

################################################################################

class TestConvertPressure(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertPressure` method.
    """
    
    def test_atm_to_bar(self):
        """
        Test the conversion of a pressure from atm to bar.
        """
        self.assertAlmostEqual(convertPressure(1.0 * pq.atm, pq.bar), 1.01325 * pq.bar, 4)

    def test_Pa_to_bar(self):
        """
        Test the conversion of a pressure from Pa to bar.
        """
        self.assertAlmostEqual(convertPressure(100000. * pq.Pa, pq.bar), 1.0 * pq.bar, 4)

    def test_torr_to_bar(self):
        """
        Test the conversion of a pressure from torr to bar.
        """
        self.assertAlmostEqual(convertPressure(760. * pq.torr, pq.bar), 1.01325  * pq.bar, 4)

    def test_psi_to_bar(self):
        """
        Test the conversion of a pressure from Pa to bar.
        """
        self.assertAlmostEqual(convertPressure(14.7 * pq.psi, pq.bar), 1.01325  * pq.bar, 3)

    def test_kPa_to_bar(self):
        """
        Test the conversion of a pressure from kPa to bar.
        """
        self.assertAlmostEqual(convertPressure(100. * pq.kPa, pq.bar), 1.0 * pq.bar, 4)

    def test_MPa_to_bar(self):
        """
        Test the conversion of a pressure from MPa to bar.
        """
        self.assertAlmostEqual(convertPressure(0.1 * pq.MPa, pq.bar), 1.0 * pq.bar, 4)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
