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

class TestConvertEnergy(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertEnergy` method.
    """
    
    def test_J_to_Jpermol(self):
        """
        Test the conversion of an energy from J to J/mol.
        """
        self.assertAlmostEqual(convertEnergy(1.0 / constants.Na * pq.J, pq.J / pq.mol), 1.0 * pq.J / pq.mol, 3)

    def test_Jpermol_to_J(self):
        """
        Test the conversion of an energy from J/mol to J.
        """
        self.assertAlmostEqual(convertEnergy(1.0 * pq.J / pq.mol, pq.J), 1.0 / constants.Na * pq.J, 3)

    def test_kJ_to_kcal(self):
        """
        Test the conversion of an energy from kJ to kcal.
        """
        self.assertAlmostEqual(convertEnergy(1.0 * kJ, kcal), 1.0 / 4.184 * kcal, 3)

    def test_kcal_to_kJ(self):
        """
        Test the conversion of an energy from kJ to kcal.
        """
        self.assertAlmostEqual(convertEnergy(1.0 * kcal, kJ), 4.184 * kJ, 3)

################################################################################

class TestConvertHeatCapacity(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertHeatCapacity` method.
    """
    
    def test_JperK_to_JpermolJ(self):
        """
        Test the conversion of a heat capacity from J/K to J/mol*K.
        """
        self.assertAlmostEqual(convertHeatCapacity(1.0 / constants.Na * pq.J / pq.K, pq.J / pq.mol / pq.K), 1.0 * pq.J / pq.mol / pq.K, 3)

    def test_JpermolK_to_JperK(self):
        """
        Test the conversion of a heat capacity from J/mol*K to J/K.
        """
        self.assertAlmostEqual(convertHeatCapacity(1.0 * pq.J / pq.mol / pq.K, pq.J / pq.K), 1.0 / constants.Na * pq.J / pq.K, 3)

    def test_kJperK_to_kcalperK(self):
        """
        Test the conversion of a heat capacity from kJ/K to kcal/K.
        """
        self.assertAlmostEqual(convertHeatCapacity(1.0 * kJ / pq.K, kcal / pq.K), 1.0 / 4.184 * kcal / pq.K, 3)

    def test_kcalperK_to_kJperK(self):
        """
        Test the conversion of a heat capacity from kJ/K to kcal/K.
        """
        self.assertAlmostEqual(convertHeatCapacity(1.0 * kcal / pq.K, kJ / pq.K), 4.184 * kJ / pq.K, 3)

################################################################################

class TestConvertFrequency(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertFrequency` method.
    """
    
    def test_Hz_to_GHz(self):
        """
        Test the conversion of a frequency from Hz to GHz.
        """
        self.assertAlmostEqual(convertFrequency(1.0 * pq.Hz, pq.GHz), 1.0e-9 * pq.GHz, 3)

    def test_GHz_to_Hz(self):
        """
        Test the conversion of a frequency from GHz to Hz.
        """
        self.assertAlmostEqual(convertFrequency(1.0 * pq.GHz, pq.Hz), 1.0e9 * pq.Hz, 3)

    def test_Hz_to_wavenumber(self):
        """
        Test the conversion of a frequency from Hz to cm^-1.
        """
        self.assertAlmostEqual(convertFrequency(1.0 * pq.Hz, pq.wavenumber), 1.0 / (constants.c * 100.) * pq.wavenumber, 3)

    def test_wavenumber_to_Hz(self):
        """
        Test the conversion of a frequency from cm^-1 to Hz.
        """
        self.assertAlmostEqual(convertFrequency(1.0 * pq.wavenumber, pq.Hz), (constants.c * 100.) * pq.Hz, 3)

################################################################################

class TestConvertMass(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertMass` method.
    """
    
    def test_g_to_gpermol(self):
        """
        Test the conversion of a mass from g to g/mol.
        """
        self.assertAlmostEqual(convertMass(1.0 / constants.Na * pq.g, pq.g / pq.mol), 1.0 * pq.g / pq.mol, 3)

    def test_gpermol_to_g(self):
        """
        Test the conversion of a mass from g/mol to g.
        """
        self.assertAlmostEqual(convertMass(1.0 * pq.g / pq.mol, pq.g), 1.0 / constants.Na * pq.g, 3)

    def test__to_amu(self):
        """
        Test the conversion of a mass from g to amu.
        """
        self.assertAlmostEqual(convertMass(1.0 * constants.Na * pq.u, pq.g), 1.0 * pq.g, 3)

    def test_amu_to_g(self):
        """
        Test the conversion of a mass from amu to g.
        """
        self.assertAlmostEqual(convertMass(1.0 / constants.Na * pq.g, pq.u), 1.0 * pq.u, 3)

    def test_g_to_kg(self):
        """
        Test the conversion of a mass from g to kg.
        """
        self.assertAlmostEqual(convertMass(1.0 * pq.g, pq.kg), 0.001 * pq.kg, 3)

    def test_kg_to_g(self):
        """
        Test the conversion of a mass from kg to g.
        """
        self.assertAlmostEqual(convertMass(1.0 * pq.kg, pq.g), 1000. * pq.g, 3)

################################################################################

class TestConvertInertia(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertInertia` method.
    """
    
    def test_amuperangstrom2_to_kgperm2(self):
        """
        Test the conversion of a moment of inertia from amu*angstrom^2 to 
        kg*m^2.
        """
        self.assertAlmostEqual(convertInertia(1.0 * constants.Na * 1e23 * pq.amu * pq.angstrom**2, pq.kg * pq.m**2), 1.0 * pq.kg * pq.m**2, 4)

    def test_kgperm2_to_amuperangstrom2(self):
        """
        Test the conversion of a moment of inertia from kg*m^2 to 
        amu*angstrom^2.
        """
        self.assertAlmostEqual(convertInertia(1.0 / constants.Na / 1e23 * pq.kg * pq.m**2, pq.amu * pq.angstrom**2), 1.0 * pq.amu * pq.angstrom**2, 4)

################################################################################

class TestConvertRateCoefficient(unittest.TestCase):
    """
    Contains unit tests of the :func:`convertRateCoefficient` method.
    """
    
    def test_cm3permols_to_m3permols(self):
        """
        Test the conversion of a rate coefficient from cm^3/(mol*s) to 
        m^3/(mol*s).
        """
        self.assertAlmostEqual(convertRateCoefficient(1.0e6 * pq.cm**3/pq.mol/pq.s, pq.m**3/pq.mol/pq.s), 1.0 * pq.m**3/pq.mol/pq.s, 4)

    def test_m3permols_to_cm3permols(self):
        """
        Test the conversion of a rate coefficient from m^3/(mol*s) to 
        cm^3/(mol*s).
        """
        self.assertAlmostEqual(convertRateCoefficient(1.0e-6 * pq.m**3/pq.mol/pq.s, pq.cm**3/pq.mol/pq.s), 1.0 * pq.cm**3/pq.mol/pq.s, 4)

    def test_cm3permols_to_cm3permolecules(self):
        """
        Test the conversion of a rate coefficient from cm^3/(mol*s) to 
        cm^3/(molecule*s).
        """
        self.assertAlmostEqual(convertRateCoefficient(constants.Na * pq.cm**3/pq.mol/pq.s, pq.cm**3/molecule/pq.s), 1.0 * pq.cm**3/molecule/pq.s, 4)

    def test_cm3permolecules_to_cm3permols(self):
        """
        Test the conversion of a rate coefficient from cm^3/(molecule*s) to 
        cm^3/(mol*s).
        """
        self.assertAlmostEqual(convertRateCoefficient(1.0 / constants.Na * pq.cm**3/molecule/pq.s, pq.cm**3/pq.mol/pq.s), 1.0 * pq.cm**3/pq.mol/pq.s, 4)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
