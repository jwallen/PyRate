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
This module contains functions for unit conversions of various types of 
physical quantities, which provide error checking and enable some additional
unit conversions not supported by the :mod:`quantities` package.
"""

import quantities as pq

import pyrate.constants as constants

################################################################################

# Explicity set the default units to SI
pq.set_default_units('si')

# These units are not defined by the quantities package, but occur frequently
# in data handled by PyRate, so we define them manually
kcal = pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
kJ = pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
kmol = pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol')
molecule = pq.UnitQuantity('molecule', pq.mol/constants.Na, symbol='molecule')
molecules = pq.UnitQuantity('molecules', pq.mol/constants.Na, symbol='molecules')

################################################################################

# Allowed simplified temperature dimensionality
TEMPERATURE_DIMENSIONS = [pq.K.simplified.dimensionality]

def convertTemperature(quantity, units):
    """
    Convert a given `quantity` with units of temperature to the given `units` 
    of temperature. A :class:`ValueError` is raised if this conversion is not
    successful.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a Quantity object with units of temperature.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    
    inputDimensionality = quantity.units.dimensionality
    outputDimensionality = units.dimensionality

    if inputDimensionality.simplified not in TEMPERATURE_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for temperature.'.format(quantity.units))
    if outputDimensionality.simplified not in TEMPERATURE_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for temperature.'.format(units))
    
    if inputDimensionality == pq.degF.dimensionality:
        quantity = quantity + 459.67 * pq.degF
    elif inputDimensionality == pq.degC.dimensionality:
        quantity = quantity + 273.15 * pq.degC
    
    quantity = quantity.rescale(units)
    
    if outputDimensionality == pq.degF.dimensionality:
        quantity = quantity - 459.67 * pq.degF
    elif outputDimensionality == pq.degC.dimensionality:
        quantity = quantity - 273.15 * pq.degC
        
    return quantity

################################################################################

# Allowed pressure units
PRESSURE_DIMENSIONS = [pq.Pa.simplified.dimensionality]

def convertPressure(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    pressure -- to the given `units` of pressure. A :class:`ValueError`
    is raised if this conversion is not successful.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a Quantity object with units of temperature.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    
    inputDimensionality = quantity.units.dimensionality
    outputDimensionality = units.dimensionality

    if inputDimensionality.simplified not in PRESSURE_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for pressure.'.format(quantity.units))
    if outputDimensionality.simplified not in PRESSURE_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for pressure.'.format(units))
        
    return quantity.rescale(units)

################################################################################

# Allowed energy units
ENERGY_DIMENSIONS = [pq.J.simplified.dimensionality]
MOLAR_ENERGY_DIMENSIONS = [(pq.J / pq.mol).simplified.dimensionality]

def convertEnergy(quantity, units, unittype='energy'):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    energy -- to the given `units` of energy. A :class:`ValueError` is raised
    if this conversion is not successful. This function can handle conversion
    between intensive (molar) and extensive energies by multiplying or dividing
    by an Avogadro number as appropriate.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a Quantity object with units of {1}.'.format(quantity, unittype))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    inputDimensionality = quantity.units.dimensionality.simplified
    outputDimensionality = units.dimensionality.simplified
    
    if inputDimensionality not in ENERGY_DIMENSIONS and inputDimensionality not in MOLAR_ENERGY_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for {1}.'.format(quantity.units, unittype))
    if outputDimensionality not in ENERGY_DIMENSIONS and outputDimensionality not in MOLAR_ENERGY_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for {1}.'.format(units, unittype))
    
    if inputDimensionality in ENERGY_DIMENSIONS and outputDimensionality in MOLAR_ENERGY_DIMENSIONS:
        quantity = quantity * constants.Na / pq.mol
    elif inputDimensionality in MOLAR_ENERGY_DIMENSIONS and outputDimensionality in ENERGY_DIMENSIONS:
        quantity = quantity / (constants.Na / pq.mol)
    
    return quantity.rescale(units)

def convertEnthalpy(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    enthalpy -- to the given `units` of enthalpy. A :class:`ValueError` is raised
    if this conversion is not successful. This function can handle conversion
    between intensive (molar) and extensive enthalpies by multiplying or
    dividing by an Avogadro number as appropriate.
    """
    return convertEnergy(quantity, units, 'enthalpy')

def convertFreeEnergy(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    free energy -- to the given `units` of free energy. A :class:`ValueError` 
    is raised if this conversion is not successful. This function can handle 
    conversion between intensive (molar) and extensive free energies by
    multiplying or dividing by an Avogadro number as appropriate.
    """
    return convertEnergy(quantity, units, 'free energy')

################################################################################

# Allowed heat capacity units
HEATCAPACITY_DIMENSIONS = [(pq.J / pq.K).simplified.dimensionality]
MOLAR_HEATCAPACITY_DIMENSIONS = [(pq.J /pq.mol / pq.K).simplified.dimensionality]

def convertHeatCapacity(quantity, units, unittype='heat capacity'):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    heat capacity -- to the given `units` of heat capacity. A 
    :class:`ValueError` is raised if this conversion is not successful. This
    function can handle conversion between intensive (molar) and extensive heat
    capacities by multiplying or dividing by an Avogadro number as appropriate.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a Quantity object with units of {1}.'.format(quantity, unittype))

    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    inputDimensionality = quantity.units.dimensionality.simplified
    outputDimensionality = units.dimensionality.simplified

    if inputDimensionality not in HEATCAPACITY_DIMENSIONS and inputDimensionality not in MOLAR_HEATCAPACITY_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for {1}.'.format(quantity.units, unittype))
    if outputDimensionality not in HEATCAPACITY_DIMENSIONS and outputDimensionality not in MOLAR_HEATCAPACITY_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for {1}.'.format(units, unittype))
    
    if inputDimensionality in HEATCAPACITY_DIMENSIONS and outputDimensionality in MOLAR_HEATCAPACITY_DIMENSIONS:
        quantity = quantity * constants.Na / pq.mol
    elif inputDimensionality in MOLAR_HEATCAPACITY_DIMENSIONS and outputDimensionality in HEATCAPACITY_DIMENSIONS:
        quantity = quantity / (constants.Na / pq.mol)
    
    return quantity.rescale(units)

def convertEntropy(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    entropy -- to the given `units` of entropy. A :class:`ValueError` is 
    raised if this conversion is not successful. This function can handle 
    conversion between intensive (molar) and extensive entropies by
    multiplying or dividing by an Avogadro number as appropriate.
    """
    return convertHeatCapacity(quantity, units, 'entropy')

################################################################################

# Allowed simplified frequency dimensionality
FREQUENCY_DIMENSIONS = [(1./pq.s.simplified).dimensionality, pq.Hz.simplified.dimensionality]
EXTRA_FREQUENCY_DIMENSIONS = {
    pq.J.dimensionality: 1.0 / (constants.h * pq.J * pq.s),
    pq.wavenumber.dimensionality: constants.c * 100. / pq.s / pq.wavenumber,
    pq.eV.dimensionality: constants.e / (constants.h * pq.eV * pq.s),
    (1.0/pq.cm).dimensionality: constants.c * 100. / pq.s * pq.cm,
}

def convertFrequency(quantity, units):
    """
    Convert a given `quantity` with units of frequency to the given `units` 
    of frequency. A :class:`ValueError` is raised if this conversion is not
    successful.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a Quantity object with units of frequency.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    
    inputDimensionality = quantity.units.dimensionality
    outputDimensionality = units.dimensionality
    
    if inputDimensionality.simplified not in FREQUENCY_DIMENSIONS and inputDimensionality not in EXTRA_FREQUENCY_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for frequency.'.format(quantity.units))
    if outputDimensionality.simplified not in FREQUENCY_DIMENSIONS and outputDimensionality not in EXTRA_FREQUENCY_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for frequency.'.format(units))
    
    if inputDimensionality in EXTRA_FREQUENCY_DIMENSIONS:
        quantity = quantity * EXTRA_FREQUENCY_DIMENSIONS[inputDimensionality]
    else:
        quantity = quantity.rescale(pq.Hz)
    
    if outputDimensionality in EXTRA_FREQUENCY_DIMENSIONS:
        return quantity / EXTRA_FREQUENCY_DIMENSIONS[outputDimensionality]
    else:
        return quantity.rescale(units)

################################################################################

# Allowed mass units
MASS_DIMENSIONS = [pq.kg.simplified.dimensionality]
MOLAR_MASS_DIMENSIONS = [(pq.kg / pq.mol).simplified.dimensionality]

def convertMass(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    mass -- to the given `units` of mass. A :class:`ValueError` is raised
    if this conversion is not successful. This function can handle conversion
    between intensive (molar) and extensive masses by multiplying or dividing
    by an Avogadro number as appropriate.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a Quantity object with units of mass.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
        
    inputDimensionality = quantity.units.dimensionality.simplified
    outputDimensionality = units.dimensionality.simplified
    
    if inputDimensionality not in MASS_DIMENSIONS and inputDimensionality not in MOLAR_MASS_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for {1}.'.format(quantity.units, unittype))
    if outputDimensionality not in MASS_DIMENSIONS and outputDimensionality not in MOLAR_MASS_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for {1}.'.format(units, unittype))
    
    if inputDimensionality in MASS_DIMENSIONS and outputDimensionality in MOLAR_MASS_DIMENSIONS:
        quantity = quantity * constants.Na / pq.mol
    elif inputDimensionality in MOLAR_MASS_DIMENSIONS and outputDimensionality in MASS_DIMENSIONS:
        quantity = quantity / (constants.Na / pq.mol)
    
    return quantity.rescale(units)

################################################################################

INERTIA_DIMENSIONS = [(pq.kg * pq.m * pq.m).simplified.dimensionality]

def convertInertia(quantity, units):
    """
    Convert a given `quantity` with units of moment of inertia to the given
    `units` of moment of inertia. A :class:`ValueError` is raised if this 
    conversion is not successful.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a Quantity object with units of moment of inertia.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    
    inputDimensionality = quantity.units.dimensionality
    outputDimensionality = units.dimensionality
    
    if inputDimensionality.simplified not in INERTIA_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for moment of inertia.'.format(quantity.units))
    if outputDimensionality.simplified not in INERTIA_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for moment of inertia.'.format(units))
    
    return quantity.rescale(units)

################################################################################

# Note: There are many possible units for rate coefficients, but they should
# all produce a combination of m^3, mol, and s when simplified
RATE_COEFFICIENT_DIMENSIONS = [
    (1.0/pq.s).simplified.dimensionality, 
    (pq.m**3/(pq.mol*pq.s)).simplified.dimensionality, 
    (pq.m**6/(pq.mol**2*pq.s)).simplified.dimensionality, 
    (pq.m**9/(pq.mol**3*pq.s)).simplified.dimensionality,
]

def convertRateCoefficient(quantity, units):
    """
    Convert a given `quantity` -- a :class:`Quantity` object with units of 
    rate coefficient -- to the given `units` of rate coefficient. A 
    :class:`ValueError` is raised if this conversion is not successful.
    Conversion between molar and molecular bases are possible, as long as you
    explicitly specify the units of molecules.
    """
    if isinstance(quantity, tuple):
        quantity = pq.Quantity(quantity[0], quantity[1])
    elif not isinstance(quantity, pq.Quantity):
        raise ValueError('Invalid value "{0}" for quantity; must be a Quantity object with units of mass.'.format(quantity))
    
    if isinstance(units, str):
        units = pq.Quantity(1.0, units)
    inputDimensionality = quantity.units.dimensionality.simplified
    outputDimensionality = units.dimensionality.simplified
    
    if inputDimensionality not in RATE_COEFFICIENT_DIMENSIONS:
        raise ValueError('Invalid input units "{0}" for {1}.'.format(quantity.units, unittype))
    if outputDimensionality not in RATE_COEFFICIENT_DIMENSIONS:
        raise ValueError('Invalid output units "{0}" for {1}.'.format(units, unittype))
    
    return quantity.rescale(units)
