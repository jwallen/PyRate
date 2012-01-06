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
