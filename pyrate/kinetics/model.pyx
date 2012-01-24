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
This module contains Cython extension types that represent various models of
chemical reaction kinetics. Both pressure-independent and pressure-dependent 
models are available.
"""

import numpy
import quantities as pq
from libc.math cimport exp, log, sqrt, log10, cos, acos

cimport pyrate.constants as constants
import pyrate.units as units

################################################################################

cpdef getRateCoefficientUnitsFromReactionOrder(order):
    """
    Given a reaction `order`, return the corresponding SI units of the rate
    coefficient. These are the units that rate coefficients are stored in
    internally, as well as the units of the rate coefficient obtained using
    the ``simplified`` attribute of a :class:`Quantity` object that represents
    a rate coefficient. Raises a :class:`ValueError` if the units could not be
    determined.
    """
    if order == 0: 
        kunits = 'mol/(cm^3*s)'
    elif order == 1:
        kunits = 's^-1'
    elif order == 2:
        kunits = 'cm^3/(mol*s)'
    elif order == 3:
        kunits = 'cm^6/(mol^2*s)'
    elif order == 4:
        kunits = 'cm^9/(mol^3*s)'
    else:
        raise ValueError('Invalid reaction order {0}.'.format(order))
    return kunits

cpdef getReactionOrderFromRateCoefficientUnits(kunits):
    """
    Given a set of rate coefficient units `kunits`, return the corresponding
    reaction order. Raises a :class:`ValueError` if the reaction order could 
    not be determined.
    """
    dimensionality = pq.Quantity(1.0, kunits).simplified.dimensionality
    order = -1
    if len(dimensionality) == 3 and pq.mol in dimensionality and pq.m in dimensionality and pq.s in dimensionality:
        if dimensionality[pq.s] == -1 and dimensionality[pq.m] == -3 and dimensionality[pq.mol] == 1:
            order = 0
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 3 and dimensionality[pq.mol] == -1:
            order = 2
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 6 and dimensionality[pq.mol] == -2:
            order = 3
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 9 and dimensionality[pq.mol] == -3:
            order = 4
    elif len(dimensionality) == 1 and pq.s in dimensionality:
        if dimensionality[pq.s] == -1:
            order = 1
    if order == -1:
        raise ValueError('Invalid rate coefficient units "{0}".'.format(str(dimensionality)))
    return order

################################################################################

cdef class KineticsModel:
    """
    A base class for kinetics models, containing several attributes
    common to all models:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `order`         The reaction order (usually the number of reactants) 
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, Tmin=None, Tmax=None, Pmin=None, Pmax=None, order=-1, comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.order = order
        self.comment = comment
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KineticsModel object.
        """
        Tmin = '({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        Tmax = '({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        Pmin = '({0:g},"{1}")'.format(float(self.Pmin), str(self.Pmin.dimensionality))
        Pmax = '({0:g},"{1}")'.format(float(self.Pmax), str(self.Pmax.dimensionality))
        return 'KineticsModel(Tmin={0}, Tmax={1}, Pmin={2}, Pmax={3}, order={4:d}, comment="""{4}""")'.format(Tmin, Tmax, Pmin, Pmax, self.order, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsModel object.
        """
        return (KineticsModel, (self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.order, self.comment))

    property Tmin:
        """The minimum temperature at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Tmin, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Tmin = 0.0 
            else:
                self._Tmin = float(units.convertTemperature(value, pq.K))

    property Tmax:
        """The maximum temperature at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Tmax, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Tmax = 0.0 
            else:
                self._Tmax = float(units.convertTemperature(value, pq.K))

    property Pmin:
        """The minimum pressure at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Pmin * 1e-5, pq.bar)
        def __set__(self, value):
            if value is None or value == 0:
                self._Pmin = 0.0 
            else:
                self._Pmin = float(units.convertPressure(value, pq.Pa))

    property Pmax:
        """The maximum pressure at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Pmax * 1e-5, pq.bar)
        def __set__(self, value):
            if value is None or value == 0:
                self._Pmax = 0.0 
            else:
                self._Pmax = float(units.convertPressure(value, pq.Pa))

    cpdef bint isTemperatureValid(self, double T) except -2:
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the temperature model, or ``False`` if not. If
        the minimum and maximum temperature are not defined, ``True`` is 
        returned.
        """
        return (self._Tmin == 0.0 or self._Tmin <= T) and (self._Tmax == 0.0 or T <= self._Tmax)

    cpdef bint isPressureValid(self, double P) except -2:
        """
        Return ``True`` if the pressure `P` in Pa is within the valid pressure
        range of the kinetics model, or ``False`` if not. If the minimum and 
        maximum pressure are not defined, ``True`` is returned.
        """
        return (self._Pmin == 0.0 or self._Pmin <= P) and (self._Pmax == 0.0 or P <= self._Pmax)

    cpdef bint isPressureDependent(self) except -2:
        """
        Returns ``True`` if the kinetics model represents pressure-dependent
        kinetics or ``False`` otherwise. This method must be overloaded in the 
        derived class.
        """
        raise NotImplementedError('Unexpected call to KineticsModel.isPressureDependent(); you should be using a class derived from KineticsModel that overrides this method.')

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of cm^3, 
        mol, and s at temperature `T` in K and, if required, pressure `P` in 
        bar. This method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to KineticsModel.getRateCoefficient(); you should be using a class derived from KineticsModel that overrides this method.')
