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
This module contains Cython extension types that represent various
thermodynamic models of heat capacity.
"""

import numpy
import quantities as pq
from libc.math cimport log

cimport pyrate.constants as constants
import pyrate.units as units

################################################################################

cdef class HeatCapacityModel:
    """
    A base class for heat capacity models, containing several attributes
    common to all models:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, Tmin=None, Tmax=None, comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.comment = comment
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ThermoModel object.
        """
        Tmin = '({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        Tmax = '({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        return 'ThermoModel(Tmin={0}, Tmax={1}, comment="""{2}""")'.format(Tmin, Tmax, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a ThermoModel object.
        """
        return (HeatCapacityModel, (self.Tmin, self.Tmax, self.comment))

    property Tmin:
        """The minimum temperature at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Tmin, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Tmin = 0.0 
            else:
                self._Tmin = float(units.convertTemperature(value, "K"))

    property Tmax:
        """The maximum temperature at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Tmax, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Tmax = 0.0 
            else:
                self._Tmax = float(units.convertTemperature(value, "K"))

    cpdef bint isTemperatureValid(self, double T) except -2:
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the thermodynamic data, or ``False`` if not. If
        the minimum and maximum temperature are not defined, ``True`` is 
        returned.
        """
        return (self._Tmin == 0.0 or self._Tmin <= T) and (self._Tmax == 0.0 or T <= self._Tmax)

    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at temperature 
        `T` in K. This method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to ThermoModel.get_heat_capacity(); you should be using a class derived from ThermoModel.')

    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in J/mol at temperature `T` in K. This method must 
        be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to ThermoModel.get_enthalpy(); you should be using a class derived from ThermoModel.')

    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at temperature `T` in K. This method must 
        be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to ThermoModel.get_entropy(); you should be using a class derived from ThermoModel.')

    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in J/mol at temperature `T` in K. This 
        method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to ThermoModel.get_free_energy(); you should be using a class derived from ThermoModel.')
