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
This module contains classes that represent various models of translational
motion. Translational energies are much smaller than :math:`k_\\mathrm{B} T`
except for temperatures approaching absolute zero, so a classical treatment of
translation is more than adequate.
"""

import math
import numpy
import quantities as pq
from libc.math cimport log, sqrt

cimport pyrate.constants as constants
cimport pyrate.statmech.schrodinger as schrodinger 
import pyrate.units as units

################################################################################

cdef class IdealGasTranslation:
    """
    A statistical mechanical model of translation in an 3-dimensional infinite
    square well by an ideal gas. The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `mass`          The mass of the translating object
    `quantum`       ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    `active`        ``True`` if the mode is active in unimolecular processes, ``False`` if adiabatic
    =============== ============================================================
    
    Translational energies are much smaller than :math:`k_\\mathrm{B} T` except
    for temperatures approaching absolute zero, so a classical treatment of
    translation is more than adequate.
    """
    
    def __init__(self, mass=None, quantum=False, active=False):
        self.mass = mass
        self.quantum = quantum
        self.active = active

    property mass:
        """The mass of the translating object."""
        def __get__(self):
            return pq.Quantity(self._mass * constants.Na * 1000., pq.g / pq.mol)
        def __set__(self, value):
            if value is None or value == 0:
                self._mass = 0.0 
            else:
                self._mass = float(units.convertMass(value,"kg"))

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        IdealGasTranslation object.
        """
        return 'IdealGasTranslation(mass=({0:g},"{1}"), quantum={2!r}, active={3!r})'.format(float(self.mass), str(self.mass.dimensionality), self.quantum, self.active)

    def __reduce__(self):
        """
        A helper function used when pickling an IdealGasTranslation object.
        """
        return (IdealGasTranslation, (self.mass, self.quantum, self.active))
    
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function at the specified temperature
        `T` in K.
        """
        cdef double Q, qt 
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for IdealGasTranslation.')
        else:
            qt = ((2 * constants.pi * self._mass) / (constants.h * constants.h))**1.5 / 101325.
            Q = qt * (constants.kB * T)**2.5
        return Q
        
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the contribution to the heat capacity in J/(mol*K) due to
        translation at the specified temperature `T` in K.
        """
        cdef double Cv
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for IdealGasTranslation.')
        else:
            Cv = 1.5 * constants.R
        return Cv

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the contribution to the enthalpy in kJ/mol due to translation
        at the specified temperature `T` in K.
        """
        cdef double H
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for IdealGasTranslation.')
        else:
            H = 1.5 * 0.001 * constants.R * T
        return H

    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the contribution to the entropy in J/(mol*K) due to translation
        at the specified temperature `T` in K.
        """
        cdef double S
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for IdealGasTranslation.')
        else:
            S = (log(self.getPartitionFunction(T)) + 2.5) * constants.R
        return S

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in kJ/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double qt
        cdef numpy.ndarray sumStates
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for IdealGasTranslation.')
        elif sumStates0 is not None:
            sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist))
        else:
            Elist = Elist * 1000. / constants.Na
            qt = ((2 * constants.pi * self._mass) / (constants.h * constants.h))**1.5 / 101325.
            sumStates = qt * Elist**2.5 / (sqrt(constants.pi) * 15.0/8.0)
        return sumStates
            
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in kJ/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double qt, dE
        cdef numpy.ndarray densStates
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for IdealGasTranslation.')
        else:
            Elist = Elist * 1000. / constants.Na
            dE = Elist[1] - Elist[0]
            qt = ((2 * constants.pi * self._mass) / (constants.h * constants.h))**1.5 / 101325.
            densStates = qt * Elist**1.5 / (sqrt(constants.pi) * 0.75) * dE
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates
