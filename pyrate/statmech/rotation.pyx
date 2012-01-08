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
This module contains classes that represent various models of rotational
motion. For most molecular systems, a classical treatment of rotational motion
is sufficient since typical rotational energies are much smaller in magnitude
than :math:`k_\\mathrm{B} T`.
"""

import math
import numpy
import quantities as pq
from libc.math cimport log, sqrt

cimport pyrate.constants as constants
import pyrate.statmech.schrodinger as schrodinger 
cimport pyrate.statmech.schrodinger as schrodinger 
import pyrate.units as units

################################################################################

cdef class LinearRotor:
    """
    A statistical mechanical model of a two-dimensional (linear) rigid rotor.
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    `active`                 ``True`` if the mode is active in unimolecular processes, ``False`` if adiabatic
    ======================== ===================================================
    
    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.

    In the majority of chemical applications, the energies involved in the
    rigid rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. 
    """
    
    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None, active=False):
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.symmetry = symmetry
        self.quantum = quantum
        self.active = active

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        LinearRotor object.
        """
        return 'LinearRotor(inertia=({0:g},"{1}"), symmetry={2:d}, quantum={3!r}, active={4!r})'.format(float(self.inertia), str(self.inertia.dimensionality), self.symmetry, self.quantum, self.active)

    def __reduce__(self):
        """
        A helper function used when pickling a LinearRotor object.
        """
        return (LinearRotor, (self.inertia, self.symmetry, self.quantum, None, self.active))

    property inertia:
        """The moment of inertia of the rotor."""
        def __get__(self):
            return pq.Quantity(self._inertia * constants.Na * 1e23, pq.u * pq.angstrom * pq.angstrom)
        def __set__(self, value):
            if value is None or value == 0:
                self._inertia = 0.0 
            else:
                self._inertia = float(units.convertInertia(value, "kg*m^2"))
    
    property rotationalConstant:
        """The rotational constant of the rotor."""
        def __get__(self):
            B = constants.h / (8 * constants.pi * constants.pi * self._inertia) / (constants.c * 100.)
            return pq.Quantity(B, pq.wavenumber)
        def __set__(self, value):
            if value is None or value == 0:
                self._inertia = 0.0 
            else:
                value = float(units.convertFrequency(value, "Hz"))
                self._inertia = constants.h / (8 * constants.pi * constants.pi * value)

    cdef double getRotationalConstantEnergy(self):
        """
        Return the value of the rotational constant in J.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia)
            
    cpdef double getLevelEnergy(self, int J):
        """
        Return the energy of level `J` in J.
        """
        return self.getRotationalConstantEnergy() * J * (J + 1)
    
    cpdef int getLevelDegeneracy(self, int J):
        """
        Return the degeneracy of level `J`.
        """
        return 2 * J + 1
    
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function at the specified temperature
        `T` in K.
        """
        cdef double B, Q
        if self.quantum:
            Q = schrodinger.getPartitionFunction(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) / self.symmetry
        else:
            B = self.getRotationalConstantEnergy()
            Q = constants.kB * T / B / self.symmetry
        return Q
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the contribution to the heat capacity in J/(mol*K) due to rigid
        rotation at the specified temperature `T` in K.
        """
        cdef double Cv
        if self.quantum:
            Cv = schrodinger.getHeatCapacity(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) * constants.R
        else:
            Cv = constants.R
        return Cv

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the contribution to the enthalpy in kJ/mol due to rigid rotation
        at the specified temperature `T` in K.
        """
        cdef double H
        if self.quantum:
            H = schrodinger.getEnthalpy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) * 0.001 * constants.R * T
        else:
            H = 0.001 * constants.R * T
        return H
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the contribution to the entropy in J/(mol*K) due to rigid 
        rotation at the specified temperature `T` in K.
        """
        cdef double S
        if self.quantum:
            S = (schrodinger.getEntropy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) - log(self.symmetry)) * constants.R
        else:
            S = (numpy.log(self.getPartitionFunction(T)) + 1.0) * constants.R
        return S

    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in kJ/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double B
        cdef numpy.ndarray sumStates
        if self.quantum:
            sumStates = schrodinger.getSumOfStates(Elist * 1000. / constants.Na, self.getLevelEnergy, self.getLevelDegeneracy, 0, sumStates0) / self.symmetry
        elif sumStates0 is not None:
            sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist, densStates0=None))
        else:
            Elist = Elist * 1000. / constants.Na
            B = self.getRotationalConstantEnergy()
            sumStates = Elist / B / self.symmetry
        return sumStates
            
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in kJ/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double B, dE
        cdef numpy.ndarray densStates
        Elist = Elist * 1000. / constants.Na
        if self.quantum:
            densStates = schrodinger.getDensityOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, densStates0) / self.symmetry
        else:
            B = self.getRotationalConstantEnergy()
            dE = Elist[1] - Elist[0]
            densStates = numpy.ones_like(Elist) * dE / B / self.symmetry
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates

################################################################################

cdef class NonlinearRotor:
    """
    A statistical mechanical model of an N-dimensional nonlinear rigid rotor.
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moments of inertia of the rotor
    `rotationalConstant`     The rotational constants of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    `active`                 ``True`` if the mode is active in unimolecular processes, ``False`` if adiabatic
    ======================== ===================================================

    Note that the moments of inertia and the rotational constants are simply two
    ways of representing the same quantity; only one set of these can be 
    specified independently.
    
    In the majority of chemical applications, the energies involved in the
    rigid rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. In the
    current implementation, the quantum mechanical model has not been 
    implemented, and a :class:`NotImplementedError` will be raised if you try
    to use it. 
    """
    
    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None, active=False):
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.symmetry = symmetry
        self.quantum = quantum
        self.active = active
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        NonlinearRotor object.
        """
        inertia = '([{0}],"{1}")'.format(','.join(['{0:g}'.format(float(I)) for I in self.inertia]), str(self.inertia.dimensionality))
        return 'NonlinearRotor(inertia={0}, symmetry={1:d}, quantum={2!r}, active={3!r})'.format(inertia, self.symmetry, self.quantum, self.active)

    def __reduce__(self):
        """
        A helper function used when pickling a NonlinearRotor object.
        """
        return (NonlinearRotor, (self.inertia, self.symmetry, self.quantum, None, self.active))

    property inertia:
        """The moments of inertia of the rotor."""
        def __get__(self):
            return pq.Quantity(self._inertia * constants.Na * 1e23, pq.u * pq.angstrom * pq.angstrom)
        def __set__(self, value):
            if value is None:
                self._inertia = numpy.array([])
            else:
                self._inertia = numpy.array(units.convertInertia(value, "kg*m^2"))
    
    property rotationalConstant:
        """The rotational constants of the rotor."""
        def __get__(self):
            B = constants.h / (8 * constants.pi * constants.pi * self._inertia) / (constants.c * 100.)
            return pq.Quantity(B, pq.wavenumber)
        def __set__(self, value):
            if value is None:
                self._inertia = 0.0 
            else:
                value = numpy.array(units.convertFrequency(value, "Hz"))
                self._inertia = constants.h / (8 * constants.pi * constants.pi * value)

    cdef numpy.ndarray getRotationalConstantEnergy(self):
        """
        Return the values of the rotational constants in J.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia)

    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function at the specified temperature
        `T` in K.
        """
        cdef double Q, theta = 1.0
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            for B in self.getRotationalConstantEnergy():
                theta *= B / constants.kB
            Q = sqrt(constants.pi * T**self.inertia.shape[0] / theta) / self.symmetry
        return Q
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the contribution to the heat capacity in J/(mol*K) due to rigid
        rotation at the specified temperature `T` in K.
        """
        cdef double Cv
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            Cv = self._inertia.shape[0] * constants.R / 2.
        return Cv

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the contribution to the enthalpy in kJ/mol due to rigid rotation
        at the specified temperature `T` in K.
        """
        cdef double H
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            H = self._inertia.shape[0] * 0.001 * constants.R * T / 2.
        return H
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the contribution to the entropy in J/(mol*K) due to rigid 
        rotation at the specified temperature `T` in K.
        """
        cdef double S
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            S = (numpy.log(self.getPartitionFunction(T)) + self._inertia.shape[0] / 2.) * constants.R
        return S
    
    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in kJ/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double theta = 1.0
        cdef numpy.ndarray sumStates
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            assert len(self.inertia) == 3
            if sumStates0 is not None:
                sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist, densStates0=None))
            else:
                Elist = Elist * 1000. / constants.Na
                for B in self.getRotationalConstantEnergy():
                    theta *= B
                sumStates = 4.0/3.0 * Elist * numpy.sqrt(Elist / theta) / self.symmetry
        return sumStates
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in kJ/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double theta = 1.0, inertia, dE
        cdef numpy.ndarray densStates
        if self.quantum:
            raise NotImplementedError('Quantum mechanical model not yet implemented for NonlinearRotor.')
        else:
            assert len(self.inertia) == 3
            Elist = Elist * 1000. / constants.Na
            dE = Elist[1] - Elist[0]
            for B in self.getRotationalConstantEnergy():
                theta *= B
            densStates = 2.0 * numpy.sqrt(Elist / theta) / self.symmetry * dE
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates

################################################################################

cdef class KRotor:
    """
    A statistical mechanical model of an active K-rotor (a one-dimensional
    rigid rotor). The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    `symmetry`               The symmetry number of the rotor
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    `active`                 ``True`` if the mode is active in unimolecular processes, ``False`` if adiabatic
    ======================== ===================================================
    
    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.
    
    In the majority of chemical applications, the energies involved in the
    K-rotor place it very nearly in the classical limit at all relevant
    temperatures; therefore, the classical model is used by default. 
    """
    
    def __init__(self, inertia=None, symmetry=1, quantum=False, rotationalConstant=None, active=True):
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.symmetry = symmetry
        self.quantum = quantum
        self.active = active
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KRotor object.
        """
        return 'KRotor(inertia=({0:g},"{1}"), symmetry={2:d}, quantum={3!r}, active={4!r})'.format(float(self.inertia), str(self.inertia.dimensionality), self.symmetry, self.quantum, self.active)

    def __reduce__(self):
        """
        A helper function used when pickling a KRotor object.
        """
        return (KRotor, (self.inertia, self.symmetry, self.quantum, None, self.active))

    property inertia:
        """The moment of inertia of the rotor."""
        def __get__(self):
            return pq.Quantity(self._inertia * constants.Na * 1e23, pq.u * pq.angstrom * pq.angstrom)
        def __set__(self, value):
            if value is None or value == 0:
                self._inertia = 0.0 
            else:
                self._inertia = float(units.convertInertia(value, "kg*m^2"))
    
    property rotationalConstant:
        """The rotational constant of the rotor."""
        def __get__(self):
            B = constants.h / (8 * constants.pi * constants.pi * self._inertia) / (constants.c * 100.)
            return pq.Quantity(B, pq.wavenumber)
        def __set__(self, value):
            if value is None or value == 0:
                self._inertia = 0.0 
            else:
                value = float(units.convertFrequency(value, "Hz"))
                self._inertia = constants.h / (8 * constants.pi * constants.pi * value)

    cdef double getRotationalConstantEnergy(self):
        """
        Return the value of the rotational constant in J.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia)
            
    cpdef double getLevelEnergy(self, int J):
        """
        Return the energy of level `J` in J.
        """
        return self.getRotationalConstantEnergy() * J * J
    
    cpdef int getLevelDegeneracy(self, int J):
        """
        Return the degeneracy of level `J`.
        """
        return 1 if J == 0 else 2
        
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function at the specified temperature
        `T` in K.
        """
        cdef double Q, B
        if self.quantum:
            Q = schrodinger.getPartitionFunction(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) / self.symmetry
        else:
            B = self.getRotationalConstantEnergy()
            Q = sqrt(constants.pi * constants.kB * T / B) / self.symmetry
        return Q
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the contribution to the heat capacity in J/(mol*K) due to rigid
        rotation at the specified temperature `T` in K.
        """
        cdef double Cv
        if self.quantum:
            Cv = schrodinger.getHeatCapacity(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) * constants.R
        else:
            Cv = 0.5 * constants.R
        return Cv
    
    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the contribution to the enthalpy in kJ/mol due to rigid rotation
        at the specified temperature `T` in K.
        """
        cdef double H
        if self.quantum:
            H = schrodinger.getEnthalpy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) * 0.001 * constants.R * T
        else:
            H = 0.5 * 0.001 * constants.R * T
        return H
    
    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the contribution to the entropy in J/(mol*K) due to rigid 
        rotation at the specified temperature `T` in K.
        """
        cdef double S
        if self.quantum:
            S = (schrodinger.getEntropy(T, self.getLevelEnergy, self.getLevelDegeneracy, 0) - log(self.symmetry)) * constants.R
        else:
            S = (numpy.log(self.getPartitionFunction(T)) + 0.5) * constants.R
        return S
    
    cpdef numpy.ndarray getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in kJ/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double B
        cdef numpy.ndarray sumStates
        if self.quantum:
            sumStates = schrodinger.getSumOfStates(Elist * 1000. / constants.Na, self.getLevelEnergy, self.getLevelDegeneracy, 0, sumStates0) / self.symmetry
        elif sumStates0 is not None:
            sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist, densStates0=None))
        else:
            B = self.getRotationalConstantEnergy()
            Elist = Elist * 1000. / constants.Na
            sumStates = 2.0 * numpy.sqrt(Elist / B) / self.symmetry
        return sumStates
    
    cpdef numpy.ndarray getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in kJ/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double B, dE
        cdef int r
        cdef numpy.ndarray densStates
        Elist = Elist * 1000. / constants.Na
        if self.quantum:
            densStates = schrodinger.getDensityOfStates(Elist, self.getLevelEnergy, self.getLevelDegeneracy, 0, densStates0) / self.symmetry
        else:
            densStates = numpy.zeros(Elist.shape[0], dtype=numpy.float)
            B = self.getRotationalConstantEnergy()
            dE = Elist[1] - Elist[0]
            for r in range(Elist.shape[0]):
                if Elist[r] == 0: continue
                densStates[r] = 1.0 / sqrt(B * Elist[r]) * dE / self.symmetry
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates
