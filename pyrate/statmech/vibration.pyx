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
This module contains classes that represent various models of vibrational
motion. For most molecular systems, a quantum treatment of vibrational motion
is required since typical vibrational energies are of similar magnitude to
:math:`k_\\mathrm{B} T`.
"""

import math
import numpy
import quantities as pq
from scipy.misc import factorial
from libc.math cimport log, exp

cimport pyrate.constants as constants
cimport pyrate.statmech.schrodinger as schrodinger 
import pyrate.units as units

################################################################################

cdef class HarmonicOscillator:
    """
    A statistical mechanical model of a set of one-dimensional independent
    harmonic oscillators. The attributes are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `frequencies`   The vibrational frequencies of the oscillators
    `quantum`       ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    `active`        ``True`` if the mode is active in unimolecular processes, ``False`` if adiabatic
    =============== ============================================================
    
    In the majority of chemical applications, the energy levels of the
    harmonic oscillator are of similar magnitude to :math:`k_\\mathrm{B} T`,
    requiring a quantum mechanical treatment. Fortunately, the harmonic
    oscillator has an analytical quantum mechanical solution.
    """
    
    def __init__(self, frequencies=None, quantum=True, active=True):
        self.frequencies = frequencies
        self.quantum = quantum
        self.active = active

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        HarmonicOscillator object.
        """
        frequencies = '([{0}],"{1}")'.format(','.join(['{0:g}'.format(float(freq)) for freq in self.frequencies]), str(self.frequencies.dimensionality))
        return 'HarmonicOscillator(frequencies={0}, quantum={1!r}, active={2!r})'.format(frequencies, self.quantum, self.active)

    def __reduce__(self):
        """
        A helper function used when pickling a HarmonicOscillator object.
        """
        return (HarmonicOscillator, (self.frequencies, self.quantum, self.active))

    property frequencies:
        """The vibrational frequencies of the oscillators."""
        def __get__(self):
            return pq.Quantity(self._frequencies / (constants.c * 100.), pq.wavenumber)
        def __set__(self, value):
            if value is None:
                self._frequencies = numpy.array([])
            else:
                self._frequencies = numpy.array(units.convertFrequency(value, "Hz"))
    
    cpdef double getPartitionFunction(self, double T) except -1:
        """
        Return the value of the partition function at the specified temperature
        `T` in K.
        """
        cdef double Q = 1.0, beta = 1.0 / (constants.kB * T), freq
        if self.quantum:
            for freq in self._frequencies:
                Q *= 1.0 / (1 - exp(-beta * constants.h * freq))
        else:
            for freq in self._frequencies:
                Q *= 1.0 / (beta * constants.h * freq)
        return Q
    
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the contribution to the heat capacity in J/(mol*K) due to rigid
        rotation at the specified temperature `T` in K.
        """
        cdef double Cv = 0.0, beta = 1.0 / (constants.kB * T), freq, x, exp_x
        if self.quantum:
            for freq in self._frequencies:
                x = beta * constants.h * freq
                exp_x = exp(x)
                Cv += x * x * exp_x / (1 - exp_x) / (1 - exp_x)
        else:
            Cv = len(self._frequencies)
        return Cv * constants.R

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the contribution to the enthalpy in kJ/mol due to rigid rotation
        at the specified temperature `T` in K.
        """
        cdef double H = 0.0, beta = 1.0 / (constants.kB * T), freq, x
        if self.quantum:
            for freq in self._frequencies:
                x = beta * constants.h * freq
                H += x / (exp(x) - 1)
        else:
            H = len(self._frequencies)
        return H * 0.001 * constants.R * T

    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the contribution to the entropy in J/(mol*K) due to rigid 
        rotation at the specified temperature `T` in K.
        """
        cdef double S = 0.0, beta = 1.0 / (constants.kB * T), freq, x
        S = log(self.getPartitionFunction(T))
        if self.quantum:
            for freq in self._frequencies:
                x = beta * constants.h * freq
                S += x / (exp(x) - 1)
        else:
            S += len(self._frequencies)
        return S * constants.R

    def getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in kJ/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef double freq
        cdef int Nfreq
        cdef numpy.ndarray sumStates
        if self.quantum:
            sumStates = sumStates0
            for freq in self._frequencies:
                sumStates = schrodinger.getSumOfStates(Elist * 1000. / constants.Na, lambda n: n * constants.h * freq, lambda n: 1, 0, sumStates)
        elif sumStates0 is not None:
            sumStates = schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist))
        else:
            Nfreq = len(self._frequencies)
            sumStates = (Elist * 1000. / constants.Na)**Nfreq / factorial(Nfreq)
            for freq in self._frequencies:
                sumStates /= constants.h * freq
        return sumStates
                    
    def getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in kJ/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef double freq, dE
        cdef int Nfreq
        cdef numpy.ndarray densStates
        Elist = Elist * 1000. / constants.Na
        if self.quantum:
            densStates = densStates0
            for freq in self._frequencies:
                densStates = schrodinger.getDensityOfStates(Elist, lambda n: n * constants.h * freq, lambda n: 1, 0, densStates)
        else:
            Nfreq = len(self._frequencies)
            dE = Elist[1] - Elist[0]
            densStates = Elist**(Nfreq-1) / factorial(Nfreq-1) * dE
            for freq in self._frequencies:
                densStates /= constants.h * freq
            if densStates0 is not None:
                densStates = schrodinger.convolve(densStates0, densStates)
        return densStates
