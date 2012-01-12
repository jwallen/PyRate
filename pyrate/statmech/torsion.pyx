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
This module contains classes that represent various models of torsional
motion. Torsional modes have both vibrational and rotational character, and
therefore should be treated semiclassically or quantum mechanically.
"""

import math
import numpy
import quantities as pq
import cython
from scipy.misc import factorial
from scipy.special import i0, i1, ellipk, ellipe
from libc.math cimport log, exp, sqrt, sin, cos

cimport pyrate.constants as constants
cimport pyrate.statmech.schrodinger as schrodinger 
import pyrate.units as units

################################################################################

cdef class HinderedRotor:
    """
    A statistical mechanical model of a one-dimensional hindered rotor.
    The attributes are:
    
    ======================== ===================================================
    Attribute                Description
    ======================== ===================================================
    `inertia`                The moment of inertia of the rotor
    `rotationalConstant`     The rotational constant of the rotor
    `symmetry`               The symmetry number of the rotor
    `barrier`                The barrier height of the rotor, if using the cosine potential
    `fourier`                An 2 x N array of Fourier series coefficients in J, if using the Fourier series potential
    `quantum`                ``True`` to use the quantum mechanical model, ``False`` to use the classical model
    `active`                 ``True`` if the mode is active in unimolecular processes, ``False`` if adiabatic
    ======================== ===================================================

    Note that the moment of inertia and the rotational constant are simply two
    ways of representing the same quantity; only one of these can be specified
    independently.
    """
    
    def __init__(self, inertia=None, symmetry=1, barrier=None, fourier=None, rotationalConstant=None, quantum=True, active=True):
        if inertia is not None and rotationalConstant is not None:
            raise ValueError('Only one of moment of inertia and rotational constant can be specified.')
        elif rotationalConstant is not None:
            self.rotationalConstant = rotationalConstant
        else:
            self.inertia = inertia
        self.symmetry = symmetry
        self.barrier = barrier
        self.quantum = quantum
        self.active = active
        self.fourier = fourier
 
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        HinderedRotor object.
        """
        inertia = '({0:g},"{1}")'.format(float(self.inertia), str(self.inertia.dimensionality))
        barrier = '({0:g},"{1}")'.format(float(self.barrier), str(self.barrier.dimensionality))
        if self._fourier is not None:
            fourier = self.fourier
            fourier = '([[{0}],[{1}]],"{2}")'.format(','.join(['{0:g}'.format(float(A)) for A in fourier[0,:]]), ','.join(['{0:g}'.format(float(B)) for B in fourier[1,:]]), str(fourier.dimensionality))
            return 'HinderedRotor(inertia={0}, symmetry={1:d}, barrier={2}, fourier={3}, quantum={4!r}, active={5!r})'.format(inertia, self.symmetry, barrier, fourier, self.quantum, self.active)
        else:
            return 'HinderedRotor(inertia={0}, symmetry={1:d}, barrier={2}, quantum={3!r}, active={4!r})'.format(inertia, self.symmetry, barrier, self.quantum, self.active)
            
    def __reduce__(self):
        """
        A helper function used when pickling a HinderedRotor object.
        """
        return (HinderedRotor, (self.inertia, self.symmetry, self.barrier, self.fourier, None, self.quantum, self.active))

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

    property barrier:
        """The barrier height of the cosine potential."""
        def __get__(self):
            return pq.Quantity(self._barrier * constants.Na * 0.001, "kJ/mol")
        def __set__(self, value):
            if value is None or value == 0:
                self._barrier = 0.0 
            else:
                self._barrier = float(units.convertEnergy(value, "J"))

    property fourier:
        """An 2 x N array of Fourier series coefficients in J, if using the Fourier series potential."""
        def __get__(self):
            if self._fourier is None:
                return None
            else:
                return pq.Quantity(self._fourier * constants.Na * 0.001, "kJ/mol")
        def __set__(self, value):
            if value is None:
                self._fourier = None
            else:
                self._fourier = numpy.array(units.convertEnergy(value, "J"))

    cpdef double getFrequency(self) except -1:
        """
        Return the frequency of vibration in cm^-1 corresponding to the limit of
        harmonic oscillation.
        """
        cdef double V0
        if self._fourier is not None:
            V0 = -numpy.sum(self._fourier[0,:])
        else:
            V0 = self._barrier
        return self.symmetry / 2.0 / constants.pi * sqrt(V0 / 2 / self._inertia) / (constants.c * 100)

    cpdef numpy.ndarray solveSchrodingerEquation(self, int Nbasis=401):
        """
        Solves the one-dimensional time-independent Schrodinger equation to 
        determine the energy levels of a one-dimensional hindered rotor with a
        Fourier series potential using `Nbasis` basis functions. For the
        purposes of this function it is usually sufficient to use 401 basis
        functions (the default). Returns the energy eigenvalues of the
        Hamiltonian matrix in J.
        """
        cdef int M, m, row, n
        cdef numpy.ndarray H
        cdef numpy.ndarray[numpy.float64_t,ndim=2] fourier
        cdef double A
    
        # The number of terms to use is 2*M + 1, ranging from -m to m inclusive
        if Nbasis % 2 == 0:
            M = Nbasis / 2
        else:
            M = (Nbasis - 1) / 2
        
        if self._fourier is not None:
            fourier = self._fourier / 2.0
            A = numpy.sum(self._fourier[0,:])
        else:
            assert self.symmetry >= 1 and self.symmetry <= 5
            fourier = numpy.zeros((2,5), numpy.float64)
            fourier[0,self.symmetry-1] = self._barrier / 2.0
            A = self._barrier / 2.0
        
        # Populate Hamiltonian matrix
        H = numpy.zeros((2*M+1,2*M+1), numpy.complex64)
        row = 0
        for m in range(-M, M+1):
            H[row,row] = A + constants.h * constants.h * m * m / (8 * constants.pi * constants.pi * self._inertia)
            for n in range(fourier.shape[0]):
                if row-n-1 > -1:    H[row,row-n-1] = complex(fourier[0,n], - fourier[1,n])
                if row+n+1 < 2*M+1: H[row,row+n+1] = complex(fourier[0,n], fourier[1,n])
            row += 1
        # The overlap matrix is the identity matrix, i.e. this is a standard
        # eigenvalue problem
        
        # Find the eigenvalues and eigenvectors of the Hamiltonian matrix
        E, V = numpy.linalg.eigh(H)
        self.energies = E - numpy.min(E)
        
        # Return the eigenvalues
        return self.energies

    cpdef double getPotential(self, double phi):
        """
        Return the values of the hindered rotor potential :math:`V(\\phi)`
        in kJ/mol at the angles `phi` in radians.
        """
        cdef double V = 0.0
        cdef int k
        if self._fourier is not None:
            for k in range(self._fourier.shape[1]):
                V += self._fourier[0,k] * cos((k+1) * phi) + self._fourier[1,k] * sin((k+1) * phi)
            V -= numpy.sum(self._fourier[0,:])
        else:
            V = 0.5 * self._barrier * (1 - cos(self.symmetry * phi))
        return V * 0.001 * constants.Na

    cdef double getRotationalConstantEnergy(self):
        """
        Return the value of the rotational constant in J.
        """
        return constants.hbar * constants.hbar / (2 * self._inertia)
    
    cpdef double getPartitionFunction(self, double T, bint semiclassical=False) except -1:
        """
        Return the value of the partition function at the specified temperature
        `T` in K. For the cosine potential, you can optionally apply the
        semi-classical correction factor by setting the `semiclassical`
        parameter to ``True`` (default is ``False``).
        """
        cdef double frequency, x, z, Q
        cdef double beta = 1.0 / (constants.R * T), V, dphi
        cdef int k
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            return numpy.sum(numpy.exp(-self.energies / constants.kB / T)) / self.symmetry
        elif self._fourier is not None:
            # Fourier series data found, so use it
            # Numerically evaluate the configuration integral
            Q = 0.0
            dphi = constants.pi/32.
            for phi in numpy.arange(0, 2*constants.pi, dphi):
                Q += exp(-beta * self.getPotential(phi) * 1000.) * dphi
            Q *= sqrt(constants.kB * T / (4 * constants.pi * self.getRotationalConstantEnergy())) / self.symmetry
        else:
            # No Fourier data, so use the cosine potential data
            frequency = self.getFrequency() * constants.c * 100
            x = constants.h * frequency / (constants.kB * T)
            z = 0.5 * self._barrier / (constants.kB * T)
            Q = sqrt(constants.pi * constants.kB * T / self.getRotationalConstantEnergy()) / self.symmetry * exp(-z) * i0(z)
        
        if semiclassical:
            Q *= x / (1 - exp(-x))
        
        return Q 
  
    cpdef double getHeatCapacity(self, double T) except -100000000:
        """
        Return the contribution to the heat capacity in J/(mol*K) due to rigid
        rotation at the specified temperature `T` in K.
        """
        cdef double frequency, x, z, exp_x, one_minus_exp_x, BB
        cdef double Tlow, Thigh, logQlow, logQhigh, logQ
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            E = self.energies
            e_kT = numpy.exp(-E / constants.kB / T)
            return (numpy.sum(E*E*e_kT) * numpy.sum(e_kT) - numpy.sum(E*e_kT)**2) / (constants.kB*T*T * numpy.sum(e_kT)**2) * constants.Na
        elif self._fourier is not None:
            # Fourier series data found, so use it
            Tlow = T * 0.999
            Thigh = T * 1.001
            logQlow = log(self.getPartitionFunction(Tlow))
            logQhigh = log(self.getPartitionFunction(Thigh))
            logQ = log(self.getPartitionFunction(T))
            return constants.R * T * T * (logQhigh - 2 * logQ + logQlow) / ((Thigh - T) * (T - Tlow)) + 2 * constants.R * T * (logQhigh - logQlow) / (Thigh - Tlow)
        else:
            # No Fourier data, so use the cosine potential data
            frequency = self.getFrequency() * constants.c * 100
            x = constants.h * frequency / (constants.kB * T)
            z = 0.5 * self._barrier / (constants.kB * T)
            exp_x = exp(x)
            one_minus_exp_x = 1.0 - exp_x
            BB = i1(z) / i0(z)
            return (x * x * exp_x / one_minus_exp_x / one_minus_exp_x - 0.5 + z * (z - BB - z * BB * BB)) * constants.R

    cpdef double getEnthalpy(self, double T) except 100000000:
        """
        Return the contribution to the enthalpy in kJ/mol due to rigid rotation
        at the specified temperature `T` in K.
        """
        cdef double Tlow, Thigh
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            E = self.energies
            e_kT = numpy.exp(-E / constants.kB / T)
            return numpy.sum(E*e_kT) / numpy.sum(e_kT) * constants.Na * 0.001
        else:
            # No Fourier data, so use the cosine potential data
            Tlow = T * 0.999
            Thigh = T * 1.001
            return (T *
                (log(self.getPartitionFunction(Thigh)) -
                log(self.getPartitionFunction(Tlow))) /
                (Thigh - Tlow)) * 0.001 * constants.R * T

    cpdef double getEntropy(self, double T) except -100000000:
        """
        Return the contribution to the entropy in J/(mol*K) due to rigid 
        rotation at the specified temperature `T` in K.
        """
        cdef double Tlow, Thigh
        
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            E = self.energies
            e_kT = numpy.exp(-E / constants.R / T)
            return constants.R * numpy.log(self.getPartitionFunction(T)) + numpy.sum(E*e_kT) / (T * numpy.sum(e_kT))
        else:
            # No Fourier data, so use the cosine potential data
            Tlow = T * 0.999
            Thigh = T * 1.001
            return (log(self.getPartitionFunction(Thigh)) +
                T * (log(self.getPartitionFunction(Thigh)) -
                log(self.getPartitionFunction(Tlow))) /
                (Thigh - Tlow)) * constants.R

    @cython.boundscheck(False)
    def getSumOfStates(self, numpy.ndarray Elist, numpy.ndarray sumStates0=None):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in kJ/mol above the ground state. If an initial sum of states 
        `sumStates0` is given, the rotor sum of states will be convoluted into
        these states.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] sumStates, _Elist = Elist
        cdef double q1f, pre, V0
        cdef int i
            
        if sumStates0 is not None:
            return schrodinger.convolve(sumStates0, self.getDensityOfStates(Elist))
        elif self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            energies = self.energies
            energy = lambda n: energies[n] if n < energies.shape[0] else 0
            return schrodinger.getSumOfStates(Elist * 1000. / constants.Na, energy, lambda n: 1, 0, sumStates0) / self.symmetry
        elif self._fourier is not None:
            # Fourier series data found, so use it
            raise NotImplementedError
        else:
            # No Fourier data, so use the cosine potential data
            _Elist = _Elist * 1000. / constants.Na
            sumStates = numpy.zeros_like(_Elist)
            q1f = sqrt(2 * constants.pi * self._inertia / constants.hbar / constants.hbar) / self.symmetry
            V0 = self._barrier
            pre = 4.0 * q1f * sqrt(V0 / (constants.pi * constants.pi * constants.pi))
            # The following is only valid in the classical limit
            # Note that ellipk(1) = infinity, so we must skip that value
            for i in range(len(_Elist)):
                if _Elist[i] < V0:
                    sumStates[i] = pre * (ellipe(_Elist[i] / V0) - (1 - _Elist[i] / V0) * ellipk(_Elist[i] / V0))
                elif _Elist[i] > V0:
                    sumStates[i] = pre * sqrt(_Elist[i] / V0) * ellipe(V0 / _Elist[i])
            return sumStates
            
    def getDensityOfStates(self, numpy.ndarray Elist, numpy.ndarray densStates0=None):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` in kJ/mol above the ground state. If an initial density
        of states `densStates0` is given, the rotor density of states will be
        convoluted into these states.
        """
        cdef numpy.ndarray[numpy.float64_t,ndim=1] densStates, _Elist = Elist
        cdef double q1f, pre, V0
        cdef int i
        
        _Elist = _Elist * 1000. / constants.Na
        if self.quantum:
            if self.energies is None: self.solveSchrodingerEquation()
            energies = self.energies
            energy = lambda n: energies[n] if n < energies.shape[0] else 0
            return schrodinger.getDensityOfStates(_Elist, energy, lambda n: 1, 0, densStates0) / self.symmetry
        elif self._fourier is not None:
            # Fourier series data found, so use it
            raise NotImplementedError
        else:
            # No Fourier data, so use the cosine potential data
            densStates = numpy.zeros_like(Elist)
            q1f = sqrt(2 * constants.pi * self._inertia / constants.hbar / constants.hbar) / self.symmetry
            V0 = self._barrier
            pre = 2.0 * q1f / sqrt(constants.pi * constants.pi * constants.pi * V0)
            # The following is only valid in the classical limit
            # Note that ellipk(1) = infinity, so we must skip that value
            for i in range(len(_Elist)):
                if _Elist[i] < V0:
                    densStates[i] = pre * ellipk(_Elist[i] / V0)
                elif _Elist[i] > V0:
                    densStates[i] = pre * sqrt(V0 / _Elist[i]) * ellipk(V0 / _Elist[i])
            densStates *= _Elist[1] - _Elist[0]
            if densStates0 is not None:
                return schrodinger.convolve(densStates0, densStates)
            else:
                return densStates
