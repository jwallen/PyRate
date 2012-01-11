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
This module provides the :class:`StatMech` class, used for working with
statistical mechanical models of a molecular system involving multiple degrees
of freedom.
"""

import math
import numpy

import pyrate.constants as constants
import pyrate.units as units

################################################################################

class StatMech:
    """
    A class for working with statistical mechanical models of a molecular
    system. The attributes are:
    
    =================== ========================================================
    Attribute           Description
    =================== ========================================================
    `modes`             A list of the molecular degrees of freedom
    `E0`                The relative ground-state electronic energy, including zero-point energies
    `spinMultiplicity`  The degeneracy of the electronic ground state
    `opticalIsomers`    The number of optical isomers
    =================== ========================================================
    
    Note that `E0` and `spinMultiplicity` reflect the electronic mode of the
    molecular system.    
    """
    
    def __init__(self, modes=None, E0=None, spinMultiplicity=1, opticalIsomers=1):
        self.modes = modes or []
        if E0 is not None:
            self.E0 = float(units.convertEnergy(E0, "kJ/mol"))
        else:
            self.E0 = 0.0
        self.spinMultiplicity = spinMultiplicity
        self.opticalIsomers = opticalIsomers

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        StatMech object.
        """
        E0 = '({0:g},"{1}")'.format(float(self.E0), str(self.E0.dimensionality))
        return 'StatMech(modes={0!r}, E0={1}, spinMultiplicity={2:d}, opticalIsomers={3:d})'.format(self.modes, E0, self.spinMultiplicity, self.opticalIsomers)

    @property
    def E0(self):
        """The relative ground-state electronic energy, including zero-point energies."""
        return pq.Quantity(self._E0 * constants.Na * 0.001, "kJ/mol")    
    @E0.setter
    def E0(self, value):
        if value is None or value == 0:
            self._E0 = 0.0 
        else:
            self._E0 = float(units.convertEnergy(value, "J"))

    def getPartitionFunction(self, T, active=True, adiabatic=True):
        """
        Return the total value of the partition function of the system at the
        specified temperature `T` in K. The `active` and `adiabatic` flags 
        control whether the returned density of states includes the active 
        and/or adiabatic modes.
        """
        Q = 1.0
        for mode in self.modes:
            if (mode.active and active) or (not mode.active and adiabatic):
                Q *= mode.getPartitionFunction(T)
        return Q * self.spinMultiplicity * self.opticalIsomers

    def getHeatCapacity(self, T):
        """
        Return the total heat capacity in J/(mol*K) of the system at the 
        specified temperature `T` in K.
        """
        Cp = constants.R
        for mode in self.modes:
            Cp += mode.getHeatCapacity(T)
        return Cp

    def getEnthalpy(self, T):
        """
        Return the total enthalpy in J/(mol*K) of the system at the specified 
        temperature `T` in K.
        """
        H = 0.001 * constants.R * T
        for mode in self.modes:
            H += mode.getEnthalpy(T)
        return H
    
    def getEntropy(self, T):
        """
        Return the total entropy in J/(mol*K) of the system at the specified 
        temperature `T` in K.
        """
        S = constants.R
        for mode in self.modes:
            S += mode.getEntropy(T)
        return S

    def getSumOfStates(self, Elist, active=True, adiabatic=True, checkEnergySpacing=True):
        """
        Return the sum of states :math:`N(E)` at the specified energies `Elist`
        in kJ/mol above the ground state. The `active` and `adiabatic` flags 
        control whether the returned density of states includes the active 
        and/or adiabatic modes. (By default all modes are included.) The
        `checkEnergySpacing` flag is used to turn off a check of the provided
        energy spacing, which normally warns the user if the energy spacing is
        too large.
        """
        sumStates = None
        
        # Warn if using a large energy spacing and convolving modes using the
        # Beyer-Swinehart algorithm, which requires a small energy spacing to
        # be accurate
        if checkEnergySpacing and any([mode.quantum and mode.active for mode in self.modes]):
            dE = Elist[1] - Elist[0]
            if dE > 0.1:
                print('RuntimeWarning: Energy spacing of {0:g} kJ/mol detected; an energy spacing smaller than 0.1 kJ/mol is strongly recommended when computing the sum of states for quantum modes.'.format(float(dE)))
         
        # Classical modes added first since they give a smooth sum of states
        # Adding quantum modes after classical modes still results in a smooth sum of states
        # Only include active and/or adiabatic modes if requested
        if active and adiabatic: 
            quantum_active = [(False, True), (False, False), (True, True), (True, False)]
        elif active: 
            quantum_active = [(False, True), (True, True)]
        elif adiabatic: 
            quantum_active = [(False, False), (True, False)]
        else: 
            quantum_active = []
        for q, a in quantum_active:
            for mode in self.modes:
                if mode.quantum == q and mode.active == a:
                    sumStates = mode.getSumOfStates(Elist, sumStates0=sumStates)
                    
        return sumStates * self.spinMultiplicity * self.opticalIsomers
        
    def getDensityOfStates(self, Elist, active=True, adiabatic=True, checkEnergySpacing=True):
        """
        Return the density of states :math:`\\rho(E) \\ dE` at the specified
        energies `Elist` above the ground state. The `active` and `adiabatic`
        flags control whether the returned density of states includes the
        active and/or adiabatic modes. (By default all modes are included.) The
        `checkEnergySpacing` flag is used to turn off a check of the provided
        energy spacing, which normally warns the user if the energy spacing is
        too large.
        """
        densStates = None
        
        # Warn if using a large energy spacing and convolving modes using the
        # Beyer-Swinehart algorithm, which requires a small energy spacing to
        # be accurate
        if checkEnergySpacing and any([mode.quantum and mode.active for mode in self.modes]):
            dE = Elist[1] - Elist[0]
            if dE > 0.1:
                print('RuntimeWarning: Energy spacing of {0:g} kJ/mol detected; an energy spacing smaller than 0.1 kJ/mol is strongly recommended when computing the density of states for quantum modes.'.format(float(dE)))
         
        # Classical modes added first since they give a smooth density of states
        # Adding quantum modes after classical modes still results in a smooth density of states
        # Only include active and/or adiabatic modes if requested
        if active and adiabatic: 
            quantum_active = [(False, True), (False, False), (True, True), (True, False)]
        elif active: 
            quantum_active = [(False, True), (True, True)]
        elif adiabatic: 
            quantum_active = [(False, False), (True, False)]
        else: 
            quantum_active = []
        for q, a in quantum_active:
            for mode in self.modes:
                if mode.quantum == q and mode.active == a:
                    densStates = mode.getDensityOfStates(Elist, densStates0=densStates)
                    
        return densStates * self.spinMultiplicity * self.opticalIsomers
