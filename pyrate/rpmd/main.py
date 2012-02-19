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
This module contains the main execution functionality for various ring polymer
molecular dynamics (RPMD) calculations.
"""

import math
import numpy
import logging
import quantities as pq

import pyrate.units as units
import pyrate.constants as constants

from ._main import *

################################################################################

class RPMDError(Exception):
    """
    An exception raised when an error occurs during an RPMD simulation. Pass a
    string describing the circumstances of the exceptional behavior.
    """
    pass

################################################################################

class RPMD:
    """
    A representation of a ring polymer molecular dynamics (RPMD) job for
    computing gas-phase chemical reaction rates. The attributes are:
    
    =========================== ================================================
    Attribute                   Description
    =========================== ================================================
    `labels`                    The labels to use for each atom in the molecular system
    `mass`                      The mass of each atom in the molecular system
    `Natoms`                    The number of atoms in the molecular system
    `reactants`                 The dividing surface near the reactants, as a :class:`Reactants` object
    `transitionState`           The dividing surface near the transition state, as a :class:`TransitionState` object
    `potential`                 A function that computes the potential and forces for a given position
    --------------------------- ------------------------------------------------
    `beta`                      The reciprocal temperature of the RPMD simulation
    `dt`                        The time step to use in the RPMD simulation
    `Nbeads`                    The number of beads per atom in the RPMD simulation
    `xi_current`                The current value of the reaction coordinate
    `mode`                      A flag indicating the type of RPMD calculation currently underway (1 = umbrella, 2 = recrossing)
    =========================== ================================================
    
    """

    def __init__(self, labels, reactants, transitionState, potential):
        """
        Initialize an RPMD object. The `mass` of each atom should be given in
        g/mol, while the `Rinf` value should be given in angstroms. (They will
        be converted to atomic units.)
        """
        self.mass = reactants.mass
        self.labels = labels
        self.Natoms = len(self.mass)
        self.reactants = reactants
        self.transitionState = transitionState
        self.potential = potential
        
        self.beta = 0
        self.dt = 0
        self.Nbeads = 0
        self.xi_current = 0
        self.mode = 0

    def computeTransmissionCoefficient(self, T, Nbeads, dt, 
                                       equilibrationTime,
                                       xi_current,
                                       parentEvolutionTime,
                                       childrenPerSampling,
                                       childEvolutionTime,
                                       childSamplingTime,
                                       saveParentTrajectory=False, 
                                       saveChildTrajectories=False):
        """
        Return the transmission coefficient by the recrossing method. In this
        approach, a constrained RPMD simulation is initiated in the presence of
        an Andersen thermostat to generate a series of independent 
        configurations with centroids on the transition state dividing surface
        :math:`\\bar{s}_1(\\mathbf{q}) = 0`. For each of these "parent"
        configurations, a set of "child" trajectories is spawned by sampling
        from a Maxwell-Boltzmann distribution, with each trajectory evolving
        in time without the dividing surface constraint or Andersen thermostat.
        The transmission coefficient is then computed via
        
        .. math::
        
            \\kappa^{(n)}(s_1) = \\lim_{t \\rightarrow \\infty} 
                \\left< \\bar{f}_{s_1}(\\mathbf{q})^{-1} \\bar{v}_{s_1}(\\mathbf{p}, \\mathbf{q}) h \\left[ \\bar{s}_1(\\mathbf{q}_t) \\right] \\right>
        
        where the bracketed quantity is averaged over a large number of the
        child trajectories spawned from a large number of parent configurations.
        
        The `saveParentTrajectory` and `saveChildTrajectories` flags enable
        saving of the parent trajectory and/or a sampling of child trajectories
        as XYZ data files for later visualization in programs such as VMD.
        This is off by default because it is very slow. 
        """
        
        # Set the parameters for the RPMD calculation
        self.beta = 4.35974417e-18 / (constants.kB * T)
        self.dt = dt / 2.418884326505e-5
        self.Nbeads = Nbeads
        equilibrationTime /= 2.418884326505e-5
        parentEvolutionTime /= 2.418884326505e-5
        childEvolutionTime /= 2.418884326505e-5
        childSamplingTime /= 2.418884326505e-5
        geometry = self.transitionState.geometry[:,:,0]
        self.xi_current = xi_current
        self.mode = 2
        
        logging.info('*****************************')
        logging.info('RPMD transmission coefficient')
        logging.info('*****************************')
        logging.info('')
        
        equilibrationSteps = int(round(equilibrationTime / self.dt))
        parentEvolutionSteps = int(round(parentEvolutionTime / self.dt))
        childEvolutionSteps = int(round(childEvolutionTime / self.dt))
        childSamplingSteps = int(round(childSamplingTime / self.dt))
        
        logging.info('Parameters')
        logging.info('==========')
        logging.info('Temperature                             = {0:g} K'.format(T))
        logging.info('Number of beads                         = {0:d}'.format(Nbeads))
        logging.info('Reaction coordinate                     = {0:g}'.format(xi_current))
        logging.info('Time step                               = {0:g} ps'.format(self.dt * 2.418884326505e-5))
        logging.info('Length of parent trajectory             = {0:g} ps ({1:d} steps)'.format(parentEvolutionSteps * self.dt * 2.418884326505e-5, parentEvolutionSteps))
        logging.info('Initial parent equilibration time       = {0:g} ps ({1:d} steps)'.format(equilibrationSteps * self.dt * 2.418884326505e-5, equilibrationSteps))
        logging.info('Frequency of child trajectory sampling  = {0:g} ps ({1:d} steps)'.format(childSamplingSteps * self.dt * 2.418884326505e-5, childSamplingSteps))
        logging.info('Length of child trajectories            = {0:g} ps ({1:d} steps)'.format(childEvolutionSteps * self.dt * 2.418884326505e-5, childEvolutionSteps))
        logging.info('Number of children per sampling         = {0:d}'.format(childrenPerSampling))
        logging.info('')
        
        kappa = compute_transmission_coefficient(
            self.beta, Nbeads, self.dt, xi_current,
            geometry, self.mass, self.potential,
            equilibrationSteps, parentEvolutionSteps, childSamplingSteps,
            childrenPerSampling, childEvolutionSteps, saveParentTrajectory,
            saveChildTrajectories,
            self.reactants.Rinf, self.reactants.massFractions, self.reactants.reactant1Atoms, self.reactants.reactant2Atoms,
            self.transitionState.formingBonds, self.transitionState.formingBondLengths,
            self.transitionState.breakingBonds, self.transitionState.breakingBondLengths,
        )
        
        return kappa
    
    def evolve(self, p, q, V, dVdq, xi, dxi, d2xi, constrain=False):
        xi = numpy.array(xi, order='F')
        evolve(p, q, V, dVdq, xi, dxi, d2xi, constrain,
            self.beta, self.dt, self.reactants.mass, self.potential,
            self.reactants.Rinf, self.reactants.massFractions, self.reactants.reactant1Atoms, self.reactants.reactant2Atoms,
            self.transitionState.formingBonds, self.transitionState.formingBondLengths,
            self.transitionState.breakingBonds, self.transitionState.breakingBondLengths,
            self.xi_current, self.mode)
        return xi
    
    def getReactionCoordinate(self, centroid):
        """
        Return the value and gradient of the reaction coordinate at the
        given `centroid`.
        """
        return get_reaction_coordinate(centroid, 
            self.reactants.Rinf, self.reactants.massFractions, self.reactants.reactant1Atoms, self.reactants.reactant2Atoms,
            self.transitionState.formingBonds, self.transitionState.formingBondLengths,
            self.transitionState.breakingBonds, self.transitionState.breakingBondLengths,
            self.xi_current, self.mode)
        
    def getPotential(self, q, xi, dxi, d2xi):
        """
        Return the potential and the force of the system of ring polymers at
        the given position `q`.
        """
        V, dVdq = get_potential(q, xi, dxi, d2xi)
        return V, dVdq

    def getCentroid(self, position):
        """
        Return the centroid of each ring polymer in the system of ring polymers
        at a given `position`.
        """
        centroid = numpy.zeros((3,self.Natoms))
        for i in range(3):
            for j in range(self.Natoms):
                for k in range(self.Nbeads):
                    centroid[i,j] += position[i,j,k]
                centroid[i,j] /= self.Nbeads
        
        return centroid
