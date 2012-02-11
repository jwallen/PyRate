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
This module contains classes and functions for working with chemical reactions.
"""

import quantities as pq

cimport pyrate.constants as constants
import pyrate.units as units

################################################################################

cdef class TransitionState:
    """
    A chemical transition state, representing a first-order saddle point on a
    potential energy surface. The attributes are:

    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `label`                 An identifying string label for the transition state
    `statmech`              The statistical mechanics model for the transition state
    `frequency`             The negative frequency of the first-order saddle point
    ======================= ====================================================

    """

    def __init__(self, label='', statmech=None, frequency=None):
        self.label = label
        self.statmech = statmech
        self.frequency = frequency

    def __repr__(self):
        """
        Return a string representation of the species.
        """
        frequency = '({0:g},"{1}")'.format(float(self.frequency), str(self.frequency.dimensionality))
        return 'TransitionState(label={0!r}, statmech={1!r}, frequency={2})'.format(self.label, self.statmech, frequency)

    def __reduce__(self):
        """
        A helper function used when pickling a TransitionState object.
        """
        return (TransitionState, (self.label, self.statmech, self.frequency))

    property frequency:
        """The negative frequency of the first-order saddle point."""
        def __get__(self):
            return pq.Quantity(self._frequency, pq.wavenumber)
        def __set__(self, value):
            if value is None or value == 0:
                self._frequency = 0.0 
            else:
                self._frequency = float(units.convertFrequency(value, pq.wavenumber))

################################################################################

cdef class Reaction:
    """
    A chemical reaction. The attributes are:
    
    ======================= ====================================================
    Attribute               Description
    ======================= ====================================================
    `reactants`             The list of reactant species
    `products`              The list of product species
    `kinetics`              The kinetics model to use for the reaction
    `reversible`            ``True`` if the reaction is reversible, ``False`` if not
    `transitionState`       Information about the transition state
    ======================= ====================================================
    
    """
    
    def __init__(self, reactants=None, products=None, kinetics=None, reversible=True, transitionState=None):
        self.reactants = reactants or []
        self.products = products or []
        self.kinetics = kinetics
        self.reversible = reversible
        self.transitionState = transitionState

    def __str__(self):
        """
        Return a string representation of the reaction, in the form 'A + B <=> C + D'.
        """
        arrow = ' <=> '
        if not self.reversible: arrow = ' -> '
        return arrow.join([' + '.join([str(s) for s in self.reactants]), ' + '.join([str(s) for s in self.products])])

    def __repr__(self):
        """
        Return a string representation of the reaction.
        """
        return 'Reaction(reactants={0!r}, products={1!r}, kinetics={2!r}, reversible={3!r}, transitionState={4!r})'.format(self.reactants, self.products, self.kinetics, self.reversible, self.transitionState)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Reaction, (self.reactants, self.products, self.kinetics, self.reversible, self.transitionState))
