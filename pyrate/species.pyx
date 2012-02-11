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
This module contains classes and functions for working with chemical species.
"""

import quantities as pq

cimport pyrate.constants as constants
import pyrate.units as units

################################################################################

cdef class LennardJones:
    """
    A set of Lennard-Jones collision parameters. The attributes are:

    =============== =================== ========================================
    Attribute       Type                Description
    =============== =================== ========================================
    `sigma`         :class:`Quantity`   Distance at which the inter-particle potential is zero
    `epsilon`       :class:`Quantity`   Depth of the potential well
    =============== =================== ========================================
    
    """

    def __init__(self, sigma=None, epsilon=None):
        self.sigma = sigma
        self.epsilon = epsilon

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        sigma = '({0:g},"{1}")'.format(float(self.sigma), str(self.sigma.dimensionality))
        epsilon = '({0:g},"{1}")'.format(float(self.epsilon), str(self.epsilon.dimensionality))
        return 'LennardJones(sigma={0}, epsilon={1})'.format(sigma, epsilon)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (LennardJones, (self.sigma, self.epsilon))
    
    property sigma:
        """The distance at which the inter-particle potential is zero."""
        def __get__(self):
            return pq.Quantity(self._sigma * 1e10, pq.angstrom)
        def __set__(self, value):
            if value is None or value == 0:
                self._sigma = 0.0 
            else:
                self._sigma = float(units.convertLength(value, "m"))

    property epsilon:
        """The depth of the potential well."""
        def __get__(self):
            return pq.Quantity(self._epsilon / constants.kB, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._epsilon = 0.0
            else:
                if isinstance(value, tuple):
                    value = pq.Quantity(value[0], value[1])
                if isinstance(value, pq.Quantity) and value.units == pq.K:
                    self._epsilon = float(value) * constants.kB
                else:
                    self._epsilon = float(units.convertEnergy(value, "J"))
