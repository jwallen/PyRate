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
This script contains unit tests of the :mod:`pyrate.reaction` module.
"""

import unittest
import math

from pyrate.reaction import *
from pyrate.species import Species
from pyrate.kinetics import Arrhenius

################################################################################

class TestTransitionState(unittest.TestCase):
    """
    Contains unit tests of the :class:`TransitionState` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.label = 'TS1'
        self.statmech = None
        self.frequency = pq.Quantity(1000, "cm^-1")
        self.transitionState = TransitionState(
            label = self.label, 
            statmech = self.statmech, 
            frequency = self.frequency,
        )
        
    def test_label(self):
        """
        Test that the TransitionState label property was properly set.
        """
        self.assertEqual(self.transitionState.label, self.label)

    def test_frequency(self):
        """
        Test that the TransitionState frequency property was properly set.
        """
        self.assertAlmostEqual(self.transitionState.frequency / self.frequency, 1.0, 6, '{0} != {1} within 6 places'.format(self.transitionState.frequency, self.frequency))

    def test_pickle(self):
        """
        Test that a TransitionState object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        transitionState = cPickle.loads(cPickle.dumps(self.transitionState))
        self.assertEqual(self.transitionState.label, transitionState.label)
        self.assertEqual(self.transitionState.statmech, transitionState.statmech)
        self.assertAlmostEqual(self.transitionState.frequency / transitionState.frequency, 1.0, 6, '{0} != {1} within 6 places'.format(self.transitionState.frequency, transitionState.frequency))
   
    def test_repr(self):
        """
        Test that a TransitionState object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('transitionState = {0!r}'.format(self.transitionState))
        self.assertEqual(self.transitionState.label, transitionState.label)
        self.assertEqual(self.transitionState.statmech, transitionState.statmech)
        self.assertAlmostEqual(self.transitionState.frequency / transitionState.frequency, 1.0, 6, '{0} != {1} within 6 places'.format(self.transitionState.frequency, transitionState.frequency))

################################################################################

class TestReaction(unittest.TestCase):
    """
    Contains unit tests of the :class:`Reaction` class.
    """

    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.reactants = [Species(label='C2H4'), Species(label='H')]
        self.products = [Species(label='C2H5')]
        self.reversible = True
        self.kinetics = Arrhenius(A=pq.Quantity(1.0e12,"cm^3/(mol*s)"), n=0.5, Ea=pq.Quantity(10.0,"kcal/mol"), T0=1*pq.K)
        self.transitionState = None
        self.reaction = Reaction(
            reactants = self.reactants,
            products = self.products,
            reversible = self.reversible,
            kinetics = self.kinetics,
            transitionState = self.transitionState,
        )
        
    def test_reactants(self):
        """
        Test that the Reaction reactants property was properly set.
        """
        self.assertEqual(self.reaction.reactants, self.reactants)

    def test_reactants(self):
        """
        Test that the Reaction products property was properly set.
        """
        self.assertEqual(self.reaction.products, self.products)

    def test_reversible(self):
        """
        Test that the Reaction reversible property was properly set.
        """
        self.assertEqual(self.reaction.reversible, self.reversible)

    def test_kinetics(self):
        """
        Test that the Reaction kinetics property was properly set.
        """
        self.assertAlmostEqual(self.reaction.kinetics.A / self.kinetics.A, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.A, self.kinetics.A))
        self.assertAlmostEqual(self.reaction.kinetics.n / self.kinetics.n, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.n, self.kinetics.n))
        self.assertAlmostEqual(self.reaction.kinetics.Ea / self.kinetics.Ea, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.Ea, self.kinetics.Ea))
        self.assertAlmostEqual(self.reaction.kinetics.T0 / self.kinetics.T0, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.T0, self.kinetics.T0))

    def test_pickle(self):
        """
        Test that a TransitionState object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        reaction = cPickle.loads(cPickle.dumps(self.reaction))
        self.assertEqual(self.reaction.reactants[0].label, reaction.reactants[0].label)
        self.assertEqual(self.reaction.reactants[1].label, reaction.reactants[1].label)
        self.assertEqual(self.reaction.products[0].label, reaction.products[0].label)
        self.assertEqual(self.reaction.reversible, reaction.reversible)
        self.assertAlmostEqual(self.reaction.kinetics.A / reaction.kinetics.A, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.A, reaction.kinetics.A))
        self.assertAlmostEqual(self.reaction.kinetics.n / reaction.kinetics.n, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.n, reaction.kinetics.n))
        self.assertAlmostEqual(self.reaction.kinetics.Ea / reaction.kinetics.Ea, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.Ea, reaction.kinetics.Ea))
        self.assertAlmostEqual(self.reaction.kinetics.T0 / reaction.kinetics.T0, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.T0, reaction.kinetics.T0))
   
    def test_repr(self):
        """
        Test that a TransitionState object can be successfully reconstructed
        from its repr() output with no loss of information.
        """
        exec('reaction = {0!r}'.format(self.reaction))
        self.assertEqual(self.reaction.reactants[0].label, reaction.reactants[0].label)
        self.assertEqual(self.reaction.reactants[1].label, reaction.reactants[1].label)
        self.assertEqual(self.reaction.products[0].label, reaction.products[0].label)
        self.assertEqual(self.reaction.reversible, reaction.reversible)
        self.assertAlmostEqual(self.reaction.kinetics.A / reaction.kinetics.A, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.A, reaction.kinetics.A))
        self.assertAlmostEqual(self.reaction.kinetics.n / reaction.kinetics.n, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.n, reaction.kinetics.n))
        self.assertAlmostEqual(self.reaction.kinetics.Ea / reaction.kinetics.Ea, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.Ea, reaction.kinetics.Ea))
        self.assertAlmostEqual(self.reaction.kinetics.T0 / reaction.kinetics.T0, 1.0, 6, '{0} != {1} within 6 places'.format(self.reaction.kinetics.T0, reaction.kinetics.T0))

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
