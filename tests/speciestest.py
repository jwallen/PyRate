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
This script contains unit tests of the :mod:`pyrate.species` module.
"""

import unittest
import math

from pyrate.species import *

################################################################################

class TestLennardJones(unittest.TestCase):
    """
    Contains unit tests of the :class:`LennardJones` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.sigma = pq.Quantity(3.70, "angstrom")
        self.epsilon = pq.Quantity(94.9, "K")
        self.lennardJones = LennardJones(
            sigma = self.sigma,
            epsilon = self.epsilon,
        )

    def test_sigma(self):
        """
        Test that the Lennard-Jones sigma property was properly set.
        """
        self.assertAlmostEqual(self.lennardJones.sigma / self.sigma, 1.0, 6, '{0} != {1} within 6 places'.format(self.lennardJones.sigma, self.sigma))

    def test_epsilon(self):
        """
        Test that the Lennard-Jones sigma property was properly set.
        """
        self.assertAlmostEqual(self.lennardJones.epsilon / self.epsilon, 1.0, 6, '{0} != {1} within 6 places'.format(self.lennardJones.epsilon, self.epsilon))

    def test_pickle(self):
        """
        Test that a Lennard-Jones object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        lennardJones = cPickle.loads(cPickle.dumps(self.lennardJones))
        self.assertAlmostEqual(self.lennardJones.sigma / lennardJones.sigma, 1.0, 6, '{0} != {1} within 6 places'.format(self.lennardJones.sigma, lennardJones.sigma))
        self.assertAlmostEqual(self.lennardJones.epsilon / lennardJones.epsilon, 1.0, 6, '{0} != {1} within 6 places'.format(self.lennardJones.epsilon, lennardJones.epsilon))
    
    def test_repr(self):
        """
        Test that a Lennard-Jones object can be successfully reconstructed from
        its repr() output with no loss of information.
        """
        exec('lennardJones = {0!r}'.format(self.lennardJones))
        self.assertAlmostEqual(self.lennardJones.sigma / lennardJones.sigma, 1.0, 6, '{0} != {1} within 6 places'.format(self.lennardJones.sigma, lennardJones.sigma))
        self.assertAlmostEqual(self.lennardJones.epsilon / lennardJones.epsilon, 1.0, 6, '{0} != {1} within 6 places'.format(self.lennardJones.epsilon, lennardJones.epsilon))

################################################################################

class TestSpecies(unittest.TestCase):
    """
    Contains unit tests of the :class:`Species` class.
    """
    
    def setUp(self):
        """
        A function run before each unit test in this class.
        """
        self.label = 'N2'
        self.thermo = None
        self.statmech = None
        self.lennardJones = LennardJones(sigma=pq.Quantity(3.70, "angstrom"), epsilon=pq.Quantity(94.9, "K"))
        self.molecularWeight = pq.Quantity(28.01, "amu")
        self.collisionModel = None
        self.species = Species(
            label = self.label, 
            thermo = self.thermo, 
            statmech = self.statmech, 
            lennardJones = self.lennardJones, 
            molecularWeight = self.molecularWeight, 
            collisionModel = self.collisionModel,
        )
        
    def test_label(self):
        """
        Test that the Species label property was properly set.
        """
        self.assertEqual(self.species.label, self.label)

    def test_molecularWeight(self):
        """
        Test that the Species molecularWeight property was properly set.
        """
        self.assertAlmostEqual(self.species.molecularWeight / self.molecularWeight, 1.0, 6, '{0} != {1} within 6 places'.format(self.species.molecularWeight, self.molecularWeight))

    def test_pickle(self):
        """
        Test that a Species object can be successfully pickled and
        unpickled with no loss of information.
        """
        import cPickle
        species = cPickle.loads(cPickle.dumps(self.species))
        self.assertEqual(self.species.label, species.label)
        self.assertEqual(self.species.thermo, species.thermo)
        self.assertEqual(self.species.statmech, species.statmech)
        self.assertAlmostEqual(self.species.lennardJones.sigma / species.lennardJones.sigma, 1.0, 6, '{0} != {1} within 6 places'.format(self.species.lennardJones.sigma, species.lennardJones.sigma))
        self.assertAlmostEqual(self.species.lennardJones.epsilon / species.lennardJones.epsilon, 1.0, 6, '{0} != {1} within 6 places'.format(self.species.lennardJones.epsilon, species.lennardJones.epsilon))
        self.assertAlmostEqual(self.species.molecularWeight / species.molecularWeight, 1.0, 6, '{0} != {1} within 6 places'.format(self.species.molecularWeight, species.molecularWeight))
        self.assertEqual(self.species.collisionModel, species.collisionModel)
   
    def test_repr(self):
        """
        Test that a Species object can be successfully reconstructed from its
        repr() output with no loss of information.
        """
        exec('species = {0!r}'.format(self.species))
        self.assertEqual(self.species.label, species.label)
        self.assertEqual(self.species.thermo, species.thermo)
        self.assertEqual(self.species.statmech, species.statmech)
        self.assertAlmostEqual(self.species.lennardJones.sigma / species.lennardJones.sigma, 1.0, 6, '{0} != {1} within 6 places'.format(self.species.lennardJones.sigma, species.lennardJones.sigma))
        self.assertAlmostEqual(self.species.lennardJones.epsilon / species.lennardJones.epsilon, 1.0, 6, '{0} != {1} within 6 places'.format(self.species.lennardJones.epsilon, species.lennardJones.epsilon))
        self.assertAlmostEqual(self.species.molecularWeight / species.molecularWeight, 1.0, 6, '{0} != {1} within 6 places'.format(self.species.molecularWeight, species.molecularWeight))
        self.assertEqual(self.species.collisionModel, species.collisionModel)

################################################################################

if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
