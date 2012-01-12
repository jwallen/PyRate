*****************************
Rotational degrees of freedom
*****************************

For most molecular systems, a classical treatment of rotational motion is
sufficient since typical rotational energies are much smaller in magnitude
than :math:`k_\mathrm{B} T`.

Linear rotors
=============

The external rotation of diatomic and other linear molecules is often well
described using a linear rigid rotor.

The time-independent Schrödinger equation for a linear rigid rotor is

.. math:: 

    \hat{H} \Psi(\theta, \phi)
    = 
    -\frac{\hbar^2}{2I} \left[ \frac{1}{\sin \theta} \frac{\partial}{\partial \theta} \left( \sin \theta \frac{\partial}{\partial \theta} \right) + \frac{1}{\sin^2 \theta} \frac{\partial^2}{\partial^2 \phi} \right] \Psi(\theta, \phi)
    =
    E \Psi(\theta, \phi)

where :math:`I` is the moment of inertia of the molecule. The solution of the
above gives energy levels

.. math:: E_\ell = \frac{\hbar^2}{2I} \ell \left( \ell + 1 \right) \hspace{2em} \ell = 0, 1, 2, \ldots

with corresponding degeneracies

.. math:: g_\ell = 2 \ell + 1

The corresponding wavefunction solutions are the spherical harmonics 
:math:`Y_\ell^m(\theta, \phi)`, which are described by two quantum numbers:
:math:`\ell = 0, 1, 2, \ldots` and 
:math:`m = -\ell, -\ell + 1, \ldots, \ell - 1, \ell`.

In the majority of chemical systems, the moment of inertia :math:`I` is such
that the energies :math:`E_\ell` are much smaller than :math:`k_\mathrm{B} T`
at temperatures of practical interest. As a result, the contributions to
various thermodynamic and statistical mechanical properties from the linear
rigid rotor are usually evaluated in the classical limit.

The classical partition function for a linear rotor is

.. math:: Q(T) = \frac{8 \pi^2 I k_\mathrm{B} T}{\sigma h^2}

where :math:`\sigma` is the symmetry number of the rotation, which is used to
avoid multiple counting of indistinguishable states that arise due to symmetries
in the molecule. Since the classical partition function has temperature 
dependence of :math:`T^{5/2}`, the sum and density of states are easy to 
evaluate by inverse Laplace transform:

.. math:: N(E) = \frac{8 \pi^2 I}{\sigma h^2} E

.. math:: \rho(E) = \frac{8 \pi^2 I}{\sigma h^2}

The contributions to the heat capacity, enthalpy, and entropy by translation
can also be derived from the partition function:

.. math:: \frac{C_\mathrm{v}(T)}{R} = 1

.. math:: \frac{H(T)}{RT} = 1

.. math:: \frac{S(T)}{R} = \ln Q(T) + 1

A linear rigid rotor is thus characterized by its moment of inertia :math:`I`
and symmetry number :math:`\sigma`. An alternative representation for the
moment of inertia is as a *rotational constant* :math:`B`:

.. math:: B = \frac{\hbar^2}{2I}

.. autoclass:: pyrate.statmech.LinearRotor
    :members:

Nonlinear rotors
================

The external rotation of polyatomic nonlinear molecules can be described by
a similar procedure as for linear molecules. Here we call the model a nonlinear
rigid rotor to distinguish it from a linear rotor. The defining characteristic
of a nonlinear molecule is that it has three moments of inertia :math:`I_x`,
:math:`I_y`, and :math:`I_z` instead of just one. Although we will not use 
this distinction in ChemPy, there are several types of nonlinear molecules
depending on the relative values of math:`I_x`, :math:`I_y`, and :math:`I_z`:
spherical top, symmetric top, asymmetric top, etc. The equations that follow
work for any of these models.

A quantum treatment of the rotational motion of nonlinear molecules is quite
involved. Since rotational motion can usually be adequately described by
classical mechanics, the quantum mechanics is beyond the scope of this 
discussion. Here we will simply state the classical results, first for the
partition function:

.. math:: Q(T) = \frac{\sqrt{\pi}}{\sigma} \left( \frac{8 \pi^2 k_\mathrm{B} T}{h^2} \right)^{3/2} \sqrt{I_\mathrm{A} I_\mathrm{B} I_\mathrm{C}}

Note the difference in temperature dependence compared to the linear rotor.
As usual, the sum and density of states can be easily determined via inverse
Laplace transform of the partition function, with energy dependences that are
again different from the corresponding expressions for the linear rotor:

.. math:: N(E) = \frac{\sqrt{\pi}}{\sigma} \left( \frac{8 \pi^2}{h^2} \right)^{3/2} \sqrt{I_\mathrm{A} I_\mathrm{B} I_\mathrm{C}} \frac{E^{3/2}}{(3/2)!}

.. math:: \rho(E) = \frac{\sqrt{\pi}}{\sigma} \left( \frac{8 \pi^2}{h^2} \right)^{3/2} \sqrt{I_\mathrm{A} I_\mathrm{B} I_\mathrm{C}} \frac{E^{1/2}}{\frac{1}{2}!}

Finally, the contributions to the heat capacity, enthalpy, and entropy reflect
that we are modeling three modes of rotational motion (instead of two, as for
the linear rotor):

.. math:: \frac{C_\mathrm{v}(T)}{R} = \frac{3}{2}

.. math:: \frac{H(T)}{RT} = \frac{3}{2}

.. math:: \frac{S(T)}{R} = \ln Q(T) + \frac{3}{2}

The nonlinear rigid rotor is characterized by its three principal moments of
inertia math:`I_x`, :math:`I_y`, and :math:`I_z` and by its symmetry number
:math:`\sigma`. As with linear rotors, the moments of inertia are sometimes
expressed in terms of rotational constants.

.. autoclass:: pyrate.statmech.NonlinearRotor
    :members:

K-rotors
========

Some situations require the consideration of a so-called "active K-rotor",
which is typically modeled as a one-dimensional rigid rotor.

The solution of the time-independent Schrödinger equation for a one-dimensional
rigid rotor gives energy levels described by a single quantum number :math:`m`

.. math:: E_m = \frac{\hbar^2}{2I} m^2 \hspace{2em} m = 0, 1, 2, \ldots

with corresponding degeneracies

.. math:: 

    g_m = \begin{cases}
    1 & m = 0 \\
    2 & \mathrm{otherwise}
    \end{cases}

As with the linear and nonlinear rigid rotors, the typical rotational energy
levels are usually closely spaced compared to  :math:`k_\mathrm{B} T`
at temperatures of practical interest. As a result, the contributions to
various thermodynamic and statistical mechanical properties from the linear
rigid rotor are usually evaluated in the classical limit.

The classical partition function for a one-dimensional rotor is

.. math:: Q(T) = \frac{1}{\sigma} \sqrt{\frac{8 \pi I k_\mathrm{B} T}{h^2}}

where :math:`\sigma` is the symmetry number of the rotation, which is used to
avoid multiple counting of indistinguishable states that arise due to symmetries
in the molecule. Since the classical partition function has temperature 
dependence of :math:`T^{1/2}`, the sum and density of states are easy to 
evaluate by inverse Laplace transform:

.. math:: N(E) = \frac{2}{\sigma} \sqrt{\frac{8 \pi^2 I E}{\sigma h^2}}

.. math:: \rho(E) = \frac{1}{\sigma} \sqrt{\frac{8 \pi^2 I}{h^2 E}}

The contributions to the heat capacity, enthalpy, and entropy by translation
can also be derived from the partition function:

.. math:: \frac{C_\mathrm{v}(T)}{R} = \frac{1}{2}

.. math:: \frac{H(T)}{RT} = \frac{1}{2}

.. math:: \frac{S(T)}{R} = \ln Q(T) + \frac{1}{2}

A one-dimensional rigid rotor is thus characterized by its moment of inertia
:math:`I` (or rotational constant :math:`B`) and symmetry number 
:math:`\sigma`. 

.. autoclass:: pyrate.statmech.KRotor
    :members:
