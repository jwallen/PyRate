********************************
Translational degrees of freedom
********************************

Translational energies are much smaller than :math:`k_\mathrm{B} T` except for
temperatures approaching absolute zero, so a classical treatment of translation
is more than adequate.

Ideal gas translation
=====================

In this section we will consider the translational motion of an ideal gas in
three dimensions. The time-independent Schrödinger equation for this system is
that of a particle in a three-dimensional box of size :math:`L_x \times L_y \times L_z`:

.. math:: 

    \hat{H} \Psi(x,y,z)
    = 
    \left[ -\frac{\hbar^2}{2m} \left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + \frac{\partial^2}{\partial z^2} \right) + V(\mathbf{x,y,z}) \right] \Psi(x,y,z)
    =
    E \Psi(x,y,z)

where :math:`m` is the mass of the ideal gas molecule and the potential is
the infinite square well potential:

.. math::

    V(\mathbf{r}) = \begin{cases}
    0 & 0 < x < L_x, 0 < y < L_y, 0 < z < L_z \\
    \infty & \mathrm{otherwise}
    \end{cases}

The solution of the above gives energy states described by three quantum numbers
:math:`n_x`, :math:`n_y`, and :math:`n_z`:

.. math:: E_{(n_x,n_y,n_z)} = \frac{\hbar^2}{2m} \left( \frac{n_x^2 \pi^2}{L_x^2} + \frac{n_y^2 \pi^2}{L_y^2} + \frac{n_z^2 \pi^2}{L_z^2} \right) \hspace{2em} n_x, n_y, n_z = 1, 2, \ldots

For atomic and molecular masses and length scales, these energy states are
spaced very closely together relative to :math:`k_\mathrm{B} T` except at
temperatures approaching absolute zero. For this reason, we commonly evaluate
the partition function in the classical limit:

.. math:: Q(T) = \left( \frac{2 \pi m k_\mathrm{B} T}{h^2} \right)^{3/2} \frac{k_\mathrm{B} T}{P}

Note that we have also replaced the volume :math:`V = L_x L_y L_z` by the
expression :math:`k_\mathrm{B} T/P` in accordance with the ideal gas law.
Since the partition function has temperature dependence of :math:`T^{5/2}`,
the sum and density of states are easy to evaluate by inverse Laplace transform:

.. math:: N(E) = \left( \frac{2 \pi m}{h^2} \right)^{3/2} \frac{E^{5/2}}{\Gamma(7/2)} \frac{1}{P}

.. math:: \rho(E) = \left( \frac{2 \pi m}{h^2} \right)^{3/2} \frac{E^{3/2}}{\Gamma(5/2)} \frac{1}{P}

The contributions to the heat capacity, enthalpy, and entropy by translation
can also be derived from the partition function:

.. math:: \frac{C_\mathrm{v}(T)}{R} = \frac{3}{2}

.. math:: \frac{H(T)}{RT} = \frac{3}{2}

.. math:: \frac{S(T)}{R} = \ln Q(T) + \frac{3}{2} + 1

In PyRate the translation of an ideal gas is modeled using the :class:`IdealGasTranslation`
class. The only parameter required is the mass of the translating molecule.

.. autoclass:: pyrate.statmech.IdealGasTranslation
    :members:
