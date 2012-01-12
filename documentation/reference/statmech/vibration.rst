******************************
Vibrational degrees of freedom
******************************

For most molecular systems, a quantum treatment of vibrational motion
is required since typical vibrational energies are of similar magnitude to
:math:`k_\mathrm{B} T`.

Harmonic oscillators
====================

The internal vibrational motions of a molecule are often well described as a
set of independent one-dimensional harmonic oscillators.

The time-independent Schrödinger equation for a simple harmonic oscillator is

.. math:: 

    \hat{H} \Psi(\theta, \phi)
    = 
    -\frac{\hbar^2}{2 n} \frac{\partial^2}{\partial x^2} \Psi(x) + V(x) \Psi(x)
    =
    E \Psi(\theta, \phi)

where :math:`m` is the mass of the oscillator and the potential is given by

.. math:: V(x) = \frac{1}{2} m \omega^2 x^2

where :math:`\omega` is the vibrational frequency. The solution of the
above gives energy states

.. math:: E_n = \hbar \omega \left( n + \frac{1}{2} \right) \hspace{2em} n = 0, 1, 2, \ldots

which implies that the ground state energy :math:`E_0 = \hbar \omega / 2` is
nonzero. Typical molecular vibrational frequencies corresponding to energies
that are significantly spaced when compared to :math:`k_\mathrm{B} T`. This
means that a full quantum treatment of vibrations is required.

The quantum partition function is easy to evaluate by summing the energy states:

.. math:: Q(T) = \frac{1}{1 - e^{-\xi}}

Above we have introduced a dimensionless variable :math:`\xi` defined as

.. math:: \xi \equiv \frac{\hbar \omega}{k_\mathrm{B} T}

Note that, by choosing the above form for the partition function, we are
choosing that our energies are referenced to the ground-state vibrational 
energy, not to the bottom of the vibrational potential well. If the latter is
preferred, the partition function gains an additional :math:`e^{-\xi/2}` in the
numerator.

The vibrational contributions to heat capacity, enthalpy, and entropy are
easy to evaluate from the quantum partition function:

.. math:: \frac{C_\mathrm{v}(T)}{R} = \xi^2 \frac{e^{\xi}}{\left( 1 - e^{\xi} \right)^2}

.. math:: \frac{H(T)}{RT} = \frac{\xi}{e^{\xi} - 1}

.. math:: \frac{S(T)}{R} = \left[ - \ln \left(1 - e^{-\xi} \right) + \frac{\xi}{e^{\xi} - 1} \right]

The sum and density of states for a quantum harmonic oscillator are summations
involving the Heaviside step function :math:`\Theta(x)` or Dirac delta function
:math:`\delta(x)`, respectively. As with the partition function, these have
been referenced such that the zero of energy resides at the vibrational ground
state, not the bottom of the vibrational potential well:

.. math:: N(E) = \sum_{n=0}^\infty \Theta(E - n \omega)

.. math:: \rho(E) = \sum_{n=0}^\infty \delta(E - n \omega)

Fortunately, there is a very efficient algorithm for numerically convolving a
harmonic oscillator into an existing density or sum of states: the 
Beyer-Swinehart algorithm. This is much more efficient than evaluating the sum
or density of states using the above and convolving in a separate step, and
if the input sum or density of states is smooth, then the output will also be
smooth.

The presence of closed-form analytical results or efficient algorithms for 
working with the quantum mechanical model suggest that there is no need for
the classical model. Nonetheless, we provide it here for completeness,
recognizing that it gives the same results as the quantum treatment at very
high temperatures:

.. math:: Q(T \rightarrow \infty) \rightarrow \frac{k_\mathrm{B} T}{\hbar \omega} = \frac{1}{\xi}

.. math:: \frac{C_\mathrm{v}(T \rightarrow \infty)}{R} \rightarrow 1

.. math:: \frac{H(T \rightarrow \infty)}{RT} \rightarrow 1

.. math:: \frac{S(T \rightarrow \infty)}{R} \rightarrow \ln Q(T) + 1

.. math:: N(E) \rightarrow \frac{E}{N_\mathrm{vib}! \prod_{i=1}^{N_\mathrm{vib}} \hbar \omega_i}

.. math:: \rho(E) \rightarrow \frac{1}{(N_\mathrm{vib} - 1)! \prod_{i=1}^{N_\mathrm{vib}} \hbar \omega_i}


.. autoclass:: pyrate.statmech.HarmonicOscillator
    :members:
