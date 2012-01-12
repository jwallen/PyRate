****************************
Torsional degrees of freedom
****************************

Torsional modes have both vibrational and rotational character, and
therefore should be treated semiclassically or quantum mechanically.

Hindered rotors
===============

Some internal molecular degrees of freedom reflect torsional motion: rotation
around a bond. There is often a significant potential barrier to this
rotation; in this limit the hindered rotor appears to be a harmonic oscillator.
However, if the potential barrier to rotation is small compared to 
:math:`k_\mathrm{B} T`, the hindered rotor appears to be a free rotor. Thus,
a hindered rotor has mostly vibrational character at low temperatures and
mostly rotational character at high temperatures. This can be modeled either
classically or quantum mechanically.

The time-independent Schrödinger equation for a one-dimensional hindered rotor 
is

.. math:: 

    \hat{H} \Psi(\phi)
    = 
    -\frac{\hbar^2}{2I} \frac{d^2}{d \phi^2} \Psi(\phi) + V(\phi) \Psi(\phi)
    =
    E \Psi(\phi)

where :math:`I` is the (reduced) moment of inertia for the hindered rotation.
There are two common expressions used for the potential. For simple torsions
a cosine expression of the form

.. math:: V(\phi) = \frac{1}{2} V_0 \left[1 - \cos \left( \sigma \phi \right) \right]

is used, where :math:`V_0` is the height of the potential barrier and 
:math:`\sigma` is the number of minima or maxima in one revolution of angle
:math:`\phi`, equivalent to the symmetry number of that rotor. However, many
potentials are much more complex, and require fitting to a Fourier series of
the form

.. math:: V(\phi) = A + \sum_{k=1}^C \left( a_k \cos k \phi + b_k \sin k \phi \right)

Both potentials are always expressed such that :math:`V(0) = 0`.

The solution method utilizes an orthonormal basis set expansion of the form

.. math:: \Psi(\phi) = \sum_{m=-M}^M c_m \frac{e^{im\phi}}{\sqrt{2 \pi}}

Inserting this into the Schrödinger equation leads to a standard eigenvalue
problem that can be solved numerically for the energy states :math:`E`.
Unfortunately there is no closed-form analytical expression for these energy
levels, so the partition function and all derived thermodynamic and statistical
mechanical properties must also be evaluated numerically.

In the classical limit, however, we can obtain analytical expressions for
the cosine potential.

.. math:: Q^\mathrm{class}(T) = \left( \frac{2 \pi I k_\mathrm{B} T}{h^2} \right)^{1/2} \frac{2 \pi}{\sigma} \exp \left( -\frac{V_0}{2 k_\mathrm{B} T} \right) I_0 \left( \frac{V_0}{2 k_\mathrm{B} T} \right)

Above, :math:`I_0(x)` is the modified Bessel function of order zero for argument
:math:`x` and :math:`\sigma` is the symmetry number of the rotor. Often a 
semiclassical correction of the form

.. math:: Q(T) = \frac{Q_\mathrm{vib}^\mathrm{quant}(T)}{Q_\mathrm{vib}^\mathrm{class}(T)} Q^\mathrm{class}(T)

is applied, which gives

.. math:: Q(T) = \frac{h \nu}{k_\mathrm{B} T} \frac{1}{1 - \exp \left(- h \nu / k_\mathrm{B} T \right)} \left( \frac{2 \pi I k_\mathrm{B} T}{h^2} \right)^{1/2} \frac{2 \pi}{\sigma} \exp \left( -\frac{V_0}{2 k_\mathrm{B} T} \right) I_0 \left( \frac{V_0}{2 k_\mathrm{B} T} \right)

where

.. math:: \nu = \frac{\sigma}{2 \pi} \sqrt{\frac{V_0}{2 I}}

is the effective vibrational frequency of the hindered rotor. The corresponding
sum and density of states expressions (without the semiclassical correction) are

.. math:: 

    N(E) = \begin{cases}
    \frac{4 q_\mathrm{1f} V_0^{1/2}}{\pi^{3/2}} \left[ \mathcal{E}(E / V_0) - \left(1 - \frac{E}{V_0} \right) \mathcal{K}(E / V_0) \right] & E < V_0 \\
    \frac{4 q_\mathrm{1f} E^{1/2}}{\pi^{3/2}} \mathcal{E}(V_0 / E) & E > V_0
    \end{cases}

.. math:: 

    \rho(E) = \begin{cases}
    \frac{2 q_\mathrm{1f}}{\pi^{3/2} V_0^{1/2}} \mathcal{K}(E / V_0) & E < V_0 \\
    \frac{2 q_\mathrm{1f}}{\pi^{3/2} E^{1/2}} \mathcal{K}(V_0 / E) & E > V_0
    \end{cases}

where

.. math:: q_\mathrm{1f} = \frac{\pi^{1/2}}{\sigma} \left( \frac{8 \pi^2 I}{h^2} \right)^{1/2}

and :math:`\mathcal{K}(x)` and :math:`\mathcal{E}(x)` are the complete 
elliptic integrals of the first and second kind, respectively. The 
corresponding expression for heat capacity is

.. math:: \frac{C_\mathrm{v}(T)}{R} = \frac{C_\mathrm{v,vib}(T)}{R} -\frac{1}{2} + \zeta^2 - \left[ \zeta \frac{I_1(\zeta)}{I_0(\zeta)} \right]^2 - \zeta \frac{I_1(\zeta)}{I_0(\zeta)}

where we have defined the dimensionless parameter :math:`\zeta` such that

.. math:: \zeta \equiv \frac{V_0}{2 k_\mathrm{B} T}

A one-dimensional hindered rotor is therefore characterized by its reduced
moment of inertia :math:`I`, its symmetry number :math:`\sigma`, and either
its barrier height :math:`V_0` for the cosine potential or set of Fourier
series coefficients for the Fourier series potential.

.. autoclass:: pyrate.statmech.HinderedRotor
    :members:
