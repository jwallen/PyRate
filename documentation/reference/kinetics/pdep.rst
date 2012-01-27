**********************************
Pressure-dependent kinetics models
**********************************

Pressure-dependent Arrhenius kinetics
=====================================

One approach to pressure-dependent kinetics is to let the Arrhenius parameters
themselves be a function of pressure:

.. math:: 
    k(T,P) = A(P) \left[ \frac{T}{T_0} \right]^{n(P)} \exp \left[ \frac{-E_\mathrm{a}(P)}{RT} \right]

In practice, this is usually implemented by fitting a set of Arrhenius
expressions at various pressures, then interpolating between them on a
logarithmic pressure scale. This functionality is provided by the 
:class:`PDepArrhenius` class. (Note that Chemkin refers to this as the ``PLOG``
format.)

.. autoclass:: pyrate.kinetics.PDepArrhenius
    :members:

Chebyshev polynomial
====================

A very general approach to pressure-dependent kinetics is to fit to polynomials
in temperature and pressure. One common choice for doing so is the set of
Chebyshev polynomials :math:`\phi_n(x)`, which are constructed to be orthogonal
to one another. The resulting expression for :math:`k(T,P)` is

.. math:: 
    \log k(T,P) = \sum_{t=1}^{N_T} \sum_{p=1}^{N_P} \alpha_{tp} \phi_t(\tilde{T}) \phi_p(\tilde{P})

where :math:`\alpha_{tp}` is the set of polynomial coefficients, 
:math:`\phi_n(x)` is the Chebyshev polynomial of degree :math:`n` evaluated at
:math:`x`, and

.. math:: \tilde{T} \equiv \frac{2T^{-1} - T_\mathrm{min}^{-1} - T_\mathrm{max}^{-1}}{T_\mathrm{max}^{-1} - T_\mathrm{min}^{-1}}

.. math:: \tilde{P} \equiv \frac{2 \log P - \log P_\mathrm{min} - \log P_\mathrm{max}}{\log P_\mathrm{max} - \log P_\mathrm{min}}

By choosing a suitable number of temperature and pressure terms :math:`N_T` and
:math:`N_P`, the complex behavior of the :math:`k(T,P)` rate coefficients can
be adequately approximated with a minimum of nonphysical behavior due to the
polynomial basis.

.. autoclass:: pyrate.kinetics.Chebyshev
    :members:
