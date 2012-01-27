************************************
Pressure-independent kinetics models
************************************

Arrhenius kinetics
==================

The Arrhenius equation is 

.. math::
    k(T) = A \left( \frac{T}{T_0} \right)^n \exp \left( -\frac{E_\mathrm{a}}{RT} \right)
        
where :math:`A` is the preexponential factor, :math:`T_0` is the reference
temperature, :math:`n` is the temperature exponent, and :math:`E_\mathrm{a}`
is the activation energy.

.. autoclass:: pyrate.kinetics.Arrhenius
    :members:
