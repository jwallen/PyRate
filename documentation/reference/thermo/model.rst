********************
Heat capacity models
********************

Wilhoit polynomial
==================

The Wilhoit polynomial is 

.. math::
    C_\mathrm{p}(T) = C_\mathrm{p}(0) + \left[ C_\mathrm{p}(\infty) -
    C_\mathrm{p}(0) \right] y^2 \left[ 1 + (y - 1) \sum_{i=0}^3 a_i y^i \right]
        
where :math:`y \equiv \frac{T}{T + B}` is a scaled temperature that ranges
from zero to one. This formulation has the advantage of correctly reproducing
the heat capacity behavior as :math:`T \rightarrow 0` and 
:math:`T \rightarrow \infty`. The low-temperature limit 
:math:`C_\mathrm{p}(0)` is taken to be :math:`3.5R` for linear molecules
and :math:`4R` for nonlinear molecules. The high-temperature limit 
:math:`C_\mathrm{p}(\infty)` is taken to be 
:math:`\left[ 3 N_\mathrm{atoms} - 1.5 \right] R` for linear molecules and
:math:`\left[ 3 N_\mathrm{atoms} - (2 + 0.5 N_\mathrm{rotors}) \right] R`
for nonlinear molecules, for a molecule composed of :math:`N_\mathrm{atoms}`
atoms and :math:`N_\mathrm{rotors}` internal rotors.

.. autoclass:: pyrate.thermo.Wilhoit
    :members:

NASA polynomial
===============

The NASA polynomial enables computation of the heat capacity, enthalpy, and
entropy via a set of dimensionless coefficients:

.. math::

    \frac{C_\mathrm{p}(T)}{R} = c_{-2} T^{-2} + c_{-1} T^{-1} + c_0 + c_1 T + c_2 T^2 + c_3 T^3 + c_4 T^4
        
.. math:: 

    \frac{H(T)}{RT} = - c_{-2} T^{-2} + c_{-1} T^{-1} \ln T + c_0 + \frac{1}{2} c_1 T + \frac{1}{3} c_2 T^2 + \frac{1}{4} c_3 T^3 + \frac{1}{5} c_4 T^4 + \frac{c_5}{T}
    
.. math:: 

    \frac{S(T)}{R} = -\frac{1}{2} c_{-2} T^{-2} - c_{-1} T^{-1} + c_0 \ln T + c_1 T + \frac{1}{2} c_2 T^2 + \frac{1}{3} c_3 T^3 + \frac{1}{4} c_4 T^4 + c_6

The above expressions are very fast to evaluate, as they are simple polynomial
fits. To improve the accuracy over a wide range of temperatures, the temperature
domain is often subdivided and separate NASA polynomials fit within each region.
The :class:`MultiNASA` class is provided for this purpose.

.. autoclass:: pyrate.thermo.NASA
    :members:

.. autoclass:: pyrate.thermo.MultiNASA
    :members:

Conversion between heat capacity models
=======================================

.. autofunction:: pyrate.thermo.convertWilhoitToMultiNASA

.. autofunction:: pyrate.thermo.convertMultiNASAToWilhoit
