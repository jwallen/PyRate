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
This module contains Cython extension types that represent various models of
chemical reaction kinetics. Both pressure-independent and pressure-dependent 
models are available.
"""

import numpy
import quantities as pq
from libc.math cimport exp, log, sqrt, log10, cos, acos

cimport pyrate.constants as constants
import pyrate.units as units

################################################################################

cpdef getRateCoefficientUnitsFromReactionOrder(order):
    """
    Given a reaction `order`, return the corresponding SI units of the rate
    coefficient. These are the units that rate coefficients are stored in
    internally, as well as the units of the rate coefficient obtained using
    the ``simplified`` attribute of a :class:`Quantity` object that represents
    a rate coefficient. Raises a :class:`ValueError` if the units could not be
    determined.
    """
    if order == 0: 
        kunits = 'mol/(cm^3*s)'
    elif order == 1:
        kunits = 's^-1'
    elif order == 2:
        kunits = 'cm^3/(mol*s)'
    elif order == 3:
        kunits = 'cm^6/(mol^2*s)'
    elif order == 4:
        kunits = 'cm^9/(mol^3*s)'
    else:
        raise ValueError('Invalid reaction order {0}.'.format(order))
    return kunits

cpdef getReactionOrderFromRateCoefficientUnits(kunits):
    """
    Given a set of rate coefficient units `kunits`, return the corresponding
    reaction order. Raises a :class:`ValueError` if the reaction order could 
    not be determined.
    """
    dimensionality = pq.Quantity(1.0, kunits).simplified.dimensionality
    order = -1
    if len(dimensionality) == 3 and pq.mol in dimensionality and pq.m in dimensionality and pq.s in dimensionality:
        if dimensionality[pq.s] == -1 and dimensionality[pq.m] == -3 and dimensionality[pq.mol] == 1:
            order = 0
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 3 and dimensionality[pq.mol] == -1:
            order = 2
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 6 and dimensionality[pq.mol] == -2:
            order = 3
        elif dimensionality[pq.s] == -1 and dimensionality[pq.m] == 9 and dimensionality[pq.mol] == -3:
            order = 4
    elif len(dimensionality) == 1 and pq.s in dimensionality:
        if dimensionality[pq.s] == -1:
            order = 1
    if order == -1:
        raise ValueError('Invalid rate coefficient units "{0}".'.format(str(dimensionality)))
    return order

################################################################################

cdef class KineticsModel:
    """
    A base class for kinetics models, containing several attributes
    common to all models:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `Pmin`          The minimum pressure at which the model is valid, or zero if unknown or undefined
    `Pmax`          The maximum pressure at which the model is valid, or zero if unknown or undefined
    `order`         The reaction order (usually the number of reactants) 
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, Tmin=None, Tmax=None, Pmin=None, Pmax=None, order=-1, comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.order = order
        self.comment = comment
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        KineticsModel object.
        """
        Tmin = '({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        Tmax = '({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        Pmin = '({0:g},"{1}")'.format(float(self.Pmin), str(self.Pmin.dimensionality))
        Pmax = '({0:g},"{1}")'.format(float(self.Pmax), str(self.Pmax.dimensionality))
        return 'KineticsModel(Tmin={0}, Tmax={1}, Pmin={2}, Pmax={3}, order={4:d}, comment="""{4}""")'.format(Tmin, Tmax, Pmin, Pmax, self.order, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a KineticsModel object.
        """
        return (KineticsModel, (self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.order, self.comment))

    property Tmin:
        """The minimum temperature at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Tmin, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Tmin = 0.0 
            else:
                self._Tmin = float(units.convertTemperature(value, pq.K))

    property Tmax:
        """The maximum temperature at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Tmax, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Tmax = 0.0 
            else:
                self._Tmax = float(units.convertTemperature(value, pq.K))

    property Pmin:
        """The minimum pressure at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Pmin * 1e-5, pq.bar)
        def __set__(self, value):
            if value is None or value == 0:
                self._Pmin = 0.0 
            else:
                self._Pmin = float(units.convertPressure(value, pq.Pa))

    property Pmax:
        """The maximum pressure at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Pmax * 1e-5, pq.bar)
        def __set__(self, value):
            if value is None or value == 0:
                self._Pmax = 0.0 
            else:
                self._Pmax = float(units.convertPressure(value, pq.Pa))

    cpdef bint isTemperatureValid(self, double T) except -2:
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the temperature model, or ``False`` if not. If
        the minimum and maximum temperature are not defined, ``True`` is 
        returned.
        """
        return (self._Tmin == 0.0 or self._Tmin <= T) and (self._Tmax == 0.0 or T <= self._Tmax)

    cpdef bint isPressureValid(self, double P) except -2:
        """
        Return ``True`` if the pressure `P` in Pa is within the valid pressure
        range of the kinetics model, or ``False`` if not. If the minimum and 
        maximum pressure are not defined, ``True`` is returned.
        """
        return (self._Pmin == 0.0 or self._Pmin <= P) and (self._Pmax == 0.0 or P <= self._Pmax)

    cpdef bint isPressureDependent(self) except -2:
        """
        Returns ``True`` if the kinetics model represents pressure-dependent
        kinetics or ``False`` otherwise. This method must be overloaded in the 
        derived class.
        """
        raise NotImplementedError('Unexpected call to KineticsModel.isPressureDependent(); you should be using a class derived from KineticsModel that overrides this method.')

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of cm^3, 
        mol, and s at temperature `T` in K and, if required, pressure `P` in 
        bar. This method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to KineticsModel.getRateCoefficient(); you should be using a class derived from KineticsModel that overrides this method.')

################################################################################

cdef class Arrhenius(KineticsModel):
    """
    A kinetics model based on the (modified) Arrhenius equation. The attributes
    are:

    =========== ================================================================
    Attribute   Description
    =========== ================================================================
    `A`         The preexponential factor
    `T0`        The reference temperature
    `n`         The temperature exponent
    `Ea`        The activation energy
    =========== ================================================================
    
    """
    
    def __init__(self, A=0.0, n=0.0, Ea=0.0, T0=(1.0,"K"), Tmin=None, Tmax=None, order=-1, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, order=order, comment=comment)
        self.A = A
        self.n = n
        self.Ea = Ea
        self.T0 = T0
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Arrhenius object.
        """
        A = '({0:g},"{1}")'.format(float(self.A), str(self.A.dimensionality))
        Ea = '({0:g},"{1}")'.format(float(self.Ea), str(self.Ea.dimensionality))
        T0 = '({0:g},"{1}")'.format(float(self.T0), str(self.T0.dimensionality))
        string = 'Arrhenius(A={0}, n={1:g}, Ea={2}, T0={3}'.format(A, self.n, Ea, T0)
        if self.Tmin != 0.0: string += ', Tmin=({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        if self.Tmax != 0.0: string += ', Tmax=({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an Arrhenius object.
        """
        return (Arrhenius, (self.A, self.n, self.Ea, self.T0, self.Tmin, self.Tmax, self.order, self.comment))

    property A:
        """The preexponential factor."""
        def __get__(self):
            kunits = getRateCoefficientUnitsFromReactionOrder(self.order)
            return pq.Quantity(self._A, kunits)
        def __set__(self, value):
            if value is None or value == 0:
                self._A = 0.0 
            else:
                if isinstance(value, tuple):
                    value = pq.Quantity(value[0], value[1])
                # Try to determine reaction order from units of value
                order = getReactionOrderFromRateCoefficientUnits(value.units)
                # If the kinetics already has a specified reaction order,
                # make sure the units of A are consistent with that order
                if self.order != -1 and self.order != order:
                    raise ValueError('Units of preexponential factor "{0}" with reaction order {1:d} do not match internal reaction order {2:d}.'.format(str(value.units.dimensionality), order, self.order))
                self.order = order
                kunits = getRateCoefficientUnitsFromReactionOrder(order)
                self._A = float(units.convertRateCoefficient(value, kunits))

    property T0:
        """The reference temperature."""
        def __get__(self):
            return pq.Quantity(self._T0, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._T0 = 0.0 
            else:
                self._T0 = float(units.convertTemperature(value, pq.K))

    property Ea:
        """The activation energy."""
        def __get__(self):
            return pq.Quantity(self._Ea * 0.001, "kJ/mol")
        def __set__(self, value):
            if value is None or value == 0:
                self._Ea = 0.0 
            else:
                self._Ea = float(units.convertEnergy(value, "J/mol"))

    cpdef bint isPressureDependent(self) except -2:
        """
        Return ``False`` since Arrhenius kinetics are not pressure-dependent.
        """
        return False

    cpdef double getRateCoefficient(self, double T, double P=0.0) except -1:
        """
        Return the rate coefficient in the appropriate combination of cm^3, 
        mol, and s at temperature `T` in K. 
        """
        return self._A * (T / self._T0)** self.n * exp(-self._Ea / (constants.R * T))

    cpdef changeT0(self, double T0):
        """
        Changes the reference temperature used in the exponent to `T0` in K, 
        and adjusts the preexponential factor accordingly.
        """
        self._A /= (self._T0 / T0)**self.n
        self._T0 = T0

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray klist, kunits, double T0=1, numpy.ndarray weights=None, bint threeParams=True):
        """
        Fit the Arrhenius parameters to a set of rate coefficient data `klist`
        in units of `kunits` corresponding to a set of temperatures `Tlist` in 
        K. A linear least-squares fit is used, which guarantees that the 
        resulting parameters provide the best possible approximation to the 
        data.
        """
        import numpy.linalg
        import scipy.stats
        if threeParams:
            A = numpy.zeros((len(Tlist),3), numpy.float64)
            A[:,0] = numpy.ones_like(Tlist)
            A[:,1] = numpy.log(Tlist / T0)
            A[:,2] = -1.0 / constants.R / Tlist
        else:
            A = numpy.zeros((len(Tlist),2), numpy.float64)
            A[:,0] = numpy.ones_like(Tlist)
            A[:,1] = -1.0 / constants.R / Tlist
        b = numpy.log(klist)
        if weights is not None:
            for n in range(b.size):
                A[n,:] *= weights[n]
                b[n] *= weights[n]
        x, residues, rank, s = numpy.linalg.lstsq(A,b)
        
        # Determine covarianace matrix to obtain parameter uncertainties
        count = klist.size
        cov = residues[0] / (count - 3) * numpy.linalg.inv(numpy.dot(A.T, A))
        t = scipy.stats.t.ppf(0.975, count - 3)
            
        if not threeParams:
            x = numpy.array([x[0], 0, x[1]])
            cov = numpy.array([[cov[0,0], 0, cov[0,1]], [0,0,0], [cov[1,0], 0, cov[1,1]]])
        
        self.A = (exp(x[0]), kunits)
        self.n = x[1]
        self.Ea = (x[2], "J/mol")
        self.T0 = (T0, "K")
        self.Tmin = (numpy.min(Tlist), "K")
        self.Tmax = (numpy.max(Tlist), "K")
        self.comment = 'Fitted to {0:d} data points; dA = *|/ {1:g}, dn = +|- {2:g}, dEa = +|- {3:g} kJ/mol'.format(
            len(Tlist),
            exp(sqrt(cov[0,0])),
            sqrt(cov[1,1]),
            sqrt(cov[2,2]) * 0.001,
        )
        
        return self

################################################################################

cdef class PDepArrhenius(KineticsModel):
    """
    A kinetic model of a phenomenological rate coefficient :math:`k(T,P)` where
    a set of Arrhenius kinetics are stored at a variety of pressures and
    interpolated between on a logarithmic scale. The attributes are:

    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `pressures`     The list of pressures
    `arrhenius`     The list of :class:`Arrhenius` objects at each pressure
    =============== ============================================================
    
    """

    def __init__(self, pressures=None, arrhenius=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, order=-1, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, order=order, comment=comment)
        self.pressures = pressures
        self.arrhenius = arrhenius or []

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        PDepArrhenius object.
        """
        pressures = '([{0}],"{1}")'.format(','.join(['{0:g}'.format(float(P)) for P in self.pressures]), str(self.pressures.dimensionality))
        string = 'PDepArrhenius(pressures={0}, arrhenius={1!r}'.format(pressures, self.arrhenius)
        if self.Tmin != 0.0: string += ', Tmin=({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        if self.Tmax != 0.0: string += ', Tmax=({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        if self.Pmin != 0.0: string += ', Pmin=({0:g},"{1}")'.format(float(self.Pmin), str(self.Pmin.dimensionality))
        if self.Pmax != 0.0: string += ', Pmax=({0:g},"{1}")'.format(float(self.Pmax), str(self.Pmax.dimensionality))
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a PDepArrhenius object.
        """
        return (PDepArrhenius, (self.pressures, self.arrhenius, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.order, self.comment))

    property pressures:
        """The list of pressures."""
        def __get__(self):
            return pq.Quantity(self._pressures * 1e-5, pq.bar)
        def __set__(self, value):
            if value is None:
                self._pressures = numpy.array([]) 
            else:
                self._pressures = numpy.array(units.convertPressure(value, pq.Pa))
    
    cpdef bint isPressureDependent(self) except -2:
        """
        Returns ``True`` since PDepArrhenius kinetics are pressure-dependent.
        """
        return True

    cdef getAdjacentExpressions(self, double P):
        """
        Returns the pressures and Arrhenius expressions for the pressures that
        most closely bound the specified pressure `P` in bar.
        """
        cdef Arrhenius arrh
        cdef int i, ilow, ihigh
        P = P * 1e5
        
        if P in self._pressures:
            arrh = self.arrhenius[[p for p in self._pressures].index(P)]
            return P, P, arrh, arrh
        else:
            ilow = 0; ihigh = -1
            for i in range(1, len(self._pressures)):
                if self._pressures[i] <= P:
                    ilow = i
                if self._pressures[i] > P and ihigh == -1:
                    ihigh = i
            return self._pressures[ilow], self._pressures[ihigh], self.arrhenius[ilow], self.arrhenius[ihigh]
    
    cpdef double getRateCoefficient(self, double T, double P=0) except -1:
        """
        Return the rate coefficient in the appropriate combination of cm^3, 
        mol, and s at temperature `T` in K and pressure `P` in bar.
        
        If k(P+) and k(P-) (the values on either side of the requested P) are
        zero, then zero is returned. If only k(P-)==0 then it is replaced with 
        k(P+)/1e10, and vice versa. This allows the logarithmic interpolation
        to proceed without zero-division errors. (The expression is not defined
        when one of them is zero and the other is not, so we have to assume
        something.)
        """
        cdef double Plow, Phigh, klow, khigh, k
        cdef Arrhenius alow, ahigh
        cdef int j
        
        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent PDepArrhenius.getRateCoefficient().')
        
        k = 0.0
        Plow, Phigh, alow, ahigh = self.getAdjacentExpressions(P)
        if Plow == Phigh:
            k = alow.getRateCoefficient(T)
        else:
            P = P * 1e5
            klow = alow.getRateCoefficient(T)
            khigh = ahigh.getRateCoefficient(T)
            if klow == khigh == 0.0: return 0.0
            if klow == 0: klow = khigh/1e10
            if khigh == 0: khigh = klow/1e10
            k = klow * 10**(log10(P/Plow)/log10(Phigh/Plow)*log10(khigh/klow))
        return k
    
    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K, kunits, double T0=1):
        """
        Fit the pressure-dependent Arrhenius model to a matrix of rate
        coefficient data `K` with units of `kunits` corresponding to a set of 
        temperatures `Tlist` in K and pressures `Plist` in Pa. An Arrhenius 
        model is fit at each pressure.
        """
        cdef int i
        self.pressures = (Plist,"bar")
        self.arrhenius = []
        for i in range(len(Plist)):
            arrhenius = Arrhenius().fitToData(Tlist, K[:,i], kunits, T0)
            self.arrhenius.append(arrhenius)
        return self

################################################################################

cdef class Chebyshev(KineticsModel):
    """
    A model of a phenomenological rate coefficient :math:`k(T,P)` using a
    set of Chebyshev polynomials in temperature and pressure. The attributes
    are:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `coeffs`        Matrix of Chebyshev coefficients
    `degreeT`       The number of terms in the inverse temperature direction
    `degreeP`       The number of terms in the log pressure direction
    =============== ============================================================
    
    """

    def __init__(self, coeffs=None, Tmin=None, Tmax=None, Pmin=None, Pmax=None, comment=''):
        KineticsModel.__init__(self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment)
        self.coeffs = coeffs
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Chebyshev object.
        """
        coeffs = self.coeffs
        s = '(['
        for i in range(self.degreeT):
            if i > 0: s += ', '
            s += '[{0}]'.format(','.join(['{0:g}'.format(float(coeffs[i,j])) for j in range(self.degreeP)]))
        s += '],"{0}")'.format(str(coeffs.dimensionality))
        
        string = 'Chebyshev(coeffs={0}'.format(s)
        if self.Tmin != 0.0: string += ', Tmin=({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        if self.Tmax != 0.0: string += ', Tmax=({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        if self.Pmin != 0.0: string += ', Pmin=({0:g},"{1}")'.format(float(self.Pmin), str(self.Pmin.dimensionality))
        if self.Pmax != 0.0: string += ', Pmax=({0:g},"{1}")'.format(float(self.Pmax), str(self.Pmax.dimensionality))
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Chebyshev object.
        """
        return (Chebyshev, (self.coeffs, self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment))

    property coeffs:
        """The matrix of Chebyshev coefficients."""
        def __get__(self):
            kunits = getRateCoefficientUnitsFromReactionOrder(self.order)
            return pq.Quantity(self._coeffs, kunits)
        def __set__(self, value):
            if value is None:
                self._coeffs = numpy.array([])
                self.degreeT = 0
                self.degreeP = 0
            else:
                if isinstance(value, tuple):
                    value = pq.Quantity(value[0], value[1])
                # Try to determine reaction order from units of value
                order = getReactionOrderFromRateCoefficientUnits(value.units)
                # If the kinetics already has a specified reaction order,
                # make sure the units of A are consistent with that order
                if self.order != -1 and self.order != order:
                    raise ValueError('Units of preexponential factor "{0}" with reaction order {1:d} do not match internal reaction order {2:d}.'.format(str(value.units.dimensionality), order, self.order))
                self.order = order
                kunits = getRateCoefficientUnitsFromReactionOrder(order)
                factor = float(units.convertRateCoefficient((1.0, value.units), kunits))
                self._coeffs = numpy.array(value)
                self._coeffs[0,0] += log10(factor)
                self.degreeT = self._coeffs.shape[0]
                self.degreeP = self._coeffs.shape[1]
        
    cpdef bint isPressureDependent(self):
        """
        Return ``True`` since Chebyshev polynomial kinetics are 
        pressure-dependent.
        """
        return True

    cdef double chebyshev(self, int n, double x):
        if n == 0:
            return 1
        elif n == 1:
            return x
        elif n == 2:
            return -1 + 2*x*x
        elif n == 3:
            return x * (-3 + 4*x*x)
        elif n == 4:
            return 1 + x*x*(-8 + 8*x*x)
        elif n == 5:
            return x * (5 + x*x*(-20 + 16*x*x))
        elif n == 6:
            return -1 + x*x*(18 + x*x*(-48 + 32*x*x))
        elif n == 7:
            return x * (-7 + x*x*(56 + x*x*(-112 + 64*x*x)))
        elif n == 8:
            return 1 + x*x*(-32 + x*x*(160 + x*x*(-256 + 128*x*x)))
        elif n == 9:
            return x * (9 + x*x*(-120 + x*x*(432 + x*x*(-576 + 256*x*x))))
        else:
            return cos(n * acos(x))

    cdef double getReducedTemperature(self, double T) except -1000:
        """
        Return the reduced temperature corresponding to the given temperature
        `T` in K. This maps the inverse of the temperature onto the domain 
        [-1, 1] using the `Tmin` and `Tmax` attributes as the limits.
        """
        return (2.0/T - 1.0/self._Tmin - 1.0/self._Tmax) / (1.0/self._Tmax - 1.0/self._Tmin)
    
    cdef double getReducedPressure(self, double P) except -1000:
        """
        Return the reduced pressure corresponding to the given pressure
        `P` in bar. This maps the logarithm of the pressure onto the domain 
        [-1, 1] using the `Pmin` and `Pmax` attributes as the limits.
        """
        return (2.0*log10(P*1e5) - log10(self._Pmin) - log10(self._Pmax)) / (log10(self._Pmax) - log10(self._Pmin))
    
    cpdef double getRateCoefficient(self, double T, double P=0) except -1:
        """
        Return the rate coefficient in the appropriate combination of cm^3, 
        mol, and s at temperature `T` in K and pressure `P` in bar by 
        evaluating the Chebyshev expression.
        """
        cdef double Tred, Pred, k
        cdef int i, j, t, p
        
        if P == 0:
            raise ValueError('No pressure specified to pressure-dependent Chebyshev.getRateCoefficient().')

        k = 0.0
        Tred = self.getReducedTemperature(T)
        Pred = self.getReducedPressure(P)
        for t in range(self.degreeT):
            for p in range(self.degreeP):
                k += self._coeffs[t,p] * self.chebyshev(t, Tred) * self.chebyshev(p, Pred)
        return 10.0**k

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K,
        str kunits, int degreeT, int degreeP, double Tmin, double Tmax, double Pmin, double Pmax):
        """
        Fit a Chebyshev kinetic model to a set of rate coefficients `K`, which
        is a matrix corresponding to the temperatures `Tlist` in K and pressures
        `Plist` in bar. `degreeT` and `degreeP` are the degree of the 
        polynomials in temperature and pressure, while `Tmin`, `Tmax`, `Pmin`,
        and `Pmax` set the edges of the valid temperature and pressure ranges
        in K and bar, respectively.
        """
        cdef int nT = len(Tlist), nP = len(Plist)
        cdef list Tred, Pred
        cdef int t1, p1, t2, p2
        cdef double T, P

        self.degreeT = degreeT; self.degreeP = degreeP

        # Set temperature and pressure ranges
        self.Tmin = (Tmin,"K")
        self.Tmax = (Tmax,"K")
        self.Pmin = (Pmin,"bar")
        self.Pmax = (Pmax,"bar")

        # Calculate reduced temperatures and pressures
        Tred = [self.getReducedTemperature(T) for T in Tlist]
        Pred = [self.getReducedPressure(P) for P in Plist]

        # Create matrix and vector for coefficient fit (linear least-squares)
        A = numpy.zeros((nT*nP, degreeT*degreeP), numpy.float64)
        b = numpy.zeros((nT*nP), numpy.float64)
        for t1, T in enumerate(Tred):
            for p1, P in enumerate(Pred):
                for t2 in range(degreeT):
                    for p2 in range(degreeP):
                        A[p1*nT+t1, p2*degreeT+t2] = self.chebyshev(t2, T) * self.chebyshev(p2, P)
                b[p1*nT+t1] = log10(K[t1,p1])

        # Do linear least-squares fit to get coefficients
        x, residues, rank, s = numpy.linalg.lstsq(A, b)

        # Extract coefficients
        coeffs = numpy.zeros((degreeT,degreeP), numpy.float64)
        for t2 in range(degreeT):
            for p2 in range(degreeP):
                coeffs[t2,p2] = x[p2*degreeT+t2]
        self.coeffs = (coeffs, kunits)
        
        return self
