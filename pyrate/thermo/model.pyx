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
This module contains Cython extension types that represent various
thermodynamic models of heat capacity.
"""

import numpy
import quantities as pq
from libc.math cimport log

cimport pyrate.constants as constants
import pyrate.units as units

################################################################################

cdef class HeatCapacityModel:
    """
    A base class for heat capacity models, containing several attributes
    common to all models:
    
    =============== ============================================================
    Attribute       Description
    =============== ============================================================
    `Tmin`          The minimum temperature at which the model is valid, or zero if unknown or undefined
    `Tmax`          The maximum temperature at which the model is valid, or zero if unknown or undefined
    `comment`       Information about the model (e.g. its source)
    =============== ============================================================

    """
    
    def __init__(self, Tmin=None, Tmax=None, comment=''):
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.comment = comment
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        ThermoModel object.
        """
        Tmin = '({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        Tmax = '({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        return 'ThermoModel(Tmin={0}, Tmax={1}, comment="""{2}""")'.format(Tmin, Tmax, self.comment)

    def __reduce__(self):
        """
        A helper function used when pickling a ThermoModel object.
        """
        return (HeatCapacityModel, (self.Tmin, self.Tmax, self.comment))

    property Tmin:
        """The minimum temperature at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Tmin, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Tmin = 0.0 
            else:
                self._Tmin = float(units.convertTemperature(value, "K"))

    property Tmax:
        """The maximum temperature at which the model is valid, or zero if unknown or undefined."""
        def __get__(self):
            return pq.Quantity(self._Tmax, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Tmax = 0.0 
            else:
                self._Tmax = float(units.convertTemperature(value, "K"))

    cpdef bint isTemperatureValid(self, double T) except -2:
        """
        Return ``True`` if the temperature `T` in K is within the valid
        temperature range of the thermodynamic data, or ``False`` if not. If
        the minimum and maximum temperature are not defined, ``True`` is 
        returned.
        """
        return (self._Tmin == 0.0 or self._Tmin <= T) and (self._Tmax == 0.0 or T <= self._Tmax)

    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at temperature 
        `T` in K. This method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to ThermoModel.get_heat_capacity(); you should be using a class derived from ThermoModel.')

    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in J/mol at temperature `T` in K. This method must 
        be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to ThermoModel.get_enthalpy(); you should be using a class derived from ThermoModel.')

    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at temperature `T` in K. This method must 
        be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to ThermoModel.get_entropy(); you should be using a class derived from ThermoModel.')

    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in J/mol at temperature `T` in K. This 
        method must be overloaded in the derived class.
        """
        raise NotImplementedError('Unexpected call to ThermoModel.get_free_energy(); you should be using a class derived from ThermoModel.')

################################################################################

cdef class Wilhoit(HeatCapacityModel):
    """
    A thermodynamics model based on the Wilhoit equation for heat capacity.
    """

    def __init__(self, Cp0=None, CpInf=None, a0=0.0, a1=0.0, a2=0.0, a3=0.0, H0=None, S0=None, B=None, Tmin=None, Tmax=None, comment=''):
        HeatCapacityModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.Cp0 = Cp0
        self.CpInf = CpInf
        self.B = B
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.H0 = H0
        self.S0 = S0
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        Wilhoit object.
        """
        Cp0 = '({0:g},"{1}")'.format(float(self.Cp0), str(self.Cp0.dimensionality))
        CpInf = '({0:g},"{1}")'.format(float(self.CpInf), str(self.CpInf.dimensionality))
        B = '({0:g},"{1}")'.format(float(self.B), str(self.B.dimensionality))
        H0 = '({0:g},"{1}")'.format(float(self.H0), str(self.H0.dimensionality))
        S0 = '({0:g},"{1}")'.format(float(self.S0), str(self.S0.dimensionality))
        string = 'Wilhoit(Cp0={0}, CpInf={1}, a0={2:g}, a1={3:g}, a2={4:g}, a3={5:g}, H0={6}, S0={7}, B={8}'.format(Cp0, CpInf, self.a0, self.a1, self.a2, self.a3, H0, S0, B)
        if self._Tmin != 0.0: string += ', Tmin=({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        if self._Tmax != 0.0: string += ', Tmax=({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a Wilhoit object.
        """
        return (Wilhoit, (self.Cp0, self.CpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.B, self.Tmin, self.Tmax, self.comment))

    property Cp0:
        """The (constant-pressure) heat capacity at zero temperature."""
        def __get__(self):
            return pq.Quantity(self._Cp0, pq.J / pq.mol / pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._Cp0 = 0.0 
            else:
                self._Cp0 = float(units.convertHeatCapacity(value, "J/(mol*K)"))

    property CpInf:
        """The (constant-pressure) heat capacity at infinite temperature."""
        def __get__(self):
            return pq.Quantity(self._CpInf, pq.J / pq.mol / pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._CpInf = 0.0 
            else:
                self._CpInf = float(units.convertHeatCapacity(value, "J/(mol*K)"))

    property B:
        """The scaled temperature coefficient."""
        def __get__(self):
            return pq.Quantity(self._B, pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._B = 0.0 
            else:
                self._B = float(units.convertTemperature(value, "K"))

    property H0:
        """The enthalpy integration constant."""
        def __get__(self):
            return pq.Quantity(self._H0 * 0.001, "kJ/mol")
        def __set__(self, value):
            if value is None or value == 0:
                self._H0 = 0.0 
            else:
                self._H0 = float(units.convertEnthalpy(value, "J/mol"))

    property S0:
        """The entropy integration constant."""
        def __get__(self):
            return pq.Quantity(self._S0, pq.J / pq.mol / pq.K)
        def __set__(self, value):
            if value is None or value == 0:
                self._S0 = 0.0 
            else:
                self._S0 = float(units.convertEntropy(value, "J/(mol*K)"))

    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperature `T` in K.
        """
        cdef double y
        y = T / (T + self._B)
        return self._Cp0 + (self._CpInf - self._Cp0) * y * y * (
            1 + (y - 1) * (self.a0 + y * (self.a1 + y * (self.a2 + y * self.a3))) 
        )
            
    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in kJ/mol at the specified temperature `T` in K.
        """
        cdef double y, H
        y = T / (T + self._B)
        H = self._H0 + self._Cp0 * T - (self._CpInf-self._Cp0) * T * (
            y * y * ((3 * self.a0 + self.a1 + self.a2 + self.a3) / 6. + 
                (4 * self.a1 + self.a2 + self.a3) * y / 12. + 
                (5 * self.a2 + self.a3) * y * y / 20. + 
                self.a3 * y * y * y / 5.) + 
            (2 + self.a0 + self.a1 + self.a2 + self.a3) * (y / 2. - 1 + (1.0 / y - 1.) * log(self._B + T))
        )
        return 0.001 * H
    
    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        cdef double y, logT, logy
        y = T / (T + self._B)
        logT = log(T)
        logy = log(y)
        return self._S0 + self._CpInf * logT - (self._CpInf - self._Cp0) * (
            logy + y * (1 + y * (self.a0 / 2. + y * (self.a1 / 3. + y * (self.a2 / 4. + y * self.a3 / 5.))))
        )
    
    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in kJ/mol at the specified temperature
        `T` in K.
        """
        return self.getEnthalpy(T) - 0.001 * T * self.getEntropy(T)
    
    def __residual(self, B, Tdata, Cpdata, linear, Nfreq, Nrotors, H298, S298):
        # The residual corresponding to the fitToData() method
        # Parameters are the same as for that method
        cdef double res = 0.0, diff
        cdef int i
        self.fitToDataForConstantB(Tdata, Cpdata, linear, Nfreq, Nrotors, H298, S298, B)
        # Objective function is linear least-squares
        for i in range(Cpdata.shape[0]):
            diff = self.getHeatCapacity(Tdata[i]) - Cpdata[i]
            res += diff * diff
        return res
    
    def fitToData(self, 
                  numpy.ndarray[numpy.float64_t, ndim=1] Tdata, 
                  numpy.ndarray[numpy.float64_t, ndim=1] Cpdata, 
                  bint linear, int Nfreq, int Nrotors, 
                  double H298, double S298, double B0=500.0):
        """
        Fit a Wilhoit model to the data points provided, allowing the 
        characteristic temperature `B` to vary so as to improve the fit. This
        procedure requires an optimization, using the ``fminbound`` function
        in the ``scipy.optimize`` module. The data consists of a set
        of heat capacity points `Cpdata` in J/mol*K at a given set of 
        temperatures `Tdata` in K, along with the enthalpy `H298` in kJ/mol and
        entropy `S298` in J/mol*K at 298 K. The linearity of the molecule, 
        number of vibrational frequencies, and number of internal rotors 
        (`linear`, `Nfreq`, and `Nrotors`, respectively) is used to set the 
        limits at zero and infinite temperature.
        """
        self._B = B0
        import scipy.optimize
        scipy.optimize.fminbound(self.__residual, 300.0, 3000.0, args=(Tdata, Cpdata, linear, Nfreq, Nrotors, H298, S298))
        return self
    
    def fitToDataForConstantB(self, 
                              numpy.ndarray[numpy.float64_t, ndim=1] Tdata, 
                              numpy.ndarray[numpy.float64_t, ndim=1] Cpdata, 
                              bint linear, int Nfreq, int Nrotors,
                              double H298, double S298, double B):
        """
        Fit a Wilhoit model to the data points provided using a specified value
        of the characteristic temperature `B`. The data consists of a set
        of heat capacity points `Cpdata` in J/mol*K at a given set of 
        temperatures `Tdata` in K, along with the enthalpy `H298` in J/mol and
        entropy `S298` in J/mol*K at 298 K. The linearity of the molecule, 
        number of vibrational frequencies, and number of internal rotors 
        (`linear`, `Nfreq`, and `Nrotors`, respectively) is used to set the 
        limits at zero and infinite temperature.
        """
        cdef numpy.ndarray[numpy.float64_t, ndim=1] b, x
        cdef numpy.ndarray[numpy.float64_t, ndim=2] A
        cdef double y
        cdef int i, j
        
        if Nfreq == 0:
            # Monatomic species
            assert Nrotors == 0
            self._Cp0 = 2.5 * constants.R
            self._CpInf = 2.5 * constants.R
            self._B = B
            self.a0 = 0.0
            self.a1 = 0.0
            self.a2 = 0.0
            self.a3 = 0.0
    
        else:
            # Polyatomic species
    
            # Set the Cp(T) limits as T -> and T -> infinity
            self._Cp0 = 3.5 * constants.R if linear else 4.0 * constants.R
            self._CpInf = self._Cp0 + (Nfreq + 0.5 * Nrotors) * constants.R
            
            # What remains is to fit the polynomial coefficients (a0, a1, a2, a3)
            # This can be done directly - no iteration required
            A = numpy.empty((Cpdata.shape[0],4), numpy.float64)
            b = numpy.empty(Cpdata.shape[0], numpy.float64)
            for i in range(Cpdata.shape[0]):
                y = Tdata[i] / (Tdata[i] + B)
                for j in range(4):
                    A[i,j] = (y*y*y - y*y) * y**j
                b[i] = ((Cpdata[i] - self._Cp0) / (self._CpInf - self._Cp0) - y*y)
            x, residues, rank, s = numpy.linalg.lstsq(A, b)
            
            self._B = float(B)
            self.a0 = float(x[0])
            self.a1 = float(x[1])
            self.a2 = float(x[2])
            self.a3 = float(x[3])

        self._H0 = 0.0; self._S0 = 0.0
        self._H0 = (H298 - self.getEnthalpy(298.15)) * 1000.
        self._S0 = S298 - self.getEntropy(298.15)

        return self

    cdef double integral_T0(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} \\ dT'
        
        evaluated at the given temperature `T` in kK. The implementation
        differs from that given in the Yelvington thesis for enthalpy by a 
        parameter-dependent (but temperature-independent) constant.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double y, y2, logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0, self._CpInf, self._B, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        y2 = y * y
        logBplusT = log(B + T)
        result = Cp0*T - (CpInf-Cp0)*T*(y2*((3*a0 + a1 + a2 + a3)/6. + (4*a1 + a2 + a3)*y/12. + (5*a2 + a3)*y2/20. + a3*y2*y/5.) + (2 + a0 + a1 + a2 + a3)*( y/2. - 1 + (1/y-1)*logBplusT))
        return result
    
    cdef double integral_TM1(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} (T')^{-1} \\ dT'
        
        evaluated at the given temperature `T` in kK. The implementation
        differs from that given in the Yelvington thesis for entropy by a 
        parameter-dependent (but temperature-independent) constant.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double y, logy, logT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0, self._CpInf, self._B, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        logy = log(y)
        logT = log(T)
        result = CpInf*logT-(CpInf-Cp0)*(logy+y*(1+y*(a0/2+y*(a1/3 + y*(a2/4 + y*a3/5)))))
        return result
    
    cdef double integral_T1(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} T' \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0, self._CpInf, self._B, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = ( (2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*T + (CpInf*T**2)/2. + (a3*B**7*(-Cp0 + CpInf))/(5.*(B + T)**5) + ((a2 + 6*a3)*B**6*(Cp0 - CpInf))/(4.*(B + T)**4) -
            ((a1 + 5*(a2 + 3*a3))*B**5*(Cp0 - CpInf))/(3.*(B + T)**3) + ((a0 + 4*a1 + 10*(a2 + 2*a3))*B**4*(Cp0 - CpInf))/(2.*(B + T)**2) -
            ((1 + 3*a0 + 6*a1 + 10*a2 + 15*a3)*B**3*(Cp0 - CpInf))/(B + T) - (3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(Cp0 - CpInf)*logBplusT)
        return result
    
    cdef double integral_T2(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} (T')^2 \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0, self._CpInf, self._B, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = ( -((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(Cp0 - CpInf)*T) + ((2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*T**2)/2. + (CpInf*T**3)/3. + (a3*B**8*(Cp0 - CpInf))/(5.*(B + T)**5) -
            ((a2 + 7*a3)*B**7*(Cp0 - CpInf))/(4.*(B + T)**4) + ((a1 + 6*a2 + 21*a3)*B**6*(Cp0 - CpInf))/(3.*(B + T)**3) - ((a0 + 5*(a1 + 3*a2 + 7*a3))*B**5*(Cp0 - CpInf))/(2.*(B + T)**2) +
            ((1 + 4*a0 + 10*a1 + 20*a2 + 35*a3)*B**4*(Cp0 - CpInf))/(B + T) + (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(Cp0 - CpInf)*logBplusT)
        return result
    
    cdef double integral_T3(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} (T')^3 \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0, self._CpInf, self._B, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = ( (4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(Cp0 - CpInf)*T + ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-Cp0 + CpInf)*T**2)/2. + ((2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*T**3)/3. +
            (CpInf*T**4)/4. + (a3*B**9*(-Cp0 + CpInf))/(5.*(B + T)**5) + ((a2 + 8*a3)*B**8*(Cp0 - CpInf))/(4.*(B + T)**4) - ((a1 + 7*(a2 + 4*a3))*B**7*(Cp0 - CpInf))/(3.*(B + T)**3) +
            ((a0 + 6*a1 + 21*a2 + 56*a3)*B**6*(Cp0 - CpInf))/(2.*(B + T)**2) - ((1 + 5*a0 + 15*a1 + 35*a2 + 70*a3)*B**5*(Cp0 - CpInf))/(B + T) -
            (5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(Cp0 - CpInf)*logBplusT)
        return result
    
    cdef double integral_T4(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} (T')^4 \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0, self._CpInf, self._B, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = ( -((5 + 10*a0 + 20*a1 + 35*a2 + 56*a3)*B**4*(Cp0 - CpInf)*T) + ((4 + 6*a0 + 10*a1 + 15*a2 + 21*a3)*B**3*(Cp0 - CpInf)*T**2)/2. +
            ((3 + 3*a0 + 4*a1 + 5*a2 + 6*a3)*B**2*(-Cp0 + CpInf)*T**3)/3. + ((2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*T**4)/4. + (CpInf*T**5)/5. + (a3*B**10*(Cp0 - CpInf))/(5.*(B + T)**5) -
            ((a2 + 9*a3)*B**9*(Cp0 - CpInf))/(4.*(B + T)**4) + ((a1 + 8*a2 + 36*a3)*B**8*(Cp0 - CpInf))/(3.*(B + T)**3) - ((a0 + 7*(a1 + 4*(a2 + 3*a3)))*B**7*(Cp0 - CpInf))/(2.*(B + T)**2) +
            ((1 + 6*a0 + 21*a1 + 56*a2 + 126*a3)*B**6*(Cp0 - CpInf))/(B + T) + (6 + 15*a0 + 35*a1 + 70*a2 + 126*a3)*B**5*(Cp0 - CpInf)*logBplusT)
        return result
    
    cdef double integral2_T0(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\left[ \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} \\right]^2 \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0, self._CpInf, self._B, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T)
        result = (CpInf**2*T - (a3**2*B**12*(Cp0 - CpInf)**2)/(11.*(B + T)**11) + (a3*(a2 + 5*a3)*B**11*(Cp0 - CpInf)**2)/(5.*(B + T)**10) -
            ((a2**2 + 18*a2*a3 + a3*(2*a1 + 45*a3))*B**10*(Cp0 - CpInf)**2)/(9.*(B + T)**9) + ((4*a2**2 + 36*a2*a3 + a1*(a2 + 8*a3) + a3*(a0 + 60*a3))*B**9*(Cp0 - CpInf)**2)/(4.*(B + T)**8) -
            ((a1**2 + 14*a1*(a2 + 4*a3) + 2*(14*a2**2 + a3 + 84*a2*a3 + 105*a3**2 + a0*(a2 + 7*a3)))*B**8*(Cp0 - CpInf)**2)/(7.*(B + T)**7) +
            ((3*a1**2 + a2 + 28*a2**2 + 7*a3 + 126*a2*a3 + 126*a3**2 + 7*a1*(3*a2 + 8*a3) + a0*(a1 + 6*a2 + 21*a3))*B**7*(Cp0 - CpInf)**2)/(3.*(B + T)**6) -
            (B**6*(Cp0 - CpInf)*(a0**2*(Cp0 - CpInf) + 15*a1**2*(Cp0 - CpInf) + 10*a0*(a1 + 3*a2 + 7*a3)*(Cp0 - CpInf) + 2*a1*(1 + 35*a2 + 70*a3)*(Cp0 - CpInf) +
             2*(35*a2**2*(Cp0 - CpInf) + 6*a2*(1 + 21*a3)*(Cp0 - CpInf) + a3*(5*(4 + 21*a3)*Cp0 - 21*(CpInf + 5*a3*CpInf)))))/(5.*(B + T)**5) +
            (B**5*(Cp0 - CpInf)*(14*a2*Cp0 + 28*a2**2*Cp0 + 30*a3*Cp0 + 84*a2*a3*Cp0 + 60*a3**2*Cp0 + 2*a0**2*(Cp0 - CpInf) + 10*a1**2*(Cp0 - CpInf) +
             a0*(1 + 10*a1 + 20*a2 + 35*a3)*(Cp0 - CpInf) + a1*(5 + 35*a2 + 56*a3)*(Cp0 - CpInf) - 15*a2*CpInf - 28*a2**2*CpInf - 35*a3*CpInf - 84*a2*a3*CpInf - 60*a3**2*CpInf))/
             (2.*(B + T)**4) - (B**4*(Cp0 - CpInf)*((1 + 6*a0**2 + 15*a1**2 + 32*a2 + 28*a2**2 + 50*a3 + 72*a2*a3 + 45*a3**2 + 2*a1*(9 + 21*a2 + 28*a3) + a0*(8 + 20*a1 + 30*a2 + 42*a3))*Cp0 -
             (1 + 6*a0**2 + 15*a1**2 + 40*a2 + 28*a2**2 + 70*a3 + 72*a2*a3 + 45*a3**2 + a0*(8 + 20*a1 + 30*a2 + 42*a3) + a1*(20 + 42*a2 + 56*a3))*CpInf))/(3.*(B + T)**3) +
            (B**3*(Cp0 - CpInf)*((2 + 2*a0**2 + 3*a1**2 + 9*a2 + 4*a2**2 + 11*a3 + 9*a2*a3 + 5*a3**2 + a0*(5 + 5*a1 + 6*a2 + 7*a3) + a1*(7 + 7*a2 + 8*a3))*Cp0 -
             (2 + 2*a0**2 + 3*a1**2 + 15*a2 + 4*a2**2 + 21*a3 + 9*a2*a3 + 5*a3**2 + a0*(6 + 5*a1 + 6*a2 + 7*a3) + a1*(10 + 7*a2 + 8*a3))*CpInf))/(B + T)**2 -
            (B**2*((2 + a0 + a1 + a2 + a3)**2*Cp0**2 - 2*(5 + a0**2 + a1**2 + 8*a2 + a2**2 + 9*a3 + 2*a2*a3 + a3**2 + 2*a0*(3 + a1 + a2 + a3) + a1*(7 + 2*a2 + 2*a3))*Cp0*CpInf +
             (6 + a0**2 + a1**2 + 12*a2 + a2**2 + 14*a3 + 2*a2*a3 + a3**2 + 2*a1*(5 + a2 + a3) + 2*a0*(4 + a1 + a2 + a3))*CpInf**2))/(B + T) +
            2*(2 + a0 + a1 + a2 + a3)*B*(Cp0 - CpInf)*CpInf*logBplusT)
        return result
    
    cdef double integral2_TM1(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\left[ \\frac{C_\\mathrm{p}^\\mathrm{Wilhoit}(T')}{R} \\right]^2 (T')^{-1} \\ dT'
        
        evaluated at the given temperature `T` in kK.
        """
        cdef double Cp0, CpInf, B, a0, a1, a2, a3
        cdef double logBplusT, result
        Cp0, CpInf, B, a0, a1, a2, a3 = self._Cp0, self._CpInf, self._B, self.a0, self.a1, self.a2, self.a3
        logBplusT = log(B + T); logT = log(T)
        result = ( (a3**2*B**11*(Cp0 - CpInf)**2)/(11.*(B + T)**11) - (a3*(2*a2 + 9*a3)*B**10*(Cp0 - CpInf)**2)/(10.*(B + T)**10) +
            ((a2**2 + 16*a2*a3 + 2*a3*(a1 + 18*a3))*B**9*(Cp0 - CpInf)**2)/(9.*(B + T)**9) -
            ((7*a2**2 + 56*a2*a3 + 2*a1*(a2 + 7*a3) + 2*a3*(a0 + 42*a3))*B**8*(Cp0 - CpInf)**2)/(8.*(B + T)**8) +
            ((a1**2 + 21*a2**2 + 2*a3 + 112*a2*a3 + 126*a3**2 + 2*a0*(a2 + 6*a3) + 6*a1*(2*a2 + 7*a3))*B**7*(Cp0 - CpInf)**2)/(7.*(B + T)**7) -
            ((5*a1**2 + 2*a2 + 30*a1*a2 + 35*a2**2 + 12*a3 + 70*a1*a3 + 140*a2*a3 + 126*a3**2 + 2*a0*(a1 + 5*(a2 + 3*a3)))*B**6*(Cp0 - CpInf)**2)/(6.*(B + T)**6) +
            (B**5*(Cp0 - CpInf)*(10*a2*Cp0 + 35*a2**2*Cp0 + 28*a3*Cp0 + 112*a2*a3*Cp0 + 84*a3**2*Cp0 + a0**2*(Cp0 - CpInf) + 10*a1**2*(Cp0 - CpInf) + 2*a1*(1 + 20*a2 + 35*a3)*(Cp0 - CpInf) +
            4*a0*(2*a1 + 5*(a2 + 2*a3))*(Cp0 - CpInf) - 10*a2*CpInf - 35*a2**2*CpInf - 30*a3*CpInf - 112*a2*a3*CpInf - 84*a3**2*CpInf))/(5.*(B + T)**5) -
            (B**4*(Cp0 - CpInf)*(18*a2*Cp0 + 21*a2**2*Cp0 + 32*a3*Cp0 + 56*a2*a3*Cp0 + 36*a3**2*Cp0 + 3*a0**2*(Cp0 - CpInf) + 10*a1**2*(Cp0 - CpInf) +
            2*a0*(1 + 6*a1 + 10*a2 + 15*a3)*(Cp0 - CpInf) + 2*a1*(4 + 15*a2 + 21*a3)*(Cp0 - CpInf) - 20*a2*CpInf - 21*a2**2*CpInf - 40*a3*CpInf - 56*a2*a3*CpInf - 36*a3**2*CpInf))/
            (4.*(B + T)**4) + (B**3*(Cp0 - CpInf)*((1 + 3*a0**2 + 5*a1**2 + 14*a2 + 7*a2**2 + 18*a3 + 16*a2*a3 + 9*a3**2 + 2*a0*(3 + 4*a1 + 5*a2 + 6*a3) + 2*a1*(5 + 6*a2 + 7*a3))*Cp0 -
            (1 + 3*a0**2 + 5*a1**2 + 20*a2 + 7*a2**2 + 30*a3 + 16*a2*a3 + 9*a3**2 + 2*a0*(3 + 4*a1 + 5*a2 + 6*a3) + 2*a1*(6 + 6*a2 + 7*a3))*CpInf))/(3.*(B + T)**3) -
            (B**2*((3 + a0**2 + a1**2 + 4*a2 + a2**2 + 4*a3 + 2*a2*a3 + a3**2 + 2*a1*(2 + a2 + a3) + 2*a0*(2 + a1 + a2 + a3))*Cp0**2 -
            2*(3 + a0**2 + a1**2 + 7*a2 + a2**2 + 8*a3 + 2*a2*a3 + a3**2 + 2*a1*(3 + a2 + a3) + a0*(5 + 2*a1 + 2*a2 + 2*a3))*Cp0*CpInf +
            (3 + a0**2 + a1**2 + 10*a2 + a2**2 + 12*a3 + 2*a2*a3 + a3**2 + 2*a1*(4 + a2 + a3) + 2*a0*(3 + a1 + a2 + a3))*CpInf**2))/(2.*(B + T)**2) +
            (B*(Cp0 - CpInf)*(Cp0 - (3 + 2*a0 + 2*a1 + 2*a2 + 2*a3)*CpInf))/(B + T) + Cp0**2*logT + (-Cp0**2 + CpInf**2)*logBplusT)
        return result

################################################################################

cdef class NASA(HeatCapacityModel):
    """
    A thermodynamics model based on the NASA polynomial. Both the 
    seven-coefficient and nine-coefficient variations are supported.
    """
    
    def __init__(self, coeffs=None, Tmin=None, Tmax=None, comment=''):
        HeatCapacityModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.coeffs = coeffs
        
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'NASA('
        if self.cm2 == 0 and self.cm1 == 0:
            string += 'coeffs=[{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g}]'.format(self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
        else:
            string += 'coeffs=[{0:g},{1:g},{2:g},{3:g},{4:g},{5:g},{6:g},{7:g},{8:g}]'.format(self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6)
        string += ', Tmin=({0:g},"{1}")'.format(float(self.Tmin), str(self.Tmin.dimensionality))
        string += ', Tmax=({0:g},"{1}")'.format(float(self.Tmax), str(self.Tmax.dimensionality))
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (NASA, ([self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6], self.Tmin, self.Tmax, self.comment))

    property coeffs:
        """The set of seven or nine NASA polynomial coefficients."""
        def __get__(self):
            if self.cm2 == 0 and self.cm1 == 0:
                return numpy.array([self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6])
            else:
                return numpy.array([self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6])
        def __set__(self, value):
            if value is None:
                value = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            if len(value) == 7:
                self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = value
                self.cm2 = 0; self.cm1 = 0
            elif len(value) == 9:
                self.cm2, self.cm1, self.c0, self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = value
            else:
                raise ValueError('Invalid number of NASA polynomial coefficients; should be 7 or 9.')
    
    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperature `T` in K.
        """
        return ((self.cm2 / T + self.cm1) / T + self.c0 + T*(self.c1 + T*(self.c2 + T*(self.c3 + self.c4*T)))) * constants.R
    
    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in kJ/mol at the specified temperature `T` in K.
        """
        cdef double T2 = T * T
        cdef double T4 = T2 * T2
        cdef double H
        H = ((-self.cm2 / T + self.cm1 * log(T)) / T + self.c0 + self.c1*T/2. + self.c2*T2/3. + self.c3*T2*T/4. + self.c4*T4/5. + self.c5/T) * constants.R * T
        return 0.001 * H
    
    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        cdef double T2 = T * T
        cdef double T4 = T2 * T2
        return ((-self.cm2 / T / 2. - self.cm1) / T + self.c0*log(T) + self.c1*T + self.c2*T2/2. +
            self.c3*T2*T/3. + self.c4*T4/4. + self.c6 ) * constants.R
    
    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in kJ/mol at the specified temperature
        `T` in K.
        """
        return self.getEnthalpy(T) - 0.001 * T * self.getEntropy(T)
    
    cdef double integral2_T0(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\left[ \\frac{C_\\mathrm{p}^\\mathrm{NASA}(T')}{R} \\right]^2 \\ dT'
        
        evaluated at the given temperature `T` in K.
        """
        cdef double c0, c1, c2, c3, c4
        cdef double T2, T4, T8, result
        c0, c1, c2, c3, c4 = self.c0, self.c1, self.c2, self.c3, self.c4
        T2 = T * T
        T4 = T2 * T2
        T8 = T4 * T4
        result = (
            c0*c0*T + c0*c1*T2 + 2./3.*c0*c2*T2*T + 0.5*c0*c3*T4 + 0.4*c0*c4*T4*T +
            c1*c1*T2*T/3. + 0.5*c1*c2*T4 + 0.4*c1*c3*T4*T + c1*c4*T4*T2/3. +
            0.2*c2*c2*T4*T + c2*c3*T4*T2/3. + 2./7.*c2*c4*T4*T2*T +
            c3*c3*T4*T2*T/7. + 0.25*c3*c4*T8 +
            c4*c4*T8*T/9.
        )
        return result
    
    cdef double integral2_TM1(self, double T):
        """
        Return the value of the dimensionless integral
        
        .. math:: \\int \\left[ \\frac{C_\\mathrm{p}^\\mathrm{NASA}(T')}{R} \\right]^2 (T')^{-1} \\ dT'
        
        evaluated at the given temperature `T` in K.
        """
        cdef double c0, c1, c2, c3, c4
        cdef double T2, T4, logT, result
        c0, c1, c2, c3, c4 = self.c0, self.c1, self.c2, self.c3, self.c4
        T2 = T * T
        T4 = T2 * T2
        logT = log(T)
        result = (
            c0*c0*logT + 2*c0*c1*T + c0*c2*T2 + 2./3.*c0*c3*T2*T + 0.5*c0*c4*T4 +
            0.5*c1*c1*T2 + 2./3.*c1*c2*T2*T + 0.5*c1*c3*T4 + 0.4*c1*c4*T4*T +
            0.25*c2*c2*T4 + 0.4*c2*c3*T4*T + c2*c4*T4*T2/3. +
            c3*c3*T4*T2/6. + 2./7.*c3*c4*T4*T2*T +
            c4*c4*T4*T4/8.
        )
        return result

################################################################################

cdef class MultiNASA(HeatCapacityModel):
    """
    A set of thermodynamic parameters given by NASA polynomials. This class
    stores a list of :class:`NASA` objects in the `polynomials`
    attribute. When evaluating a thermodynamic quantity, a polynomial that
    contains the desired temperature within its valid range will be used.
    """
    
    def __init__(self, polynomials=None, Tmin=0.0, Tmax=0.0, comment=''):
        HeatCapacityModel.__init__(self, Tmin=Tmin, Tmax=Tmax, comment=comment)
        self.polynomials = polynomials or []
    
    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        MultiNASA object.
        """
        string = 'MultiNASA('
        string += 'Tmin=({0:g},"{1}"), '.format(float(self.Tmin), str(self.Tmin.dimensionality))
        string += 'Tmax=({0:g},"{1}"), '.format(float(self.Tmax), str(self.Tmax.dimensionality))
        string += 'polynomials=[{0}]'.format(','.join([repr(poly) for poly in self.polynomials]))
        if self.comment != '': string += ', comment="""{0}"""'.format(self.comment)
        string += ')'
        return string

    def __reduce__(self):
        """
        A helper function used when pickling a MultiNASA object.
        """
        return (MultiNASA, (self.polynomials, self.Tmin, self.Tmax, self.comment))

    cpdef double getHeatCapacity(self, double T) except -1000000000:
        """
        Return the constant-pressure heat capacity in J/mol*K at the
        specified temperature `T` in K.
        """
        return self.getPolynomialForTemperature(T).getHeatCapacity(T)
    
    cpdef double getEnthalpy(self, double T) except 1000000000:
        """
        Return the enthalpy in kJ/mol at the specified temperature `T` in K.
        """
        return self.getPolynomialForTemperature(T).getEnthalpy(T)
    
    cpdef double getEntropy(self, double T) except -1000000000:
        """
        Return the entropy in J/mol*K at the specified temperature `T` in K.
        """
        return self.getPolynomialForTemperature(T).getEntropy(T)
    
    cpdef double getFreeEnergy(self, double T) except 1000000000:
        """
        Return the Gibbs free energy in kJ/mol at the specified temperature
        `T` in K.
        """
        return self.getPolynomialForTemperature(T).getFreeEnergy(T)
    
    cpdef NASA getPolynomialForTemperature(self, double T):
        cdef NASA poly
        for poly in self.polynomials:
            if poly.isTemperatureValid(T): return poly
        else:
            raise ValueError("No valid NASA polynomial found for T={0:g} K".format(T))
