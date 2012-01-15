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

cimport numpy

cdef class HeatCapacityModel:
    
    cdef double _Tmin, _Tmax
    cdef public str comment
    
    cpdef bint isTemperatureValid(self, double T) except -2

    cpdef double getHeatCapacity(self, double T) except -1000000000

    cpdef double getEnthalpy(self, double T) except 1000000000

    cpdef double getEntropy(self, double T) except -1000000000

    cpdef double getFreeEnergy(self, double T) except 1000000000

################################################################################

cdef class Wilhoit(HeatCapacityModel):
    
    cdef public _Cp0, _CpInf, _B, _H0, _S0
    cdef public double a0, a1, a2, a3
    
    cpdef double getHeatCapacity(self, double T) except -1000000000

    cpdef double getEnthalpy(self, double T) except 1000000000

    cpdef double getEntropy(self, double T) except -1000000000

    cpdef double getFreeEnergy(self, double T) except 1000000000

################################################################################

cdef class NASA(HeatCapacityModel):
    
    cdef public double cm2, cm1, c0, c1, c2, c3, c4, c5, c6
    
    cpdef double getHeatCapacity(self, double T) except -1000000000

    cpdef double getEnthalpy(self, double T) except 1000000000

    cpdef double getEntropy(self, double T) except -1000000000

    cpdef double getFreeEnergy(self, double T) except 1000000000

################################################################################

cdef class MultiNASA(HeatCapacityModel):

    cdef public list polynomials

    cpdef double getHeatCapacity(self, double T) except -1000000000

    cpdef double getEnthalpy(self, double T) except 1000000000

    cpdef double getEntropy(self, double T) except -1000000000

    cpdef double getFreeEnergy(self, double T) except 1000000000

    cpdef NASA getPolynomialForTemperature(self, double T)
