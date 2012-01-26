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

################################################################################

cpdef getRateCoefficientUnitsFromReactionOrder(order)

cpdef getReactionOrderFromRateCoefficientUnits(kunits)

################################################################################

cdef class KineticsModel:
    
    cdef double _Tmin, _Tmax
    cdef double _Pmin, _Pmax
    cdef public int order
    cdef public str comment
    
    cpdef bint isTemperatureValid(self, double T) except -2

    cpdef bint isPressureValid(self, double P) except -2

    cpdef bint isPressureDependent(self) except -2

    cpdef double getRateCoefficient(self, double T, double P=?) except -1

################################################################################

cdef class Arrhenius(KineticsModel):
    
    cdef double _A, _Ea, _T0
    cdef public double n
    
    cpdef bint isPressureDependent(self) except -2

    cpdef double getRateCoefficient(self, double T, double P=?) except -1

    cpdef changeT0(self, double T0)

    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray klist, kunits, double T0=?, numpy.ndarray weights=?, bint threeParams=?)

################################################################################

cdef class PDepArrhenius(KineticsModel):
    
    cdef numpy.ndarray _pressures
    cdef public list arrhenius
    
    cpdef bint isPressureDependent(self) except -2

    cdef getAdjacentExpressions(self, double P)
    
    cpdef double getRateCoefficient(self, double T, double P=?) except -1
    
    cpdef fitToData(self, numpy.ndarray Tlist, numpy.ndarray Plist, numpy.ndarray K, kunits, double T0=?)
