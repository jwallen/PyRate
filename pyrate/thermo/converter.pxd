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
cimport cython

from pyrate.thermo.model cimport Wilhoit, NASA, MultiNASA
cimport pyrate.constants as constants

cpdef MultiNASA convertWilhoitToMultiNASA(Wilhoit wilhoit, double Tmin, double Tmax, double Tint, bint fixedTint=?, bint weighting=?, int continuity=?)

cpdef Wilhoit convertMultiNASAToWilhoit(MultiNASA multiNASA, bint linear, int Nfreq, int Nrotors)

################################################################################

cpdef Wilhoit_to_NASA(Wilhoit wilhoit, double Tmin, double Tmax, double Tint, bint weighting, int contCons)

cpdef Wilhoit_to_NASA_TintOpt(Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int contCons)

cpdef double Wilhoit_to_NASA_TintOpt_objFun(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int contCons)

cpdef double Wilhoit_to_NASA_TintOpt_objFun_NW(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int contCons)

cpdef double Wilhoit_to_NASA_TintOpt_objFun_W(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int contCons)
