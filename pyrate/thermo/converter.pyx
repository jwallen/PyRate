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
This module contains functions for converting between various heat capacity
models.
"""

import cython
import numpy
import logging
import quantities as pq

from libc.math cimport sqrt, log

from pyrate.thermo.model cimport Wilhoit, NASA, MultiNASA
import pyrate.constants as constants

################################################################################

cpdef MultiNASA convertWilhoitToMultiNASA(Wilhoit wilhoit, double Tmin, double Tmax, double Tint, bint fixedTint=False, bint weighting=True, int continuity=3):
    """
    Convert a :class:`Wilhoit` object `wilhoit` to a :class:`MultiNASA` 
    object. You must specify the minimum and maximum temperatures of the fit
    `Tmin` and `Tmax`, as well as the intermediate temperature `Tint` to use
    as the bridge between the two fitted polynomials. The remaining parameters
    can be used to modify the fitting algorithm used:
    
    * `fixedTint` - ``False`` to allow `Tint` to vary in order to improve the fit, or ``True`` to keep it fixed

    * `weighting` - ``True`` to weight the fit by :math:`T^{-1}` to emphasize good fit at lower temperatures, or ``False`` to not use weighting

    * `continuity` - The number of continuity constraints to enforce at `Tint`:

        - 0: no constraints on continuity of :math:`C_\\mathrm{p}(T)` at `Tint`

        - 1: constrain :math:`C_\\mathrm{p}(T)` to be continous at `Tint`

        - 2: constrain :math:`C_\\mathrm{p}(T)` and :math:`\\frac{d C_\\mathrm{p}}{dT}` to be continuous at `Tint`

        - 3: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, and :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}` to be continuous at `Tint`

        - 4: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}`, and :math:`\\frac{d^3 C_\\mathrm{p}}{dT^3}` to be continuous at `Tint`

        - 5: constrain :math:`C_\\mathrm{p}(T)`, :math:`\\frac{d C_\\mathrm{p}}{dT}`, :math:`\\frac{d^2 C_\\mathrm{p}}{dT^2}`, :math:`\\frac{d^3 C_\\mathrm{p}}{dT^3}`, and :math:`\\frac{d^4 C_\\mathrm{p}}{dT^4}` to be continuous at `Tint`
        
    Note that values of `continuity` of 5 or higher effectively constrain all
    the coefficients to be equal and should be equivalent to fitting only one
    polynomial (rather than two).

    Returns the fitted :class:`MultiNASA` object containing the two fitted
    :class:`NASA` objects.
    """
    cdef Wilhoit wilhoit_scaled
    cdef NASA nasa_low, nasa_high
    cdef double iseUnw, rmsUnw, iseWei, rmsWei, Hlow, Slow, Hhigh, Shigh
    cdef str rmsStr
    
    # Scale the temperatures to kK
    Tmin /= 1000.
    Tint /= 1000.
    Tmax /= 1000.
    
    # Make copy of Wilhoit data so we don't modify the original
    wilhoit_scaled = Wilhoit(wilhoit.Cp0, wilhoit.CpInf, wilhoit.a0, wilhoit.a1, wilhoit.a2, wilhoit.a3, wilhoit.H0, wilhoit.S0, wilhoit.B, Tmin=wilhoit.Tmin, Tmax=wilhoit.Tmax, comment=wilhoit.comment)
    # Rescale Wilhoit parameters
    wilhoit_scaled.Cp0 /= constants.R
    wilhoit_scaled.CpInf /= constants.R
    wilhoit_scaled.B /= 1000.
    
    #if we are using fixed Tint, do not allow Tint to float
    if fixedTint:
        nasa_low, nasa_high = Wilhoit_to_NASA(wilhoit_scaled, Tmin, Tmax, Tint, weighting, continuity)
    else:
        nasa_low, nasa_high, Tint = Wilhoit_to_NASA_TintOpt(wilhoit_scaled, Tmin, Tmax, weighting, continuity)
    iseUnw = Wilhoit_to_NASA_TintOpt_objFun(Tint, wilhoit_scaled, Tmin, Tmax, 0, continuity) #the scaled, unweighted ISE (integral of squared error)
    rmsUnw = sqrt(iseUnw/(Tmax-Tmin))
    rmsStr = '(Unweighted) RMS error = %.3f*R;'%(rmsUnw)
    if (weighting == 1):
        iseWei = Wilhoit_to_NASA_TintOpt_objFun(Tint, wilhoit_scaled, Tmin, Tmax, weighting, continuity) #the scaled, weighted ISE
        rmsWei = sqrt(iseWei/log(Tmax/Tmin))
        rmsStr = 'Weighted RMS error = %.3f*R;'%(rmsWei)+rmsStr

    #print a warning if the rms fit is worse that 0.25*R
    if (rmsUnw > 0.25 or rmsWei > 0.25):
        logging.warning("Poor Wilhoit-to-NASA fit quality: RMS error = %.3f*R" % (rmsWei if weighting == 1 else rmsUnw))

    #restore to conventional units of K for Tint and units based on K rather than kK in NASA polynomial coefficients
    Tint *= 1000.
    Tmin *= 1000.
    Tmax *= 1000.

    nasa_low.c1 /= 1000.
    nasa_low.c2 /= 1000000.
    nasa_low.c3 /= 1000000000.
    nasa_low.c4 /= 1000000000000.
    
    nasa_high.c1 /= 1000.
    nasa_high.c2 /= 1000000.
    nasa_high.c3 /= 1000000000.
    nasa_high.c4 /= 1000000000000.

    # output comment
    comment = 'NASA function fitted to Wilhoit function. ' + rmsStr + wilhoit.comment
    nasa_low.Tmin = (Tmin,"K")
    nasa_low.Tmax = (Tint,"K")
    nasa_low.comment = 'Low temperature range polynomial'
    nasa_high.Tmin = (Tint,"K")
    nasa_high.Tmax = (Tmax,"K")
    nasa_high.comment = 'High temperature range polynomial'

    #for the low polynomial, we want the results to match the Wilhoit value at 298.15K
    #low polynomial enthalpy:
    Hlow = (wilhoit.getEnthalpy(298) - nasa_low.getEnthalpy(298)) * 1000. /constants.R
    #low polynomial entropy:
    Slow = (wilhoit.getEntropy(298) - nasa_low.getEntropy(298))/constants.R

    # update last two coefficients
    nasa_low.c5 = Hlow
    nasa_low.c6 = Slow

    #for the high polynomial, we want the results to match the low polynomial value at Tint
    #high polynomial enthalpy:
    Hhigh = (nasa_low.getEnthalpy(Tint) - nasa_high.getEnthalpy(Tint)) * 1000. /constants.R
    #high polynomial entropy:
    Shigh = (nasa_low.getEntropy(Tint) - nasa_high.getEntropy(Tint))/constants.R

    # update last two coefficients
    #polynomial_high.coeffs = (b6,b7,b8,b9,b10,Hhigh,Shigh)
    nasa_high.c5 = Hhigh
    nasa_high.c6 = Shigh

    return MultiNASA(Tmin=(Tmin,"K"), Tmax=(Tmax,"K"), polynomials=[nasa_low,nasa_high], comment=comment)

################################################################################

@cython.boundscheck(False)
cpdef Wilhoit convertMultiNASAToWilhoit(MultiNASA multiNASA, bint linear, int Nfreq, int Nrotors):
    """
    Convert a :class:`MultiNASA` object `multiNASA` to a :class:`Wilhoit` 
    object. You must specify the linearity of the molecule `linear`, the number
    of vibrational modes `Nfreq`, and the number of hindered rotor modes
    `Nrotors` so the algorithm can determine the appropriate heat capacity
    limits at zero and infinite temperature.
    """
    cdef double Tmin, Tmax, dT, H298, S298
    cdef numpy.ndarray[numpy.float64_t, ndim=1] Tdata, Cpdata
    cdef int i
    
    Tmin = float(multiNASA.Tmin)
    Tmax = float(multiNASA.Tmax)
    dT = min(50.0, (Tmax - Tmin) / 100.)
    
    Tdata = numpy.arange(Tmin, Tmax, dT)
    Cpdata = numpy.zeros_like(Tdata)
    
    for i in range(Tdata.shape[0]):
        Cpdata[i] = multiNASA.getHeatCapacity(Tdata[i])
    H298 = multiNASA.getEnthalpy(298)
    S298 = multiNASA.getEntropy(298)
    
    return Wilhoit().fitToData(Tdata, Cpdata, linear, Nfreq, Nrotors, H298, S298)

################################################################################

@cython.boundscheck(False)
cpdef Wilhoit_to_NASA(Wilhoit wilhoit, double Tmin, double Tmax, double Tint, bint weighting, int contCons):
    """
    Convert a Wilhoit polynomial to a pair of NASA polynomials.
    
    :param wilhoit: The Wilhoit polynomial to convert, with dimensionless heat
                    capacity limits and scaled temperature coefficient in kK
    :param Tmin:    The minimum temperature of the low-temperature NASA 
                    polynomial, in kK
    :param Tmax:    The maximum temperature of the high-temperature NASA 
                    polynomial, in kK
    :param Tint:    The intermediate temperature dividing the low-temperature
                    and high-temperature NASA polynomials, in kK
    :param weighting: ``True`` to weight the fit by inverse temperature, ``False`` to apply no weighting
    :param contCons:  The number of continuity constraints to apply to the 
                      fitted NASA polynomials at `Tint`:
        0: no constraints on continuity of Cp(T) at Tint
        1: constrain Cp to be continuous at Tint
        2: constrain Cp and dCp/dT to be continuous at Tint
        3 (default): constrain Cp, dCp/dT, and d2Cp/dT2 to be continuous at Tint
        4: constrain Cp, dCp/dT, d2Cp/dT2, and d3Cp/dT3 to be continuous at Tint
        5: constrain Cp, dCp/dT, d2Cp/dT2, d3Cp/dT3, and d4Cp/dT4 to be continuous at Tint; note: this effectively constrains all the coefficients to be equal and should be equivalent to fitting only one polynomial (rather than two)
        
        note: 5th (and higher) derivatives of NASA Cp(T) are zero and hence will automatically be continuous at Tint by the form of the Cp(T) function
                     
    :result: The pair of NASA polynomials with scaled parameters
    """
    cdef numpy.ndarray[numpy.float64_t, ndim=2] A
    cdef numpy.ndarray[numpy.float64_t, ndim=1] b, x
    cdef double w0min, w1min, w2min, w3min, w4min, wM1min
    cdef double w0int, w1int, w2int, w3int, w4int, wM1int
    cdef double w0max, w1max, w2max, w3max, w4max, wM1max
    cdef NASA nasa_low, nasa_high
    cdef int i, j
    
    #construct (typically 13*13) symmetric A matrix (in A*x = b); other elements will be zero
    A = numpy.zeros([10+contCons,10+contCons])
    b = numpy.zeros([10+contCons])

    if weighting:
        A[0,0] = 2*log(Tint/Tmin)
        A[0,1] = 2*(Tint - Tmin)
        A[0,2] = Tint*Tint - Tmin*Tmin
        A[0,3] = 2.*(Tint*Tint*Tint - Tmin*Tmin*Tmin)/3
        A[0,4] = (Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin)/2
        A[1,4] = 2.*(Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin)/5
        A[2,4] = (Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/3
        A[3,4] = 2.*(Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/7
        A[4,4] = (Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/4
    else:
        A[0,0] = 2*(Tint - Tmin)
        A[0,1] = Tint*Tint - Tmin*Tmin
        A[0,2] = 2.*(Tint*Tint*Tint - Tmin*Tmin*Tmin)/3
        A[0,3] = (Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin)/2
        A[0,4] = 2.*(Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin)/5
        A[1,4] = (Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/3
        A[2,4] = 2.*(Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/7
        A[3,4] = (Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/4
        A[4,4] = 2.*(Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint - Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin*Tmin)/9
    A[1,1] = A[0,2]
    A[1,2] = A[0,3]
    A[1,3] = A[0,4]
    A[2,2] = A[0,4]
    A[2,3] = A[1,4]
    A[3,3] = A[2,4]

    if weighting:
        A[5,5] = 2*log(Tmax/Tint)
        A[5,6] = 2*(Tmax - Tint)
        A[5,7] = Tmax*Tmax - Tint*Tint
        A[5,8] = 2.*(Tmax*Tmax*Tmax - Tint*Tint*Tint)/3
        A[5,9] = (Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint)/2
        A[6,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint)/5
        A[7,9] = (Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint)/3
        A[8,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint)/7
        A[9,9] = (Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint)/4
    else:
        A[5,5] = 2*(Tmax - Tint)
        A[5,6] = Tmax*Tmax - Tint*Tint
        A[5,7] = 2.*(Tmax*Tmax*Tmax - Tint*Tint*Tint)/3
        A[5,8] = (Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint)/2
        A[5,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint)/5
        A[6,9] = (Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint)/3
        A[7,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint)/7
        A[8,9] = (Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint)/4
        A[9,9] = 2.*(Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax*Tmax - Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint*Tint)/9
    A[6,6] = A[5,7]
    A[6,7] = A[5,8]
    A[6,8] = A[5,9]
    A[7,7] = A[5,9]
    A[7,8] = A[6,9]
    A[8,8] = A[7,9]

    if(contCons > 0):#set non-zero elements in the 11th column for Cp(T) continuity contraint
        A[0,10] = 1.
        A[1,10] = Tint
        A[2,10] = Tint*Tint
        A[3,10] = A[2,10]*Tint
        A[4,10] = A[3,10]*Tint
        A[5,10] = -A[0,10]
        A[6,10] = -A[1,10]
        A[7,10] = -A[2,10]
        A[8,10] = -A[3,10]
        A[9,10] = -A[4,10]
        if(contCons > 1): #set non-zero elements in the 12th column for dCp/dT continuity constraint
            A[1,11] = 1.
            A[2,11] = 2*Tint
            A[3,11] = 3*A[2,10]
            A[4,11] = 4*A[3,10]
            A[6,11] = -A[1,11]
            A[7,11] = -A[2,11]
            A[8,11] = -A[3,11]
            A[9,11] = -A[4,11]
            if(contCons > 2): #set non-zero elements in the 13th column for d2Cp/dT2 continuity constraint
                A[2,12] = 2.
                A[3,12] = 6*Tint
                A[4,12] = 12*A[2,10]
                A[7,12] = -A[2,12]
                A[8,12] = -A[3,12]
                A[9,12] = -A[4,12]
                if(contCons > 3): #set non-zero elements in the 14th column for d3Cp/dT3 continuity constraint
                    A[3,13] = 6
                    A[4,13] = 24*Tint
                    A[8,13] = -A[3,13]
                    A[9,13] = -A[4,13]
                    if(contCons > 4): #set non-zero elements in the 15th column for d4Cp/dT4 continuity constraint
                        A[4,14] = 24
                        A[9,14] = -A[4,14]

    # make the matrix symmetric
    for i in range(1,10+contCons):
        for j in range(0, i):
            A[i,j] = A[j,i]

    #construct b vector
    w0int = wilhoit.integral_T0(Tint)
    w1int = wilhoit.integral_T1(Tint)
    w2int = wilhoit.integral_T2(Tint)
    w3int = wilhoit.integral_T3(Tint)
    w0min = wilhoit.integral_T0(Tmin)
    w1min = wilhoit.integral_T1(Tmin)
    w2min = wilhoit.integral_T2(Tmin)
    w3min = wilhoit.integral_T3(Tmin)
    w0max = wilhoit.integral_T0(Tmax)
    w1max = wilhoit.integral_T1(Tmax)
    w2max = wilhoit.integral_T2(Tmax)
    w3max = wilhoit.integral_T3( Tmax)
    if weighting:
        wM1int = wilhoit.integral_TM1(Tint)
        wM1min = wilhoit.integral_TM1(Tmin)
        wM1max = wilhoit.integral_TM1(Tmax)
    else:
        w4int = wilhoit.integral_T4(Tint)
        w4min = wilhoit.integral_T4(Tmin)
        w4max = wilhoit.integral_T4(Tmax)

    if weighting:
        b[0] = 2*(wM1int - wM1min)
        b[1] = 2*(w0int - w0min)
        b[2] = 2*(w1int - w1min)
        b[3] = 2*(w2int - w2min)
        b[4] = 2*(w3int - w3min)
        b[5] = 2*(wM1max - wM1int)
        b[6] = 2*(w0max - w0int)
        b[7] = 2*(w1max - w1int)
        b[8] = 2*(w2max - w2int)
        b[9] = 2*(w3max - w3int)
    else:
        b[0] = 2*(w0int - w0min)
        b[1] = 2*(w1int - w1min)
        b[2] = 2*(w2int - w2min)
        b[3] = 2*(w3int - w3min)
        b[4] = 2*(w4int - w4min)
        b[5] = 2*(w0max - w0int)
        b[6] = 2*(w1max - w1int)
        b[7] = 2*(w2max - w2int)
        b[8] = 2*(w3max - w3int)
        b[9] = 2*(w4max - w4int)

    # solve A*x=b for x (note that factor of 2 in b vector and 10*10 submatrix of A
    # matrix is not required; not including it should give same result, except
    # Lagrange multipliers will differ by a factor of two)
    import scipy.linalg
    x = scipy.linalg.solve(A,b,overwrite_a=1,overwrite_b=1)

    nasa_low = NASA(Tmin=0, Tmax=0, coeffs=[x[0], x[1], x[2], x[3], x[4], 0.0, 0.0], comment='')
    nasa_high = NASA(Tmin=0, Tmax=0, coeffs=[x[5], x[6], x[7], x[8], x[9], 0.0, 0.0], comment='')

    return nasa_low, nasa_high

cpdef Wilhoit_to_NASA_TintOpt(Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int contCons):
    """
    Convert a Wilhoit polynomial to a pair of NASA polynomials, using an
    optimization algorithm to choose the best value of the intermediate
    temperature. The parameters are the same as for the :func:`Wilhoit_to_NASA`
    function.
    """
    import scipy.optimize
    Tint = scipy.optimize.fminbound(Wilhoit_to_NASA_TintOpt_objFun, Tmin, Tmax, args=(wilhoit, Tmin, Tmax, weighting, contCons))
    Tint = float(Tint) # fminbound returns a numpy.ndarray object
    #note that we have not used any guess when using this minimization routine
    #2. determine the bi parameters based on the optimized Tint (alternatively, maybe we could have Wilhoit_to_NASA_TintOpt_objFun also return these parameters, along with the objective function, which would avoid an extra calculation)
    (nasa1, nasa2) = Wilhoit_to_NASA(wilhoit, Tmin, Tmax, Tint, weighting, contCons)
    return nasa1, nasa2, Tint

cpdef double Wilhoit_to_NASA_TintOpt_objFun(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, bint weighting, int contCons):
    """
    Evaluate the objective function used to convert a Wilhoit polynomial to a
    pair of NASA polynomials. The parameters are the same as for the 
    :func:`Wilhoit_to_NASA` function.
    """
    if (weighting == 1):
        result = Wilhoit_to_NASA_TintOpt_objFun_W(Tint, wilhoit, Tmin, Tmax, contCons)
    else:
        result = Wilhoit_to_NASA_TintOpt_objFun_NW(Tint, wilhoit, Tmin, Tmax, contCons)

    # numerical errors could accumulate to give a slightly negative result
    # this is unphysical (it's the integral of a *squared* error) so we
    # set it to zero to avoid later problems when we try find the square root.
    if result < 0:
        logging.info("Negative ISE of %f reset to zero."%(result))
        result = 0

    return result

cpdef double Wilhoit_to_NASA_TintOpt_objFun_NW(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int contCons):
    """
    Evaluate the unweighted objective function used to convert a Wilhoit 
    polynomial to a pair of NASA polynomials. The parameters are the same as 
    for the :func:`Wilhoit_to_NASA` function.
    """
    cdef NASA nasa_low, nasa_high
    cdef double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10
    cdef double qM1, q0, q1, q2, q3, result

    nasa_low, nasa_high = Wilhoit_to_NASA(wilhoit,Tmin,Tmax,Tint, 0, contCons)
    b1, b2, b3, b4, b5 = nasa_low.c0, nasa_low.c1, nasa_low.c2, nasa_low.c3, nasa_low.c4
    b6, b7, b8, b9, b10 = nasa_high.c0, nasa_high.c1, nasa_high.c2, nasa_high.c3, nasa_high.c4

    q0 = wilhoit.integral_T0(Tint)
    q1 = wilhoit.integral_T1(Tint)
    q2 = wilhoit.integral_T2(Tint)
    q3 = wilhoit.integral_T3(Tint)
    q4 = wilhoit.integral_T4(Tint)
    result = (wilhoit.integral2_T0(Tmax) - wilhoit.integral2_T0(Tmin) +
                 nasa_low.integral2_T0(Tint) - nasa_low.integral2_T0(Tmin) +
                 nasa_high.integral2_T0(Tmax) - nasa_high.integral2_T0(Tint)
                 - 2* (b6*(wilhoit.integral_T0(Tmax)-q0)+b1*(q0-wilhoit.integral_T0(Tmin))
                 +b7*(wilhoit.integral_T1(Tmax) - q1) +b2*(q1 - wilhoit.integral_T1(Tmin))
                 +b8*(wilhoit.integral_T2(Tmax) - q2) +b3*(q2 - wilhoit.integral_T2(Tmin))
                 +b9*(wilhoit.integral_T3(Tmax) - q3) +b4*(q3 - wilhoit.integral_T3(Tmin))
                 +b10*(wilhoit.integral_T4(Tmax) - q4)+b5*(q4 - wilhoit.integral_T4(Tmin))))

    return result

cpdef double Wilhoit_to_NASA_TintOpt_objFun_W(double Tint, Wilhoit wilhoit, double Tmin, double Tmax, int contCons):
    """
    Evaluate the weighted objective function used to convert a Wilhoit 
    polynomial to a pair of NASA polynomials. The parameters are the same as 
    for the :func:`Wilhoit_to_NASA` function. The weighting is by inverse
    temperature, to bias the fit towards the lower temperatures, where the
    heat capacity is changing more rapidly.
    
    If the fit is close to perfect, the result may be slightly negative due to 
    numerical errors in evaluating this integral.
    """
    cdef NASA nasa_low, nasa_high
    cdef double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10
    cdef double qM1, q0, q1, q2, q3, result
    
    nasa_low, nasa_high = Wilhoit_to_NASA(wilhoit,Tmin,Tmax,Tint, 1, contCons)
    b1, b2, b3, b4, b5 = nasa_low.c0, nasa_low.c1, nasa_low.c2, nasa_low.c3, nasa_low.c4
    b6, b7, b8, b9, b10 = nasa_high.c0, nasa_high.c1, nasa_high.c2, nasa_high.c3, nasa_high.c4

    qM1 = wilhoit.integral_TM1(Tint)
    q0 = wilhoit.integral_T0(Tint)
    q1 = wilhoit.integral_T1(Tint)
    q2 = wilhoit.integral_T2(Tint)
    q3 = wilhoit.integral_T3(Tint)
    result = (wilhoit.integral2_TM1(Tmax) - wilhoit.integral2_TM1(Tmin) +
                 nasa_low.integral2_TM1(Tint) - nasa_low.integral2_TM1(Tmin) +
                 nasa_high.integral2_TM1(Tmax) - nasa_high.integral2_TM1(Tint)
                 - 2* (b6*(wilhoit.integral_TM1(Tmax)-qM1)+b1*(qM1 - wilhoit.integral_TM1(Tmin))
                 +b7*(wilhoit.integral_T0(Tmax)-q0)+b2*(q0 - wilhoit.integral_T0(Tmin))
                 +b8*(wilhoit.integral_T1(Tmax)-q1)+b3*(q1 - wilhoit.integral_T1(Tmin))
                 +b9*(wilhoit.integral_T2(Tmax)-q2)+b4*(q2 - wilhoit.integral_T2(Tmin))
                 +b10*(wilhoit.integral_T3(Tmax)-q3)+b5*(q3 - wilhoit.integral_T3(Tmin))))

    return result
