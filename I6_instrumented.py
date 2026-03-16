# file: I6_instrumented.py
# started: February 2026

# A fully instrumented version of the 6D inverse algorithm
# for use with the publically available clifford package.
# Includes detection of (scaled) unitary and singular cases.

import clifford
#import numpy
#import numba

# quiet deprecation warnings
#from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
#import warnings
#warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
#warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

# adaptable constants
VERBOSE = True # when True, full instrumentation is provided
SMALL = 1e-12

# 6D Multivector Inverse for the publically available package 'clifford'
def I6_instrumented( A ):
    if VERBOSE: print( "A = ", str(A) )
    H = A * ~A
    if VERBOSE: print( "H = A * ~A = ", str(H) )
    if H.isScalar():
        d = H.value[0]
        if d < SMALL and d < SMALL: return None
        iA = ~A * (1.0/d)
        if VERBOSE: print( "Multivector is Scaled Unitary\n1/A = ", str(iA) )
        return iA
    
    S1 = H + 0     # S1 = H.copy()
    S1.value[0] *= -3
    S1 *= H
    if VERBOSE: print( "S1 = (H o 3) H = ", str(S1) )

    S0 = S1 + 0    # S0 = S1.copy()
    S0.value[0] *= -1
    S0 *= H
    if VERBOSE: print( "S0 = (S1 o 1) H = ", str(S0) )

    S = S0 + 0    # S = S0.copy()
    S.value[0] *= (-1/3)
    if VERBOSE: print( "S = S0 o 1/3 = ", str(S) )
    
    D = S * H
    if VERBOSE: print( "D = S * H = ", str(D) )

    d = D.value[0]
    if -SMALL < d and d < SMALL: 
        if VERBOSE: print( "Multivector has no inverse." )
        return None
    iH = S * (1.0/d)
    if VERBOSE: print( "1/H = S/D = ", str( iH ) )
    if VERBOSE: print( "H * 1/H = ", str( H * iH ) )

    iA = ~A * iH
    if VERBOSE: print( "1/A = ~A * 1/H = ", str( iA ) )
    if VERBOSE: print( "A * 1/A = ", str( A * iA ) )
    
    return iA

# E X A M P L E
if __name__ == "__main__":

    layout, blades = clifford.Cl(3,3)
    locals().update(blades)

    print( "\nA COMPACT 6D TEST" )
    A = e16 + e25 + e34 + e124 + e135 + e236 + e123456
    iA = I6_instrumented( A )
