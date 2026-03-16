# file: jones_6D_inverse
# date: February 2026
# author: DAJones

"""
Note that `copy()` is not an available method for this implementation of multivector.  
Python treats a statement like `S=H` as providing an alias `S` for the multivector `H` 
when what is needed by the algorithm is an initialized alternate multivector with the same value as `H` 
which can be modified while the value of `H` remains intact.  
The construct `S=H+0` is used to force the intention.
"""

SMALL = 1e-12
SILENT = False

# this inverse routine should behave well for any multivector of 6 or less dimensions
def jones_inverse( A ):
    H = A * ~A
    if H.isScalar(): # catch scaled unitary
        d = H.value[0]
        if -SMALL < d and d < SMALL:
            if not SILENT: print( "This multivector is singular." )
            return None
        if not SILENT: print( "This multivector is scaled unitary." )
        return ~A * (1/H.value[0])
    else:
        N = A.layout.dims
        if N > 6:
            if not SILENT: print( "This inverse is limited to 6 dimensions." )
            return None # 6D is the limit
        else:
            S = H + 0   # independent copy 
        if N > 4: # 5,6
            S.value[0] *= -3
            S *= H
            S.value[0] *= -1
            S *= H
            S.value[0] *= (-1/3)
        else:  # N <= 4
            S.value[0] *= -1
    D = S * H
    d = D.value[0]
    if -SMALL < d and d < SMALL:
        if not SILENT: print( "This multivector is singular." )
        return None
    else:
        iH = S * (1/d)
        iA = ~A * iH
        return iA

# 6D Multivector Inverse for publically available package 'clifford' through conda
def jones_6D_inverse( A ):
    B = A * ~A
    C = B + 0     # C = B.copy()
    C.value[0] *= -3
    C *= B
    C.value[0] *= -1
    C *= B
    C.value[0] *= (-1/3)
    D = C * B   # D is scalar
    d = D.value[0]
    if -SMALL < d and d < SMALL:
        return None # singular
    return ~A * (C * (1.0 / d))
