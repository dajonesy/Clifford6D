# file: acus_6D_inverse.py
# date: March 2026
# author: DAJones

# This is a 6D inverse based on the work of Acus & Dargys 

SMALL = 1e-12
SILENT = True

def bar( H ):
    return 2*H.value[0] - H

# notice that S() factors out an H/3 from the determinant
def S( H ):
    Hb = bar(H)
    return H * bar(H*H) + 2 * bar(Hb * bar(Hb*Hb))

# -----+---------------------+------------------------------+
# dim  |    determinant      | A *      G / D  == 1         |
# -----+---------------------+------------------------------+
# 5,6  | (1/3)(A*~A)S(A*~A)  | A * ~AS(H) / HS(H),  (H=A~A) |
# 3,4  | (1/3)AS(A)          | A *   S(A) / AS(A)           |
# 1,2  | A*bar(A)            | A * bar(A) / Abar(A)         |
# any  | scalar              | A *     ~A / H       (H=A~A) | 
# -----+---------------------+------------------------------+

def acus_inverse(A):
    H = A * ~A
    if H.isScalar(): # catch scaled unitary
        if not SILENT: print( "This multivector is scaled unitary." )
        G = ~A
        D = H
    else:
        N = A.layout.dims
        if N > 6:
            if not SILENT: print( "This inverse is limited to 6 dimensions." )
            return None # 6D is the limit
        elif N > 4: # 5,6
            G = ~A * S(H)
        elif N > 2: # 3,4
            G = S(A)
        else: # N == 1,2
            G = bar(A)
        D = A * G
    d = D.value[0]
    if -SMALL < d and d < SMALL:
        if not SILENT: print( "This multivector is singular." )
        return None
    else:
        return G * (1.0/d)
 
def acus_6D_inverse(A):
    H = A * ~A
    Hb = bar(H)
    S = H * bar(H*H) + 2 * bar(Hb * bar(Hb*Hb))
    G = ~A * S
    D = A * G # D is scalar
    d = D.value[0]
    if -SMALL < d and d < SMALL:
        return None
    else:
        return G * (1.0/d)
