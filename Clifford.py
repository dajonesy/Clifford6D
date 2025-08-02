# ======================================================
# Clifford Module agreeing with Sangwine implementation, and cleaned of development clutter
# ------------------------------------------------------
# 2025.07.31 DAJones: copied to AACA folder.
# 2025.08.01 DAJones: reversed the multiplication order to agree with Sangwine
# ======================================================

import math

# restrict the number of dimensions because it effects accumulator size
dimensions = 0

# Use bases() to find out how many registers are in an accumulator
def bases():
    return 1<<dimensions

# Use mask() to present the important bits in an element type
def mask():
    return (1<<dimensions)-1

# The grade of an element is so useful, we table it here
Grade = [0]

# if the number of dimensions changes, Grade[] needs to be recomputed
def ComputeGrade():
    global Grade
    Grade = [0]
    n = 0
    # to support signature, we need to include an extra dimension for "adjustment"
    while len(Grade) < (bases()<<1):
        Grade += [g+1 for g in Grade]
    return

# ---------------------------------------------------------------------
# Set the signature bit for dimensions you want to consider negative
# (signature=-1 sets up anti-Euclidean space regardless of dimension)
# a conceptual extra dimension can serve for anti-Euclidean components
signature = 0

def Initialize( d, s=0 ):
    global dimensions
    global signature
    # order is important
    dimensions = d
    signature = s
    ComputeGrade()
    return

# ======================================================
# when multiplying two elements, the rearrangement of primitives is immitated to count the swaps needed
# built to Sangwine's specification
def swapCount( bitsLeft, bitsRight ):
    # initialize the swap count with the number of common elements both left and right and negative signature
    nSwap = Grade[signature & bitsLeft & bitsRight]
    nJump = Grade[ bitsRight ] 
    # the number of high order vecs in the right pile for a lower order vec from the left to jump over
    mask = 1 # one bit at a time, beginning on the right
    commons = bitsLeft & bitsRight # mark the common elements left and right which have negative signature
    while 0 != nJump:
        if 0 != (mask & bitsRight): # found low order vec on the right
            nJump -= 1 # leave it on the right
        if 0 != (mask & bitsLeft): # found low order vec on the left
            nSwap += nJump # it moves over this many high order vecs on the right
        mask <<= 1 # next bit
    return nSwap

# ======================================================
# return a matrix of +-1 describing only the sign part of the multiplication table
# side effect printing is controlled version, which helps with documentation
SignTableVersion = 2 # turn off side effect printing
def SignTable():
    signTable = []
    for i in range(bases()):
        line = ''
        signLine = []
        for j in range(bases()):
            n = swapCount(i,j)
            if SignTableVersion == 1: line=line+(['  ',' -'][1&n])     # version 1
            if SignTableVersion == 2: line=line+(['0','1'][1&n])       # version 2 (the Boolean version)
            if SignTableVersion == 3: line=line+([' 1',' -'][1&n])     # version 3
            if SignTableVersion == 4: line=line+([' 1,','-1,'][1&n])   # version 4 (the numerical version)
            signLine.append([1,-1][1&n])
        if SignTableVersion != 0: print(line)
        signTable.append(signLine)
    return signTable
    
# ======================================================
# print a basis multiplication table for dimensions < 9
def Table():
    for i in range(bases()):
        line = ''
        for j in range(bases()):
            n = swapCount(i,j)
            line = line + ([' ','-'][n&1])
            k = i ^ j
            line = line + (format(k,'02x'))
            line = line + ' '
        print( line )
    return

# ======================================================
# Elements are primitive components of multivectors
class Element:

    def __init__(self, type, value):
        self.type = type
        self.value = value
  
    def __str__(self):
        return "({0:02x},{1:8.3f})".format(self.type, self.value)
    
    def __mul__(self, other):
        type = self.type ^ other.type
        value = self.value * other.value
        n = swapCount(self.type, other.type)
        if 0 != (1&n): # negate for odd number of swaps
            value = -value
        return Element(type, value)

# ======================================================
# an accumulator implements a multivector equating type and index
# the main Clifford purpose for this class is to ease multivector multiplication
class Accum:

    def __init__(self):
        self.Reg = [0.0]*bases()
        # these are late additions in case I need to interrogate init conditions
        self.dimensions = dimensions
        self.signature = signature
        self.bases = bases()

    FORMAT = '{0:08b} {1:16.8f}\n' 
    def __str__(self):
        s = ''
        for i in range(self.bases):
            if (abs(self.Reg[i]) < Accum.SMALL):
                continue
            s += self.FORMAT.format(i,self.Reg[i])
        return s
    
    def __add__(self, other):
        A = Accum()
        for i in range(self.bases):
            A.Reg[i] = self.Reg[i] + other.Reg[i]
        return A

    def __sub__(self, other):
        A = Accum()
        for i in range(self.bases):
            A.Reg[i] = self.Reg[i] - other.Reg[i]
        return A

    def __neg__(self):
        A = Accum()
        for i in range(self.bases):
            A.Reg[i] = - self.Reg[i]
        return A
    
    # this is self*other, not other*self
    # when one of the multivectors is sparse it is best to put it on the 
    def __mul__(self, other):
    # no inputs are modified
        A = Accum()
        for i in range(len(other.Reg)):
            otherValue = other.Reg[i]
            if 0.0 == otherValue: 
                continue
            for j in range(self.bases):
                type = i ^ j
                value = otherValue * self.Reg[j]
                if 0!=(1&swapCount(j,i)): # order is important
                    value = -value
                A.Reg[type] += value
        return A
    # other in the outer loop presumes other is more remote
    # we look at other data one at a time
    # the more local self registers are repetitivly accessed.
    # the cost of communicating the other register array may be
    # balanced with the additional local processing required
        
    SMALL = 0.00000001
    def __eq__(self, other):
        for i in range(self.bases):
            if abs(self.Reg[i] - other.Reg[i]) > self.SMALL:
                return False
        return True
        
    def addElement(self, element):
        self.Reg[element.type] += element.value
        return
    
    def clear(self):
        for i in range(self.bases):
            self.Reg[i] = 0.0
        return
    
    def copy(self):
        other = Accum()
        for i in range(self.bases):
            other.Reg[i] = self.Reg[i]
        return other

    # involutions - see Lounesto page 56
    def __invert__(self):
        other = self.copy()
        for i in range(self.bases):
#           # while this makes some kind of sense, it does not agree with Sangwine
#           j = i & ~signature # mask off the negative signature (independent bivector) elements
#           if (3 & Grade[j]) > 1 : # that is when 3&Grade[j] in [2,3]
            if (3 & Grade[i]) > 1 : # that is when 3&Grade[i] in [2,3]
                other.Reg[i] = - other.Reg[i]
        return other
                    
    def reverse(self): # ~ 'reversion' (Lounesto uses ~)
        return ~self

    def conjugate(self): # ! 'Clifford conjugation' (Lounesto uses -)
        other=self.copy()
        for i in range(self.bases):
            n = 3 & Grade[i]
            if (2==n) | (1==n):
                other.Reg[i] = - other.Reg[i]
        return other
    
    def automorph(self): # !~ or ~! 'grade involution' (Lounesto uses ^)
        other=self.copy()
        for i in range(self.bases):
            n = 3 & Grade[i]
            if (3==n) | (1==n):
                other.Reg[i] = - other.Reg[i]
        return other

    def Make3DScalar(self):
        B = self * self.reverse()
        C = B * B.conjugate()
        return C
    
    def mag(self):
        a = 0.0
        for i in range(self.bases):
            value = self.Reg[i]
            a = a + value*value
        return math.sqrt(a)

    def scale(self,scalar):
        other = Accum()
        for i in range(self.bases):
            other.Reg[i] = scalar*self.Reg[i]
        return other

    def isomorph(self):
        # assume that signature == (dimensions-1)<<1
        A = Accum2(dimensions-1)
        for i in range(A.bases):
            A.RegR[i] = self.Reg[i]
            A.RegI[i] = self.Reg[i+self.signature]
        return A
    
    def normalize(self):
        r = 0.0
        for i in range(self.bases):
            r += self.Reg[i]*self.Reg[i]
        R = math.sqrt(r)
        RR = 1/R
        for i in range(self.bases):
            self.Reg[i] *= RR
        return R

# =============================================================================

# Complex Clifford accumulator
class Accum2:

    def __init__(self, d):
        self.dimensions = d
        self.bases = 1<<d
        self.RegR = [0.0]*self.bases
        self.RegI = [0.0]*self.bases
    
    Format = '{0:08b} {1:16.8f} {2:16.8f}\n' 
    def __str__(self):
        s = ''
        for i in range(self.bases):
            s += self.Format.format(i,self.RegR[i],self.RegI[i])
        return s

    def __add__(self, other):
        A = Accum2(self.dimensions)
        for i in range(self.bases):
            A.RegR[i] = self.RegR[i] + other.RegR[i]
            A.RegI[i] = self.RegI[i] + other.RegI[i]
        return A

    def __sub__(self, other):
        A = Accum2(self.dimensions)
        for i in range(self.bases):
            A.RegR[i] = self.RegR[i] - other.RegR[i]
            A.RegI[i] = self.RegI[i] - other.RegI[i]
        return A

    def copy(self):
        A = Accum2(self.dimensions)
        for i in range(self.bases):
            A.RegR[i] = self.RegR[i]
            A.RegI[i] = self.RegI[i]
        return A

    # involutions - see Lounesto page 56
    def __invert__(self):
        other=self.copy()
        for i in range(self.bases):
            n = 3 & Grade[i]
            if (3==n) | (2==n):
                other.RegR[i] = - other.RegR[i]
                other.RegI[i] = - other.RegI[i]
        return other
    
    def reverse(self): # ~ 'reversion' (Lounesto uses ~)
        return ~self
    
    def conjugate(self): # ! 'Clifford conjugation' (Lounesto uses -)
        other=self.copy()
        for i in range(self.bases):
            n = 3 & Grade[i]
            if (2==n) | (1==n):
                other.RegR[i] = - other.RegR[i]
                other.RegI[i] = - other.RegI[i]
        return other
    
    def automorph(self): # !~ or ~! 'grade involution' (Lounesto uses ^)
        other=self.copy()
        for i in range(self.bases):
            n = 3 & Grade[i]
            if (3==n) | (1==n):
                other.RegR[i] = - other.RegR[i]
                other.RegI[i] = - other.RegI[i]
        return other

    def scale(self, scalar):
        A = Accum2(self.dimensions)
        for i in range(self.bases):
            A.RegR[i] = scalar*self.RegR[i]
            A.RegI[i] = scalar*self.RegI[i]
        return A

    def scale2(self, real, imag):
        A = Accum2(self.dimensions)
        for i in range(self.bases):
            tempR = real*self.RegR[i] - imag*self.RegI[i]
            tempI = real*self.RegI[i] + imag*self.RegR[i]
            A.RegR[i] = tempR
            A.RegI[i] = tempI
        return A

    # this is self*other, not other*self
    def __mul__(self, other):
    # no inputs are modified
        A = Accum2(self.dimensions)
        for i in range(self.bases):
            otherValueR = other.RegR[i]
            otherValueI = other.RegI[i]
            for j in range(self.bases):
                type = i ^ j
                valueR = (otherValueR * self.RegR[j]) - (otherValueI * self.RegI[j])
                valueI = (otherValueR * self.RegI[j]) + (otherValueI * self.RegR[j])
                if 0 != (1&swapCount2(i,j)): # order is important
                    valueR = -valueR
                    valueI = -valueI
                A.RegR[type] += valueR
                A.RegI[type] += valueI
        return A

    # convert from complex to real
    def isomorph(self):
        # enforce determinant externals
        global dimensions
        dimensions = self.dimensions+1
        global signature
        signature = self.bases
        A = Accum()
        for i in range(self.bases): # assume accum2.bases==bases()>>1
            A.Reg[i] = self.RegR[i]
            A.Reg[i+self.bases] = self.RegI[i]
        return A            
        
# ============================================================================= 
# back to global functions
#
# the following functions are offered: for testing, and in example
# ============================================================================= 

def printInvolutes(A):
    rA = A.reverse()   # ~A
    cA = A.conjugate() # !A
    aA = A.automorph() # ~!A
    # (A*!A)*~(A*!A) = A*(!A*~!A*~A) = A*I
    # I*A = !A*~!A*~A*A = !(~A*A)*(~A*A)
    I = A.conjugate()*A.automorph()*A.reverse() # (!A*~!A*~A)
    AI = A*I
    I = I.scale(1.0/AI.mag())
    AA = A*A
    ArA = A*A.reverse()   # A*~A
    AcA = A*A.conjugate() # A*!A
    AaA = A*A.automorph() # A*~!A
    AI = A*I              # A*(!A*~A*~!A)
    IA = I*A
    print( '                A       ~A       !A      ~!A        I      A*A     A*~A     A*!A    A*~!A      A*I      I*A')
    for n in range(len(A.Reg)):
        print( '{0:08b} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f} {6:8.3f} {7:8.3f} {8:8.3f} {9:8.3f} {10:8.3f} {11:8.3f}'.
            format(n,A.Reg[n],rA.Reg[n],cA.Reg[n],aA.Reg[n],I.Reg[n],AA.Reg[n],ArA.Reg[n],AcA.Reg[n],AaA.Reg[n],AI.Reg[n],IA.Reg[n]))
    print( 'magnitude{0:8.3f} {1:8.3f} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f} {6:8.3f} {7:8.3f} {8:8.3f} {9:8.3f} {10:8.3f}'.
       format(A.mag(),rA.mag(),cA.mag(),aA.mag(),I.mag(),AA.mag(),ArA.mag(),AcA.mag(),AaA.mag(),AI.mag(),IA.mag()))

# ============================================================================
# gross testing of these classes:
# The class Accum is most important for execution
# The class Element is more important to user interaction
def Test():
    e=Element(0,1.0)
    print("scalar unit e=Element",e)
    a=Element(1,2.0)
    print("vector a=Element",a)
    b=Element(2,3.0)
    print("vector b=Element",b)
    c=a*b
    print("bivector a*b=Element",c)
    d=b*a
    print("bivector b*a=Element",d,"demonstrating anti-commutivity of bivectors")

    A=Accum()
    print("users promise an interest in no more than", dimensions, "dimensions.")
    print("Clifford numbers of this dimension require", bases(), "registers")
    print("initially",A.Reg)
    A.addElement(a)
    A.addElement(b)
    A.addElement(d)
    A.addElement(e)
    print("accumulating 1+a+b+b*a results in",A.Reg)
    B=A.reverse()
    print("reverse:",B.Reg)
    B=A.conjugate()
    print("conjugate:",B.Reg)
    B=A.automorph()
    print("automorph:",B.Reg)
    