# -*- coding: utf-8 -*-
"""
util.py
Created on Sat Nov 24 09:40:24 2018

# 2025.07.31 DAJones: copied to AACA folder.
# 2025.08.01 DAJones: added optimized 6D container function

"""

import Clifford
#from numpy import random,array
import numpy
from math import sqrt

AccumFormat = '{0:08b} {1:15.3f}' 
Accum2Format = '{0:08b} {1:15.3f} {2:15.3f}'

# print a real multivector
def printAccum(A):
    for n in range(Clifford.bases()):
        print(AccumFormat.format(n,A.Reg[n]))

# print 2 real multivectors together
def print2Accum(A,B):
    for n in range(Clifford.bases()):
        print(Accum2Format.format(n,A.Reg[n],B.Reg[n]))

# print a complex multivector        
def printAccum2(A):
    for n in range(A.bases):
        print(Accum2Format.format(n,A.RegR[n],A.RegI[n]))

# randomly fill a multivector accumulator
def random(signature,dimensions):
#    Clifford.signature = signature
#    Clifford.dimensions = dimensions
    Clifford.Initialize(dimensions,signature)
    A = Clifford.Accum()
    samples = numpy.random.normal( 0.0, 1.0, Clifford.bases() )
    A.Reg = numpy.array(samples).tolist()
    return A

# select the grades given in v[] from the multivector accumulator S
def grade(S,v):
    T = S.copy()
    for n in range(Clifford.bases()):
        m = Clifford.Grade[n]
        if not(m in v): T.Reg[n] = 0.0
    return T

# negate the grades given in v[] in the multivector accumulator S
def involve(S,v):
    T = S.copy()
    for n in range(Clifford.bases()):
        m = Clifford.Grade[n]
        if (m in v): T.Reg[n] = - T.Reg[n]
    return T

def testClearGrade(A):
    V = [True]*(1+Clifford.dimensions)
    for n in range(Clifford.bases()):
        if 0.000001 < abs(A.Reg[n]):
            V[Clifford.Grade[n]]=False
    #for n in range(len(V)):
    #    print('{0:8n} {1:4}'.format(n,V[n]))
    return V

def testEquality(A,B):
    if A.dimensions != B.dimensions:
        return False
    if A.signature != B.signature:
        return False
    for n in range(A.bases):
        if 0.000001 < abs(A.Reg[n]-B.Reg[n]):
            return False
    return True

# return an equivalent Euclidean multivector with an extra dimension
# Clifford context is modified
def promote(A):
    signature = Clifford.signature # the presumed signature intended for A
    Clifford.dimensions += 1
    Clifford.signature = 0
    B = Clifford.Accum()
    for i in range(A.bases):
        if 0 != 1&Clifford.Grade[i&signature]:
            B.Reg[i+A.bases] = A.Reg[i]
        else:
            B.Reg[i] = A.Reg[i]
    return B

# =============================================================================
# Reduced multiplication load 6D inverse for Cl(6,0) only (not all signatures)

import InverseSupport_6D

# this function completes the optimized Euclidean 6D inverse
def I6O(A):
    G = Clifford.Accum()
    # pass just the array from A and stuff the array result into G
    G.Reg = InverseSupport_6D.I6OS(A.Reg) 
    return( ~A * G )

# =============================================================================
# Inverse routines for multivectors of particular dimension but arbitrary signature

# this algorithm was derived from 5D, so it looks like it
def I4A(A):
    '''returns the inverse of a 4D multivector'''
    B = A*~A # leaving grades 0,1,4
    C = B*involve(B,[0]) # leaving only grade 0
    # 1=C(C{0}/D)=BB{0,5}C{0}/D=A~AB{0,5}C{0}/D
    I = involve(B,[0]).scale(1/C.Reg[0]) # 1/B
    return ~A*I # 1/A

# this is the same algorithm rewritten for comparison to I6()
def I4(A):
    '''returns the inverse of a 4D multivector'''
    B = A*~A # leaving grades 0,1,4
    G = B.copy()
    G.Reg[0] *= -1
    I = G.copy()
    G *= B # G is now a scalar
    # 1 = B(I/G) = A(~AI/G)
    return ~A*I.scale(1/G.Reg[0])

def I5(A):
    '''returns the inverse of a 5D multivector'''
    B = A*~A # leaving grades 0,1,4,5
    C = B*involve(B,[0,5]) # leaving grades 0,5
    D = C*involve(C,[0]) # leaving only grade 0
    # 1=C(C{0}/D)=BB{0,5}C{0}/D=A~AB{0,5}C{0}/D
    I = involve(C,[0]).scale(1/D.Reg[0]) # 1/C
    I = involve(B,[0,5])*I # 1/B
    return ~A*I # 1/A

def I6(A):
    '''returns the inverse of a 6D multivector'''
    B = A*~A # lacking grades 2,3,6,7, leaving 0,1,4,5
    G = B.copy()
    G.Reg[0] *= -3
    G *= B
    G.Reg[0] *= -1
    G *= B
    G.Reg[0] *= (-1/3)
    I = G.copy()
    G *= B # G is now a scalar
    # 1 = B(I/G) = A(~AI/G)
    return ~A*I.scale(1/G.Reg[0])
    
