#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 07:16:41 2021

@author: leonardo
"""
from sympy import symbols, Matrix, BlockMatrix, eye
from sympy import integrate, diff
from sympy import latex, pprint
from sympy import DotProduct


[xi,eta,L,H] =symbols(['xi','eta','L','H'])
[q1,q2,q3,q4,q5,q6,q7,q8] = symbols(['q1','q2','q3','q4','q5','q6','q7','q8'])

q = Matrix([q1,q2,q3,q4,q5,q6,q7,q8])

# elemento quadrático
S1 = -2/L**2 * xi * (L/2 - xi)
S2 = eta * S1
S3 = 2/L**2 * xi * (L/2 + xi)
S4 = eta * S3
S5 = - 4/L**2 * (xi - L/2) * (xi + L/2)
S6 = eta * S5

S = Matrix([[S1,0,S2,0,S3,0,S4,0,S5,0,S6,0],
            [0,S1,0,S2,0,S3,0,S4,0,S5,0,S6]])

# elemento linear
S1 = (L/2 - xi)
S2 = eta * S1
S3 = (L/2 + xi)
S4 = eta * S3

S = 1/L * Matrix([[S1,0 ,S2,0 ,S3,0 ,S4,0],
                       [0 ,S1,0 ,S2,0 ,S3,0 ,S4]])

Sxxi = diff(S[0,:],xi)
Sxeta = diff(S[0,:],eta)
Syxi = diff(S[1,:],xi)
Syeta = diff(S[1,:],eta)


J = Matrix([[Sxxi*q,Sxeta*q],[Syxi*q,Syeta*q]])
F = J
eps = 1/2 * (F.T*F - eye(2))   

Deq = diff(eps,q)

