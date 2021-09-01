#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 07:16:41 2021

@author: leonardo
"""
from sympy import symbols, Matrix
from sympy import integrate
from sympy import latex

[xi,eta,L,H] =symbols(['xi','eta','L','H'])


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

I = integrate(integrate(S,(xi,-L/2,L/2)),(eta,-H/2,H/2)) * 4/(H*L)

out = latex(I,mat_str='bmatrix',mat_delim='')

