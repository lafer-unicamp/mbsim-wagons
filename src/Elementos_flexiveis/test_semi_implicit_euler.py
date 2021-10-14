# -*- coding: utf-8 -*-
"""
TESTE COM MÉTODO DE EULER SEMI-IMPLÍCITO

BURGERMEISTER, B.; ARNOLD, M.; ESTERL, B. DAE time integration for real-time
 applications in multi-body dynamics. ZAMM, v. 86, n. 10, 
 p. 759–771, 24 out. 2006. 


Created on Wed Oct 13 14:11:13 2021

@author: lbaru
"""


import numpy as np
import matplotlib.pyplot as plt

massa = 1.0

def spring(x,u=0):
    if x >= 0:
        return - 10.0 * x - 0.1 * u
    else:
        return - 10e4 * x


dt = 0.5e-2
tf = 3.0

t = [0]
q = []
u = []

q0 = 0
u0 = 3

q.append(q0)
u.append(u0)

while t[-1] < tf:
        
    if q[-1] == 0:
        Jq = 0
    else:
        Jq = 1/massa * spring(q[-1]) / q[-1]
    
    q.append(q[-1] + dt*u[-1])
    
    du = (dt * spring(q[-2]) + dt**2 * Jq * u[-1])
    
    u.append(u[-1] + du)
    
    t.append(t[-1] + dt)
    
def integrateEulerSemiImplicit(q0,u0,dt,M,extForce,Phi,tfinal):
    
    t = [0]
    q = []
    u = []
    
    ndof = M.shape[0]
    
    Minv = np.linalg.inv(M)
    
    while t[-1] <= tfinal:
        
        Jq =  np.dot(Minv, 
    
    
plt.subplot(2,1,1)
plt.plot(t,q)
plt.ylabel('x')
plt.subplot(2,1,2)
plt.plot(t,u)
plt.xlabel('v')

    
    
