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

massa1 = 1.0
massa2 = 1.0

def spring(x,u=0):
    
    f1 = - 30. * x[0] + 10. * (x[1] - x[0])
    f2 = - 10. * (x[1] - x[0])
    
    return np.array([f1,f2])


dt = 0.5e-3
tf = 3.

t = [0]
q = []
u = []

q0 = np.array([0.,1.])
u0 = np.array([0.,0.])

q.append(q0)
u.append(u0)

    
def evaluateJacobian(function,x0,parameter):
    
    f0 = function(x0)
    
    J = np.zeros([f0.shape[0],len(x0)])
    
    for i,col in enumerate(J.T):
        x0save = x0[i]
        x0[i] += 1e-6
        
        col[:] = (function(x0) - f0) * 1e6
        
        x0[i] = x0save
        
    return J
        
        
    
    
def integrateEulerSemiImplicit(q0,u0,dt,M,extForce,Phi,tfinal):
    
    t = [0]
    q = []
    u = []
    lam = []
    
    ndof = M.shape[0]
    ncons = Phi.shape[0]
    
    LHS = np.zeros([ndof+ncons,ndof+ncons])
    RHS = np.zeros(ndof+ncons)
    
    q.append(q0)
    u.append(u0)
    lam.append(np.zeros(ncons))
    
    while t[-1] <= tfinal:
        
        Jq =  evaluateJacobian(extForce,q[-1],0)
        
        # explicit update for positions
        q.append(q[-1] + dt*u[-1])
        
        
        LHS[0:ndof,0:ndof] = M - dt * dt * Jq
        LHS[0:ndof,ndof:] = Phi.T
        LHS[ndof:,0:ndof] = Phi
        
        RHS[0:ndof] = dt * extForce(q[-1]) + dt * dt * np.dot(Jq,u[-1])
        RHS[ndof:] = -np.dot(Phi,u[-1])
        
        incr = np.linalg.solve(LHS,RHS)
        du = incr[0:ndof]
        lam.append(incr[ndof:])
        
        # implicit update for velocities
        u.append(u[-1] + du)
        
        t.append(t[-1] + dt)
        
    return t,q,u,lam
    
M = np.diag([massa1,massa2])
    
t,q,u,lam = integrateEulerSemiImplicit(q0, u0, dt, M, spring, np.matrix([[0,1]]), tf)
    
plt.subplot(3,1,1)
plt.plot(t,q)
plt.ylabel('x')
plt.subplot(3,1,2)
plt.plot(t,u)
plt.ylabel('v')
plt.subplot(3,1,3)
plt.plot(t,lam)
plt.ylabel('force')

    
    
