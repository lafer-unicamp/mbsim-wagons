#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 07:36:03 2021

@author: leonardo
"""

from nachbagauer import elementQuadratic, node
from materials import linearElasticMaterial
from flexibleBody import flexibleBody
import numpy as np
from matplotlib import pyplot as plt

np.seterr('raise')

'''
TESTE DOS ELEMENTOS EM BIELA MANIVELA
'''



steel = linearElasticMaterial('Steel',200e9,0.3,7.85e3)
biela = flexibleBody('Biela',steel)

# comprimentos
Lmanivela = 150.0e-3
Lbiela = 300.0e-3

# seção transversal da biela
diametro = 6.0e-3

# gravidade
g = [0,-9.810*0]

# malha
nel = 2
nodesA = [0]*(nel*3-1)*4
for i in range(nel+1):
    nodesA[i]=node(Lmanivela + Lbiela * i/nel,0.0,0.0,1.0)
    
elemA = [0]*nel
for i in range(nel):
    elemA[i]=elementQuadratic(nodesA[i],nodesA[i+1], 0.876*diametro, 0.876*diametro)
    
    
biela.addElement(elemA)


'''
Manivela
'''

massaManivela = np.pi * diametro**2 / 4 * Lmanivela * steel.rho
IzManivela = massaManivela * (0.25  * diametro**2 / 4 + 0.334  * Lmanivela**2)

'''
Pistão
'''

massaPistao = massaManivela



'''
RESTRIÇõES
'''

def Phi(q,u):
    
    Phi = np.zeros([5,q.shape[0]])
    
    # PINO ENTRE MANIVELA E BIELA
    Phi[0,0] = Lmanivela*np.sin(q[0])
    Phi[1,0] = -Lmanivela*np.cos(q[0])
    Phi[0,1] = 1.
    Phi[1,2] = 1.
    
    # PINO ENTRE BIELA E PISTÃO
    Phi[2,-5] = 1.
    Phi[2,-1] = -1.
    
    # PRISMÁTICA
    Phi[3,-4] = 1.
    
    # VELOCIDADE ANGULAR CONSTANTE
    Phi[4,0] = 1.
    
    return Phi


'''
SOLVER
'''


dt = 1e-5
tf = 9 / 150

t = [0]
q = []
u = []

q0 = np.array([0.]*(2+biela.totalDof))
q0[-1] = Lbiela + Lmanivela
u0 = np.array([0.]*(2+biela.totalDof))
u0[0] = 150  # velocidade angular da manivela
u0[2] = 150 * Lmanivela
u0[6] = 150 * Lmanivela * 0.75
u0[10] = 150 * Lmanivela * 0.50
u0[14] = 150 * Lmanivela * 0.25

q.append(q0)
u.append(u0)


def forcas(q,u):
    
    forcaPesoManivela = massaManivela * g[1] * Lmanivela / 2 * np.cos(q[0])
    
    forcaPesoBiela = biela.assembleWeightVector(g).A1
    
    biela.updateDisplacements(q[1:])
    forcaElasticaBiela = biela.assembleElasticForceVector().A1
    
    forcas = np.copy(q) * 0
    
    forcas[0] = forcaPesoManivela
    forcas[1:-1] = forcaPesoBiela.T - forcaElasticaBiela
    
    return forcas


def perturbaForca(function,f0,x0,xdot0):
    qsave = q[i]
    q[i] += 1e-6
    
    rigidez = (function(x0,xdot0) - f0) * 1e6
    
    q[i] = qsave
    
    return rigidez

import multiprocessing
def evaluateJacobian(function,x0,xdot0=0,parameter=0):
    
    f0 = function(x0,xdot0)
    
    J = np.zeros([f0.shape[0],x0.shape[0]])
    
    if parameter == 0:
        q = x0
    else:
        q = xdot0
    
    for i,col in enumerate(J.T):
        qsave = q[i]
        q[i] += 1e-6
        
        col[:] = (function(x0,xdot0) - f0) * 1e6
        
        q[i] = qsave
        

        
    return J


t = []
q = []
u = []
lam = []
probe = []

def output(t_,q_,u_,lam_):
    t.append(t_)
    q.append(q_)
    u.append(u_)
    lam.append(lam_)
    
    posi = biela.plotPositions(show=False)
    posiP = np.mean(posi,axis=0).A1
    posiA = posi[0]
    posiAB = posi[-1] - posiA
    
    Ydef = posiP - (posiA + posiAB/2)
    
    probe.append(Ydef.A1)
    



def integrateEulerSemiImplicit(q0,u0,dt,M,extForce,Phi,tfinal,writeOutput,dtout=1e-3):
    
    tout = 0.
    
    tn = 0.
    qn = q0.copy()
    un = u0.copy()
    
    # constraint Jacobian
    Cqn = Phi(qn,un)
    
    ndof = M.shape[0]
    ncons = Cqn.shape[0]
    
    LHS = np.zeros([ndof+ncons,ndof+ncons])
    RHS = np.zeros(ndof+ncons)
    
    writeOutput(tout,q0,u0,np.zeros(ncons))
   
    
    Jq = evaluateJacobian(extForce,qn,un,0)
    Ju = evaluateJacobian(extForce,qn,un,1)
    
    jacoCounter = 0
    
    while tn <=tfinal:
        
        
        # explicit update for positions
        deltaq = dt*un
        qnp1 = qn + deltaq
        
        # constraint Jacobian
        Cqn = Phi(qn,un)
        Cqnp1 = Phi(qnp1,un)
        
        # assemble implicit problem
        LHS[0:ndof,0:ndof] = M - dt * Ju - dt * dt * Jq
        if Cqn.shape[0] != 0:
            LHS[0:ndof,ndof:] = Cqn.T
            LHS[ndof:,0:ndof] = Cqnp1
        
        RHS[0:ndof] = dt * extForce(qn,un) + dt * dt * np.dot(Jq,un)
        if Cqn.shape[0] != 0:
            RHS[ndof:] = -np.dot(Cqnp1,un)
            RHS[-1] += 150 
        
        incr = np.linalg.solve(LHS,RHS)
        du = incr[0:ndof]
        lamnp1 = incr[ndof:]/dt
        
        
        # implicit update for velocities
        unp1 = un + du
        
        # jacobian update using Broyden's method
        dJq = (extForce(qnp1,unp1) - extForce(qn,un) - np.dot(Jq,deltaq)) 
        
        if np.linalg.norm(dJq) > 50:
            Jq = evaluateJacobian(extForce,qnp1,unp1,0)
            #Ju = evaluateJacobian(extForce,q[-1],u[-1],1)
            if jacoCounter < 5:
                # if the number of steps without jacobian reevaluations is smaller than N, reduce time step
                dt = dt/2
                print('Reduced time step to {0:1.4e}'.format(dt))
            print('Reevaluate jacobian after {0} steps at {1:1.8f} s'.format(jacoCounter,tn))
            jacoCounter = 0
        else:
            jacoCounter += 1
            Jq += dJq / (np.dot(deltaq,deltaq)+2e-16) * deltaq.T
            #Ju += (extForce(q[-1],u[-1]) - extForce(q[-2],u[-2]) - np.dot(Ju,du)) / (np.dot(du,du)) * du.T
        
        if jacoCounter > 20:
                dt = 1.15 * dt
                print('Increased time step to {0:1.4e}'.format(dt))
        
        tnp1 = tn + dt
        
        
        if tnp1 - tout > dtout:
            writeOutput(tnp1,qnp1,unp1,lamnp1)
            tout = tnp1
            
            
        # updates states
        qn = qnp1
        un = unp1
        tn = tnp1
            
    
M = np.zeros([2+biela.totalDof,2+biela.totalDof])
M[0,0] = IzManivela
M[1:-1,1:-1] = biela.assembleMassMatrix()
M[-1,-1] = massaPistao
    
t,q,u,lam = integrateEulerSemiImplicit(q0, u0, dt, M, forcas, Phi, tf, output, dtout=2e-4)
