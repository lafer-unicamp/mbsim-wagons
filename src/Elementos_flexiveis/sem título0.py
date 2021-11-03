#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 16:06:24 2021

@author: leonardo
"""
import multiprocessing as mp
import numpy as np
from time import time
from nachbagauer import elementQuadratic, node
from flexibleBody import flexibleBody
from materials import linearElasticMaterial

print("Number of cpu : ", mp.cpu_count())

steel = linearElasticMaterial('Steel',200e9,0.3,7.85e3)
biela = flexibleBody('Biela',steel)

L = 2.

nel = 2

nodesA = [0]*(nel*3-1)*4
for i in range(nel+1):
    nodesA[i]=node(L * i/nel,0.0,0.0,1.0)
    
elemA = [0]*nel
for i in range(nel):
    elemA[i]=elementQuadratic(nodesA[i],nodesA[i+1], 0.1, 0.1)
    

biela.addElement(elemA)

def elasticForces(q,qdot):
    
    biela.updateDisplacements(q)
    f = biela.assembleElasticForceVector()
    print(elemA[0].changedStates,elemA[1].changedStates)
    
    return f

def evaluateJacobian(function,x0,xdot0=0,parameter=0,dofCalculate=None):
    
       
    f0 = function(x0,xdot0)
    
    if dofCalculate is None:
        dofCalculate = range(f0.shape[0])
    
    J = np.zeros([f0.shape[0],x0.shape[0]])
    
    if parameter == 0:
        q = x0
    else:
        q = xdot0
    
    for i in dofCalculate:
        qsave = q[i]
        q[i] += 1e-6
        
        J[:,i] = (function(x0,xdot0) - f0).A1 * 1e6
        
        q[i] = qsave
    return J

ts = time()
J = evaluateJacobian(elasticForces,np.zeros(biela.totalDof))
print(time()-ts)
M = biela.assembleTangentStiffnessMatrix()
d = M-J