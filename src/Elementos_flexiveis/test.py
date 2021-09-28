# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 15:40:45 2021

@author: lbaru
"""

from nachbagauer import elementLinear, elementQuadratic, node
from materials import linearElasticMaterial
from flexibleBody import flexibleBody
import numpy as np
from scipy.optimize import fsolve

'''
TEST PROGRAM
'''


steel = linearElasticMaterial('Steel',207e3,0.3,7.85e-6)
body = flexibleBody('Bar',steel)
body2 = flexibleBody('Bar quadratic',steel)


n = []
nq = []
nel = 2
totalLength = 2000.
for i in range(nel+1):
    n.append(node(totalLength * i/nel,0.0,0.0,1))
    nq.append(node(totalLength * i/nel,0.0,0.0,1))


e = []
eq = []
for j in range(len(n)-1):
    e.append(elementLinear(n[j],n[j+1],500,100))
    eq.append(elementQuadratic(nq[j],nq[j+1],500,100))


g = np.matrix([[0.0,-9.81]])

body.addElement(e)
body.assembleMassMatrix()
body2.addElement(eq)
body2.assembleMassMatrix()




# Constraint matrix 
simBody = body2
conDof = [0,1,2,3]
gdl = simBody.totalDof
Phi = np.zeros([len(conDof),gdl])

# fixed bc
for d in conDof:
    Phi[d,d] = 1

record  = []
def f(z):
    
    ndof = gdl
    
    x = z[0:ndof]
    lam = z[ndof:]
    
    for ele in simBody.elementList:
        for nd in ele.nodes:
            nd.q = x[nd.globalDof]
            
    Qe = simBody.assembleElasticForceVector()
    Qa = Qe*0
    Qa[-3,0] = 5.0e8 * 0.5 * 0.5 * 0.5
    
    goal = [0]*(gdl+4)
    
    goal[0:ndof] = (-Qe.T + Qa.T + np.dot(Phi.T,lam)).tolist()[0]
    
    goal[ndof:] = np.dot(Phi,x).tolist()
    
    #print(f"max(|F|) = {np.max(np.abs(goal)):.8e}")
    
    record.append(goal)
    
    return goal


z0 = [0]*(gdl+4)
z0[-3] =  - 5.0e8 * 0.5 * 0.5 * 0.5
z0[-2] = - z0[-3] * 2000
#z = opt.newton_krylov(f,z0,maxiter=40,f_tol=1e-4,verbose=True)

z, info, ier, msg = fsolve(f, z0,full_output=True)
print(msg)

xy = simBody.plotPositions()
tipDisp = xy[-1,:]
print('dx = {0:1.8f} m   | dy = {1:1.8f} m'.format(-z[-8]/1000,z[-7]/1000))
