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
nel = 4
totalLength = 2000.
for i in range(nel+1):
    n.append(node(totalLength * i/nel,0.0,0.0,1.0))
    nq.append(node(totalLength * i/nel,0.0,0.0,1.0))


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

#Kt = body.assembleTangentStiffnessMatrix()



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
    
    simBody.updateDisplacements(x)
            
    Qe = simBody.assembleElasticForceVector()
    Qa = Qe*0
    Qa[-3,0] = -5.0e5 * 0.5 * 0.5 * 0.5
    
    goal = [0]*(gdl+4)
    
    goal[0:ndof] = (-Qe.T + Qa.T - np.dot(Phi.T,lam)).tolist()[0]
    
    goal[ndof:] = np.dot(Phi,x).tolist()
    
    #print(f"max(|F|) = {np.max(np.abs(goal)):.8e}")
    
    record.append(goal)
    
    return goal

z0 = [0]*(gdl+4)
z0[-3] =  5.0e5 * 0.5 * 0.5 * 0.5
z0[-1] = - z0[-3] * 2000
#z = opt.newton_krylov(f,z0,maxiter=40,f_tol=1e-4,verbose=True)

z, info, ier, msg = fsolve(f, z0,full_output=True, col_deriv=True)
print(msg)

xy = simBody.plotPositions(show=False)
tipDisp = xy[-1,:]
gam = np.pi/2-np.arctan2(z[gdl-1]+1,z[gdl-2])
print('dx = {0:1.8e} m   | dy = {1:1.8e} m | theta = {2:1.8e}'.format(-z[-8]/1000,z[-7]/1000,gam))
