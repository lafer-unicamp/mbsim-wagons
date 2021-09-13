#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementation of the planar beam in 
K. Nachbagauer, A. Pechstein, H. Irschik, e J. Gerstmayr, “A new locking-free 
formulation for planar, shear deformable, linear and quadratic 
beam finite elements based on the absolute nodal coordinate formulation”, 
Multibody System Dynamics, vol. 26, no 3, p. 245–263, 2011,
doi: 10.1007/s11044-011-9249-8.


Created on Thu May 13 14:49:04 2021

@author: Leonardo Bartalini Baruffaldi
"""

import numpy as np

from numpy.matlib import matrix, block, eye, zeros, square
from numpy.linalg import norm, inv, det
import scipy.optimize as opt

import flexibleBody
import materials

class node(object):
    """
    finite element node with four dof
    
    Parameters
    __________
        u1 - reference coordinate x
        u2 - reference coordinate y
        up1 - reference slope relative to x
        up2 - reference slope relative to y
    """
    
    def __init__(self, u1, u2, up1 = 0, up2 = 0):
        self.q0 = np.array([[u1,u2,up1,up2]])
        self.q = 0 * self.q0
        
        
        
    @property
    def q(self):
        """nodal displacement relative to reference configuration"""
        return self._q
    @q.setter
    def q(self,dofMatrix):
        self._q = dofMatrix.reshape(1,-1)
        
    @property
    def qtotal(self):
        """nodal displacement relative to global frame"""
        return self.q + self.q0
    
    

    
    

        
        
        
        
        
         

    
########################################
class element(object):
    '''
    Generic class for finite elements
    '''
    def __init__(self):
        self.parentBody = None
    




class elementLinear(element):
    """
    Planar finite element with linear interpolation
    """
     
    def __init__(self, node1, node2, _height, _width):
        self.length = norm(node1.qtotal[0:2] - node2.qtotal[0:2])
        self.height = _height
        self.width = _width
        self.nodes = [node1,node2]
        
        self.gaussIntegrationPoints = [-0.7745967, 0, 0.7745967]
        self.gaussWeights = [0.555556, 0.888889, 0.555556]
        
        # L = self.length
        # h = self.height
        # self.mass = self.parentBody.material.rho * L * h
        # self.massMatrix = self.parentBody.material.rho * L * matrix([
        #     [  h/3,     0,         0,         0,   h/6,     0,         0,         0],
        #     [    0,   h/3,         0,         0,     0,   h/6,         0,         0],
        #     [    0,     0,   h**3/36,         0,     0,     0,   h**3/72,         0],
        #     [    0,     0,         0,   h**3/36,     0,     0,         0,   h**3/72],
        #     [  h/6,     0,         0,         0,   h/3,     0,         0,         0],
        #     [    0,   h/6,         0,         0,     0,   h/3,         0,         0],
        #     [    0,     0,   h**3/72,         0,     0,     0,   h**3/36,         0],
        #     [    0,     0,         0,   h**3/72,     0,     0,         0,   h**3/36]])
    
    @property
    def qtotal(self):
        """nodal displacement relative to global frame"""
        myq = np.empty((0,0))
        for node in self.nodes:
            myq = np.append(myq,node.qtotal)
            
        return matrix(myq)
    
    @property
    def q(self):
        myq = np.empty((0,0))
        for node in self.nodes:
            myq = np.append(myq,node.q)
            
        return matrix(myq)
    
    @property
    def q0(self):
        return block([self.nodes[0].q0, self.nodes[1].q0])
    
    @property
    def mass(self):
        return self.length * self.height * self.width * self.parentBody.material.rho
        
        
    def interpolatePosition(self,xi_, eta_):
        """
        Returns the interpolated position given the non-dimensional parameters
        xi_ and eta_. Notice that xi_ and eta_ range from -1 to 1.
        """
        
        r = self.shapeFunctionMatrix(xi_ ,eta_) * block([[self.nodes[0].qtotal.T],[self.nodes[1].qtotal.T]])
        
        return r.reshape(1,-1)

    def shapeFunctionMatrix(self, xi_, eta_):

        L = self.length
        xi = xi_ * L/2
        eta = eta_ * self.height / 2

        S1 = (L/2 - xi)
        S2 = eta * S1
        S3 = (L/2 + xi)
        S4 = eta * S3
        
        return 1/L * matrix([[S1,0 ,S2,0 ,S3,0 ,S4,0],
                       [0 ,S1,0 ,S2,0 ,S3,0 ,S4]])
    
    def shapeFunctionDerivative(self,coord,param,xi_,eta_):
        """

        Parameters
        ----------
        coord : INTEGER
            Number of the coordinate, x-0, y-1, z-2
        param : INTEGER
            Element parameter, xi-0, eta-1, zeta-2
        xi : DOUBLE
            Longitudinal position
        eta : lateral position

        Returns
        -------
        The derivative of the shape function evaluated at interest points.

        """
        
        L = self.length
        xi = xi_ * L/2
        eta = eta_ * self.height / 2
        
        S1 = (L/2 - xi)
        S3 = (L/2 + xi)
        
        # all the following must be scaled by 1/L. We do that in return
        if [coord,param] == [0,0]:
            dS =  [-1,0 ,-eta,0   ,1,0,eta,0]
        elif [coord,param] == [1,0]:
            dS =  [0 ,-1,0   ,-eta,0,1,0  ,eta]
        elif [coord,param] == [0,1]:
            dS = [0 ,0 ,S1  ,0   ,0,0,S3 ,0]
        else:
            dS = [0 ,0 ,0   ,S1  ,0,0,0  ,S3]
        
        return matrix(dS) * 1/L
    
    
    
    
    def initialJacobian(self,xi_,eta_):
        q = block([[self.nodes[0].q0.T],[self.nodes[1].q0.T]])
        return self.getJacobian(xi_,eta_,q)
    
    def inverseInitialJacobian(self,xi_,eta_):
        J0 = self.initialJacobian(xi_,eta_)
        return inv(J0)
    
    def currentJacobian(self,xi_,eta_):
        q = block([[self.nodes[0].qtotal.T],[self.nodes[1].qtotal.T]])
        return self.getJacobian(xi_,eta_,q)
    
    def getJacobian(self,xi_,eta_,q):
        '''
        

        Parameters
        ----------
        xi_ : DOUBLE
            1st ELEMENT INSTRINSIC COORDINATE [-1,1].
        eta_ : DOUBLE
            2nd ELEMENT INSTRINSIC COORDINATE [-1,1].
        q : VECTOR DOUBLE
            NODAL COORDINATES.

        Returns
        -------
        BLOCK MATRIX
            JACOBIAN CALCULATED AT xi_, eta_, under coordinates q.

        '''
  
        dSx_dxi =  self.shapeFunctionDerivative(0,0,xi_,eta_)
        dSy_dxi =  self.shapeFunctionDerivative(1,0,xi_,eta_)
        dSx_deta = self.shapeFunctionDerivative(0,1,xi_,eta_)
        dSy_deta = self.shapeFunctionDerivative(1,1,xi_,eta_)
        
        
        
        return block([[dSx_dxi * q , dSx_deta * q],
                      [dSy_dxi * q , dSy_deta * q]])
    
    def strainTensor(self,xi_,eta_):
        F = self.currentJacobian(xi_,eta_) * inv(self.initialJacobian(xi_,eta_))
        
        return 0.5 * (F.T * F - eye(2))
    
    def strainTensorDerivative(self,xi_,eta_):
        '''
        Gets the strain tensor derivative at a certain point of the element
        (given by xi_ and eta_).

        Parameters
        ----------
        xi_ : DOUBLE
            LONGITUDINAL POSITION [-1,1].
        eta_ : DOUBLE
            TRANVERSAL POSITION [-1,1].

        Returns
        -------
        deps_dq : NUMPY.NDARRAY
            3 rd order tensor of the strain tensor derivative.
            deps_dq[:,:,n] can be used to access the n-th slice of the tensor

        '''
        invJ0 = self.inverseInitialJacobian(xi_,eta_)
        
        q = block([self.nodes[0].qtotal,self.nodes[1].qtotal])
        
        dSx_dxi = self.shapeFunctionDerivative(0,0,xi_,eta_)
        dSx_deta = self.shapeFunctionDerivative(0,1,xi_,eta_)
        dSy_dxi = self.shapeFunctionDerivative(1,0,xi_,eta_)
        dSy_deta = self.shapeFunctionDerivative(1,1,xi_,eta_)
        
        U11 = np.sum([dSx_dxi.T*dSx_dxi, dSy_dxi.T*dSy_dxi],axis=0)
        U12 = np.sum([dSx_dxi.T*dSx_deta, dSy_dxi.T*dSy_deta], axis=0)
        U21 = U12.T
        U22 = np.sum([dSx_deta.T*dSx_deta, dSy_deta.T*dSy_deta],axis=0)
                     
        U = block([[2*np.dot(q,U11),np.dot(q,U12+U21)],
                   [np.dot(q,U21+U12),2*np.dot(q,U22)]])
        
        deps_dq = np.zeros((2,2,8))
        
        for m in range(8):
            deps_dq[:,:,m] = 0.5 * invJ0.T * U[:,(m,m+8)] * invJ0
        
       
        return deps_dq
    
    
    def getMassMatrix(self):
        
        # Gauss integration points
        gauss = self.gaussIntegrationPoints
        
        # Gauss weights
        w = self.gaussWeights
        
        M = 0*eye(8)
        
        for i in range(3):
            for j in range(3):
                S = self.shapeFunctionMatrix(gauss[i],gauss[j])
                M = M + S.T*S * w[i] * w[j]
                
        """we have to multiply by the length and height because
        calculations are carried out on non-dimensional coordinates [-1,1]
        """        
        return self.parentBody.material.rho * M * self.length * self.height * self.width / 4
    
    def stressTensorByPosition(self,xi_,eta_):
        return self.stressTensor(self.strainTensor(xi_, eta_))
    
    def stressTensor(self,strainTensor,split=False):
        '''
        Calcultes the stress tensor given a strain tensor and 
        supposing the material is linear elastic

        Parameters
        ----------
        strainTensor : MATRIX
            2nd order strain tensor.

        Returns
        -------
        2nd order stress tensor
        
        TODO (2021.08.23):
        -------
        This method should be placed inside the MATERIAL class, because it
        depends on the constitutive model. I'll wait, however, until the model
        is full 3D
        
        Reference:
        Lai, W. Michael, David Rubin, e Erhard Krempl. 
        Introduction to continuum mechanics. 
        4th ed. Amsterdam ; Boston: Butterworth-Heinemann/Elsevier, 2010. pp 208

        '''
        
        T = zeros([2,2])
        Tc = T.copy()
        delta = eye(2)
        
        
        # gets Lamè constants from material
        mu = self.parentBody.material.mu
        lam = self.parentBody.material.lam
        
        # this is the trace of the strain tensor
        e = strainTensor.diagonal().sum()
        
        # regular stress tensor
        if not split:
            for i in range(2):
                for j in range(2):
                    T[i,j] = lam*e*delta[i,j] + 2*mu*strainTensor[i,j]
                
                
        # locking-free stress tensor
        else:
            E = self.parentBody.material.E
            nu = self.parentBody.material.nu
            ks = 10*(1+nu)/(12+11*nu)
        
            T[0,0] = E * strainTensor[0,0]
            T[1,1] = E * strainTensor[1,1]
            T[0,1] = 2 * ks * mu * strainTensor[0,1]
            T[1,0] = T[0,1]
            
            Tc[0,0] = nu * ( nu * strainTensor[0,0] + strainTensor[1,1])
            Tc[1,1] = nu * ( strainTensor[0,0] + nu * strainTensor[1,1])
            Tc *=  E / ( 1- nu*nu )
                
                
        return T, Tc
        
        
    
    def getNodalElasticForces(self):
        # Gauss integration points
        gauss = self.gaussIntegrationPoints
        
        # Gauss weights
        w = self.gaussWeights
        
        # beam geometry
        L = self.length
        H = self.height
        W = self.width
        
        Qe = zeros([8,1])
        
        # Gaussian quadrature
        
        for p in range(3):
            'length quadrature'
            for b in range(3):
                'heigth quadrature'
                detJ0 = det(self.initialJacobian(gauss[p], gauss[b]))
                deps_dq = self.strainTensorDerivative(gauss[p], gauss[b])
                T, Tc = self.stressTensor(self.strainTensor(gauss[p], gauss[b]),split=True)
                for m in range(8):
                    for i in range(2):
                        for j in range(2):
                            Qe[m,0] += deps_dq[i,j,m]*T[i,j] * detJ0 * w[p] * w[b]
                            if b == 1:
                                Qe[m,0] += 2 * deps_dq[i,j,m]*Tc[i,j] * detJ0 * w[p]
                            
        
        # TODO: check errors here
        return Qe * W * L * H / 4
    
    
    def getWeightNodalForces(self,grav):
        L = self.length
        H = self.height
        W = self.width
        Qg =  L * H * W * 0.25 * grav * matrix([
            [2, 0, 0, 0, 2, 0, 0, 0],
            [0, 2, 0, 0, 0, 2, 0, 0]])*eye(8)*self.parentBody.material.rho
        
        return Qg
                            
        
       
                



'''
TEST PROGRAM
'''


steel = materials.linearElasticMaterial('Steel',207e9,0.3,7.85e3)
body = flexibleBody.flexibleBody('Bar',steel)


n = []
nel = 1
totalLength = 2.
for i in range(nel+1):
    n.append(node(totalLength * i/nel,0.0,0.0,1.0))


e = []
for j in range(len(n)-1):
    e.append(elementLinear(n[j],n[j+1],0.5,0.1))


g = matrix([[0.0,-9.81]])

body.addElement(e)
body.assembleMassMatrix()


# Constraint matrix 
conDof = [0,1,2,3]
gdl = 4*len(n)
Phi = zeros([len(conDof),gdl])

for d in conDof:
    Phi[d,d] = 1


def f(z):
    
    ndof = (nel + 1)*4
    
    x = z[0:ndof]
    lam = z[ndof:]
    
    dof = 0
    for node in n:
        node.q = np.array(x[dof:dof+4])
        dof += 4
        
    Qe = body.assembleElasticForceVector()
    Qg = body.assembleWeightVector(g)*0
    Qg[-3,0] = 500000. * 0.5 * 0.5 * 0.5
    
    goal = [0]*((nel+1)*4+4)
    
    goal[0:ndof] = (-Qe.T + Qg.T + np.dot(Phi.T,lam)).tolist()[0]
    
    goal[ndof:] = np.dot(Phi,z[0:ndof]).tolist()[0]
    
    print(f"max(|F|) = {np.max(np.abs(goal)):.8e}")
    
    return goal


z0 = [0]*((nel+1)*4+4)
z0[-3] = -500000. * 0.5 * 0.5 * 0.5
z0[-2] = - z0[-3] * 2
'''z = opt.newton_krylov(f,
                  z0,
                  maxiter=40,
                  f_tol=1e-4,
                  verbose=True)'''
z, info, ier, msg = opt.fsolve(f, z0,full_output=True, xtol=1e-12)

xy = body.plotPositions()

tipDisp = xy[-1,1]
lastElemIncl = np.arctan((xy[-1,1]-xy[-2,1])/(xy[-1,0]-xy[-2,0]))





    

            