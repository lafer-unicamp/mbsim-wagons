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

from numpy.matlib import matrix, block, eye, zeros
from numpy.linalg import norm, inv, det

#import flexibleBody
#import materials

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
    
    def __init__(self, u1=0, u2=0, up1 = 0, up2 = 1):
        self.q0 = [u1,u2,up1,up2]
        self.q = [0.0]*4
        self.globalDof = [0,1,2,3]
        
            
        
    @property
    def q(self):
        """nodal displacement relative to reference configuration"""
        return self._q
    @q.setter
    def q(self,dofMatrix):
        self._q = dofMatrix
        
    @property
    def qtotal(self):
        """nodal displacement relative to global frame"""
        return np.array(self.q) + np.array(self.q0)
    
    

    
    

        
        
        
        
        
         

    
########################################
class beamANCFelement(object):
    '''
    Base class for finite elements
    '''
    def __init__(self):
        self.parentBody = None
        self.length = 1.0
        self.height = 1.0
        self.width = 1.0
        
       
    
    @property
    def q0(self):
        q0 = []
        for nd in self.nodes:
            q0 = q0 + nd.q0
        return q0
    
    @property
    def gaussIntegrationPoints(self):
        return {1:[0],
                2:[-0.57735027,0.57735027],
                3:[-0.7745967, 0, 0.7745967]}
    
    @property
    def gaussWeights(self): 
        return {1:[2],
                2:[1,1],
                3:[0.555556, 0.888889, 0.555556]}
    
    @property
    def qtotal(self):
        """nodal position relative to global frame"""
        myq = np.empty((0,))
        for node in self.nodes:
            myq = np.append(myq,node.qtotal)  
        return myq
    
    @property
    def q(self):
        '''nodal displacement relative to global frame'''
        myq = np.empty((0,))
        for node in self.nodes:
            myq = np.append(myq,node.q)
            
        return myq
    
    @property
    def globalDof(self):
        gd = []
        for nd in self.nodes:
            gd.extend(nd.globalDof)
        return gd
    
    @property
    def mass(self):
        return self.length * self.height * self.width * self.parentBody.material.rho
    
    
    
    def interpolatePosition(self,xi_, eta_):
        """
        Returns the interpolated position given the non-dimensional parameters
        xi_ and eta_. Notice that xi_ and eta_ range from -1 to 1.
        """
        
        r = np.dot(self.shapeFunctionMatrix(xi_ ,eta_), self.qtotal.T)
        
        return r.reshape(1,-1)
    
    
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
        
        
        
        return matrix([[np.dot(dSx_dxi.A1,q) , np.dot(dSx_deta.A1,q)],
                      [np.dot(dSy_dxi.A1,q) , np.dot(dSy_deta.A1,q)]])
    
    
    
    def initialJacobian(self,xi_,eta_):
        return self.getJacobian(xi_,eta_,self.q0)
    
    def inverseInitialJacobian(self,xi_,eta_):
        J0 = self.initialJacobian(xi_,eta_)
        return inv(J0)
    
    def currentJacobian(self,xi_,eta_):
        return self.getJacobian(xi_,eta_,self.qtotal)
    
    def getMassMatrix(self):
        
        # Gauss integration points
        gauss = self.gaussIntegrationPoints[3]
        npoints = len(gauss)
        
        # Gauss weights
        w = self.gaussWeights[3]
        
        M = 0*eye(len(self.q))
        
        for i in range(npoints):
            for j in range(npoints):
                S = self.shapeFunctionMatrix(gauss[i],gauss[j])
                M = M + S.T*S * w[i] * w[j]
                
        """we have to multiply by the length and height because
        calculations are carried out on non-dimensional coordinates [-1,1]
        """        
        return self.parentBody.material.rho * M * self.length * self.height * self.width / 4
    
    
    def stressTensorByPosition(self,xi_,eta_):
        return self.parentBody.material.stressTensor(self.strainTensor(xi_, eta_))
    
    def strainTensor(self,xi_,eta_):
        F = self.currentJacobian(xi_,eta_) * inv(self.initialJacobian(xi_,eta_))
        
        return 0.5 * (F.T * F - eye(2))
    
    
    def shapeFunctionDerivative(self,xi_,eta_):
        raise NotImplementedError()
    
    
    
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
        
        q = self.qtotal
        
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
        
        ndof = len(q)
        deps_dq = np.zeros((2,2,ndof))
        
        for m in range(ndof):
            deps_dq[:,:,m] = 0.5 * invJ0.T * U[:,(m,m+ndof)] * invJ0
        
        # L = self.length
        # eta = eta_*self.height/2
        # xi = xi_*L/2
        
        # deps_dq[:,:,0] = np.array([[-1.0*(-eta*q[0,2]/L + eta*q[0,6]/L - q[0,0]/L + q[0,4]/L)/L,
        #                       -(0.5*q[0,2]*(L/2 - xi)/L + 0.5*q[0,6]*(L/2 + xi)/L)/L],
        #                      [-(0.5*q[0,2]*(L/2 - xi)/L + 0.5*q[0,6]*(L/2 + xi)/L)/L,
        #                       0]])
        # deps_dq[:,:,1] = np.array([[-1.0*(-eta*q[0,3]/L + eta*q[0,7]/L - q[0,1]/L + q[0,5]/L)/L,
        #                         -(0.5*q[0,3]*(L/2 - xi)/L + 0.5*q[0,7]*(L/2 + xi)/L)/L],
        #                        [-(0.5*q[0,3]*(L/2 - xi)/L + 0.5*q[0,7]*(L/2 + xi)/L)/L, 
        #                         0]])
        # deps_dq[:,:,2] = np.array([[-1.0*eta*(-eta*q[0,2]/L + eta*q[0,6]/L - q[0,0]/L + q[0,4]/L)/L, 
        #                         -eta*(0.5*q[0,2]*(L/2 - xi)/L + 0.5*q[0,6]*(L/2 + xi)/L)/L + 0.5*(L/2 - xi)*(-eta*q[0,2]/L + eta*q[0,6]/L - q[0,0]/L + q[0,4]/L)/L],
        #                        [-eta*(0.5*q[0,2]*(L/2 - xi)/L + 0.5*q[0,6]*(L/2 + xi)/L)/L + 0.5*(L/2 - xi)*(-eta*q[0,2]/L + eta*q[0,6]/L - q[0,0]/L + q[0,4]/L)/L,
        #                         1.0*(L/2 - xi)*(q[0,2]*(L/2 - xi)/L + q[0,6]*(L/2 + xi)/L)/L]])
        # deps_dq[:,:,3] = np.array([[-1.0*eta*(-eta*q[0,3]/L + eta*q[0,7]/L - q[0,1]/L + q[0,5]/L)/L, 
        #                             -eta*(0.5*q[0,3]*(L/2 - xi)/L + 0.5*q[0,7]*(L/2 + xi)/L)/L + 0.5*(L/2 - xi)*(-eta*q[0,3]/L + eta*q[0,7]/L - q[0,1]/L + q[0,5]/L)/L],
        #                            [-eta*(0.5*q[0,3]*(L/2 - xi)/L + 0.5*q[0,7]*(L/2 + xi)/L)/L + 0.5*(L/2 - xi)*(-eta*q[0,3]/L + eta*q[0,7]/L - q[0,1]/L + q[0,5]/L)/L,
        #                             1.0*(L/2 - xi)*(q[0,3]*(L/2 - xi)/L + q[0,7]*(L/2 + xi)/L)/L]])
        # deps_dq[:,:,4] = np.array([[1.0*(-eta*q[0,2]/L + eta*q[0,6]/L - q[0,0]/L + q[0,4]/L)/L,
        #                             (0.5*q[0,2]*(L/2 - xi)/L + 0.5*q[0,6]*(L/2 + xi)/L)/L],
        #                            [(0.5*q[0,2]*(L/2 - xi)/L + 0.5*q[0,6]*(L/2 + xi)/L)/L,
        #                             0]])
        # deps_dq[:,:,5] = np.array([[1.0*(-eta*q[0,3]/L + eta*q[0,7]/L - q[0,1]/L + q[0,5]/L)/L,
        #                             (0.5*q[0,3]*(L/2 - xi)/L + 0.5*q[0,7]*(L/2 + xi)/L)/L],
        #                            [(0.5*q[0,3]*(L/2 - xi)/L + 0.5*q[0,7]*(L/2 + xi)/L)/L,
        #                             0]])
        # deps_dq[:,:,6] = np.array([[1.0*eta*(-eta*q[0,2]/L + eta*q[0,6]/L - q[0,0]/L + q[0,4]/L)/L,
        #                             eta*(0.5*q[0,2]*(L/2 - xi)/L + 0.5*q[0,6]*(L/2 + xi)/L)/L + 0.5*(L/2 + xi)*(-eta*q[0,2]/L + eta*q[0,6]/L - q[0,0]/L + q[0,4]/L)/L],
        #                            [eta*(0.5*q[0,2]*(L/2 - xi)/L + 0.5*q[0,6]*(L/2 + xi)/L)/L + 0.5*(L/2 + xi)*(-eta*q[0,2]/L + eta*q[0,6]/L - q[0,0]/L + q[0,4]/L)/L, 
        #                             1.0*(L/2 + xi)*(q[0,2]*(L/2 - xi)/L + q[0,6]*(L/2 + xi)/L)/L]])
        # deps_dq[:,:,7] = np.array([[1.0*eta*(-eta*q[0,3]/L + eta*q[0,7]/L - q[0,1]/L + q[0,5]/L)/L, 
        #                             eta*(0.5*q[0,3]*(L/2 - xi)/L + 0.5*q[0,7]*(L/2 + xi)/L)/L + 0.5*(L/2 + xi)*(-eta*q[0,3]/L + eta*q[0,7]/L - q[0,1]/L + q[0,5]/L)/L], 
        #                            [eta*(0.5*q[0,3]*(L/2 - xi)/L + 0.5*q[0,7]*(L/2 + xi)/L)/L + 0.5*(L/2 + xi)*(-eta*q[0,3]/L + eta*q[0,7]/L - q[0,1]/L + q[0,5]/L)/L, 
        #                             1.0*(L/2 + xi)*(q[0,3]*(L/2 - xi)/L + q[0,7]*(L/2 + xi)/L)/L]])
        
       
        return deps_dq 
    
    
    
    def getNodalElasticForces(self):
        # Gauss integration points
        gaussL = self.gaussIntegrationPoints[3]
        gaussH = self.gaussIntegrationPoints[2]
        
        # Gauss weights
        wL = self.gaussWeights[3]
        wH = self.gaussWeights[2]
        
        # beam geometry
        L = self.length
        H = self.height
        W = self.width
        
        ndof = len(self.q)
        Qe = zeros([ndof,1])
        
        # regular Gaussian quadrature
        
        # for p in range(npoints):
        #     'length quadrature'
        #     for b in range(npoints):
        #         'heigth quadrature'
        #         detJ0 = det(self.initialJacobian(gaussH[p], gaussH[b]))
        #         deps_dq = self.strainTensorDerivative(gaussH[p], gaussH[b])
        #         T, Tc = self.stressTensor(self.strainTensor(gaussH[p], gaussH[b]),split=True)
        #         for m in range(8):
        #             for i in range(2):
        #                 for j in range(2):
        #                     Qe[m,0] += deps_dq[i,j,m]*T[i,j] * detJ0 * w3[p] * w3[b] * L/2 * H/2
        #                     if b == 1:
        #                         Qe[m,0] += deps_dq[i,j,m]*Tc[i,j] * detJ0 * w3[p] * H * L/2
                                
                                
        # selective reduced integration
        for p in range(len(gaussL)):
            'length quadrature'
            for b in range(len(gaussH)):
                'heigth quadrature'
                detJ0 = det(self.initialJacobian(gaussL[p], gaussH[b]))
                deps_dq = self.strainTensorDerivative(gaussL[p], gaussH[b])
                T, Tc = self.parentBody.material.stressTensor(self.strainTensor(gaussL[p], gaussH[b]),split=True)
                for m in range(ndof):
                    for i in range(2):
                        for j in range(2):
                            Qe[m,0] += deps_dq[i,j,m]*T[i,j] * detJ0 * wL[p] * wH[b]
            # end of height quadrature
            detJ0 = det(self.initialJacobian(gaussL[p], 0))
            deps_dq = self.strainTensorDerivative(gaussL[p], 0)
            T, Tc = self.parentBody.material.stressTensor(self.strainTensor(gaussL[p], 0),split=True)
            for m in range(ndof):
                    for i in range(2):
                        for j in range(2):
                            Qe[m,0] += 2 * deps_dq[i,j,m]*Tc[i,j] * detJ0 * wL[p]
   

        return Qe * W * L * H / 4
    
    def getWeightNodalForces(self,grav):
        L = self.length
        H = self.height
        W = self.width
        Qg =  L * H * W * 0.25 * grav * matrix([
            [2, 0, 0, 0, 2, 0, 0, 0],
            [0, 2, 0, 0, 0, 2, 0, 0]])*eye(len(self.q))*self.parentBody.material.rho
        
        return Qg
    



#%%
class elementLinear(beamANCFelement):
    """
    Planar finite element with linear interpolation
    """
     
    def __init__(self, node1, node2, _height, _width):
        self.length = norm(node1.qtotal[0:2] - node2.qtotal[0:2])
        self.height = _height
        self.width = _width
        self.nodes = [node1,node2]
        
        
    

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
    
    
    
    
    
    
    
class elementQuadratic(beamANCFelement):
    """
    Planar finite element with quadratic interpolation
    """
     
    def __init__(self, node1, node2, _height, _width):
        self.length = norm(node1.qtotal[0:2] - node2.qtotal[0:2])
        self.height = _height
        self.width = _width
        intermediateNode = node()
        intermediateNode.q0 = [(a+b)*0.5 for a,b in zip(node1.q0,node2.q0)]
        self.nodes = [node1,intermediateNode,node2]
        
  
    def shapeFunctionMatrix(self, xi_, eta_):
        '''
        Shape functions respect the order of the nodes: 1, intermediate, 2
        '''
        eta = eta_ * self.height / 2

        S1 = - xi_/2 * (1-xi_)
        S2 = eta * S1
        S3 = xi_/2 * (1+xi_)
        S4 = eta * S3
        S5 = 1 - xi_*xi_
        S6 = eta*S5
        
        return matrix([[S1, 0 ,S2, 0 , S5, 0 , S6, 0 , S3, 0 ,S4 ,0],
                       [0 , S1, 0, S2, 0 , S5, 0 , S6, 0 ,S3 ,0  ,S4]])
    
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
        
        # all the following must be scaled by 1/L^2. We do that in return
        if [coord,param] == [0,0]:
            s1 = -L + 4*xi
            s2 = L + 4*xi
            dS =  [s1, 0, eta*s1, 0, -8*xi, 0, -8*xi*eta, 0, s2, 0, eta*s2, 0]
        elif [coord,param] == [1,0]:
            s1 = -L + 4*xi
            s2 = L + 4*xi
            dS =  [0, s1, 0, eta*s1, 0, -8*xi, 0, -8*eta*xi, 0, s2, 0,  eta*s2]
        elif [coord,param] == [0,1]:
            dS = [0, 0, xi*(-L + 2*xi), 0, 0, 0, L**2 - 4*xi**2, 0, 0, 0, xi*(L + 2*xi), 0]
        else:
            dS = [0 ,0 ,0, xi*(-L + 2*xi),0 ,0 ,0 ,L**2 - 4*xi**2, 0, 0, 0, xi*(L + 2*xi)]
        
        return matrix(dS) * 1/L**2





    

            