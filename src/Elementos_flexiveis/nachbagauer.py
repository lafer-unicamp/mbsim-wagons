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
from numpy import dot
from numpy.linalg import norm, inv


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
    
    def __init__(self, u1=0., u2=0., up1 = 0., up2 = 1.):
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
        #self.qtotal = np.array(dofMatrix) + np.array(self.q0)
        
    def updateqTotal(self):
        self.qtotal = np.array(self.q) + np.array(self.q0)
        
    
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
                3:[-0.7745967, 0, 0.7745967],
                4:[-0.86113631,-0.33998104,0.33998104,0.86113631]}
    
    @property
    def gaussWeights(self): 
        return {1:[2],
                2:[1,1],
                3:[0.555556, 0.888889, 0.555556],
                4:[0.33998104,0.6521451549,0.6521451549,0.33998104]}
    
    @property
    def qtotal(self):
        """nodal position relative to global frame"""
        myq = []
        for node in self.nodes:
            myq.extend(node.qtotal.tolist())  
        return myq
    
    @property
    def q(self):
        '''nodal displacement relative to global frame'''
        myq = []
        for node in self.nodes:
            myq.extend(node.q)
            
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
        
        r = dot(self.shapeFunctionMatrix(xi_ ,eta_), self.qtotal)
        
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
        
  
        dSx_dxi, dSx_deta, dSy_dxi, dSy_deta = self.shapeFunctionDerivative(xi_,eta_)
        
        
        
        return matrix([[dot(dSx_dxi,q) , dot(dSx_deta,q)],
                      [dot(dSy_dxi,q) , dot(dSy_deta,q)]])
    
    def saveInitialJacobian(self):
        
        jaco = []       
        
        jaco.append([self.initialJacobian(-1,1),self.initialJacobian(1,1)])
        jaco.append([self.initialJacobian(-1,-1),self.initialJacobian(1,-1)])
        
        invJaco = np.linalg.inv(jaco)
        
        detJaco = np.linalg.det(jaco)
        
        if jaco[0][0].all() == jaco[0][1].all() and jaco[0][1].all() == jaco[1][0].all():
            constant = True
        
        return jaco, invJaco, detJaco, constant
        
    
    def loadInitialJacobian(self,xi_,eta_):
        
        if self.isJ0constant:
            J = self.J0[0][0]
            detJ = self.detJ0[0][0]
            invJ = self.invJ0[0][0]
            
        return J,detJ,invJ

    
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
    
    def getTangentStiffnessMatrix(self):
        '''
        Finite difference approximation to the tangent stiffness matrix

        Returns
        -------
        Kte : Numpy matrix
            Tangent stiffness matrix.

        '''
        
        ndof = len(self.globalDof)
        
        Kte = np.zeros([ndof,ndof])
        
        col = 0
        for nd in self.nodes:
            for i,curDof in enumerate(nd.q):
                savePos = curDof
                Q0 = self.getNodalElasticForces()
                nd.q[i] += 1e-6
                Kte[:,col] = (self.getNodalElasticForces() - Q0) * 1e6
                nd.q[i] = savePos
                col += 1
                
        return Kte
                
    
    
    def stressTensorByPosition(self,xi_,eta_,split=True):
        return self.parentBody.material.stressTensor(self.strainTensor(xi_, eta_),split)
    
    def strainTensor(self,xi_,eta_):
        F = self.currentJacobian(xi_,eta_) * self.loadInitialJacobian(xi_,eta_)[2]
        
        return 0.5 * (dot(F.T,F) - np.matrix([[1.,0.],[0.,1.]]))
    
    
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
        invJ0 = self.loadInitialJacobian(xi_, eta_)[2]
        
        q = self.qtotal
        
        dSx_dxi, dSx_deta, dSy_dxi, dSy_deta = self.shapeFunctionDerivative(xi_,eta_)
        
        
        # TODO: separate into threads
        U11 = np.sum([np.outer(dSx_dxi,dSx_dxi), np.outer(dSy_dxi,dSy_dxi)],axis=0)
        U12 = np.sum([np.outer(dSx_dxi,dSx_deta), np.outer(dSy_dxi,dSy_deta)], axis=0)
        U21 = U12.T
        U22 = np.sum([np.outer(dSx_deta,dSx_deta), np.outer(dSy_deta,dSy_deta)],axis=0)
            
        qU12 = dot(q,U12+U21)       
        U = block([[2*dot(q,U11),qU12],
                   [qU12,2*dot(q,U22)]])
        
        ndof = len(q)
        deps_dq = np.zeros((2,2,ndof))
        
        for m in range(ndof):
            #TODO errors when invJ0 == eye
            #deps_dq[:,:,m] = 0.5 * invJ0.T * U[:,(m,m+ndof)] * invJ0
            deps_dq[:,:,m] = 0.5 * dot(dot(invJ0.T,U[:,(m,m+ndof)]),invJ0)
        
       
        return deps_dq 
    
    
    
    def getNodalElasticForces(self):
        # Gauss integration points
        gaussL = self.gaussIntegrationPoints[3]
        gaussH = self.gaussIntegrationPoints[3]
        
        # Gauss weights
        wL = self.gaussWeights[3]
        wH = self.gaussWeights[3]
        
        # beam geometry
        L = self.length
        H = self.height
        W = self.width
        
        ndof = len(self.q)
        Qe = np.asarray([0.]*ndof)                                
        # selective reduced integration
        for p in range(len(gaussL)):
            'length quadrature'
            for b in range(len(gaussH)):
                'heigth quadrature'
                detJ0 = self.loadInitialJacobian(gaussL[p], gaussH[b])[1]
                deps_dq = self.strainTensorDerivative(gaussL[p], gaussH[b])
                T, Tc = self.parentBody.material.stressTensor(self.strainTensor(gaussL[p], gaussH[b]),split=True)
                Tweight = T * detJ0 * wL[p] * wH[b]
                for m in range(ndof):
                    Qe[m] += np.multiply(deps_dq[:,:,m],Tweight).sum()
            # end of height quadrature
            detJ0 = self.loadInitialJacobian(gaussL[p], 0)[1]
            deps_dq = self.strainTensorDerivative(gaussL[p], 0)
            T, Tc = self.parentBody.material.stressTensor(self.strainTensor(gaussL[p], 0),split=True)
            TcWeight = Tc * detJ0 * wL[p]
            for m in range(ndof):
                Qe[m] += np.multiply(deps_dq[:,:,m],TcWeight).sum()
                
   

        return Qe * W * L * H / 4
    
    
    



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
        self.J0, self.invJ0, self.detJ0, self.isJ0constant = self.saveInitialJacobian()
        
    

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
    
    def shapeFunctionDerivative(self,xi_,eta_):
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

        dSxxi =  np.array([-1,0 ,-eta,0   ,1,0,eta,0])/L

        dSyxi =  np.array([0 ,-1,0   ,-eta,0,1,0  ,eta])/L

        dSxeta = np.array([0 ,0 ,S1  ,0   ,0,0,S3 ,0])/L

        dSyeta = np.array([0 ,0 ,0   ,S1  ,0,0,0  ,S3])/L
        
        return dSxxi,dSxeta, dSyxi, dSyeta
    
    
    def getWeightNodalForces(self,grav):
        L = self.length
        H = self.height
        W = self.width
        Qg =  L * H * W * 0.25 * dot(grav,matrix([
            [2, 0, 0, 0, 2, 0, 0, 0],
            [0, 2, 0, 0, 0, 2, 0, 0]]))*eye(len(self.q))*self.parentBody.material.rho
        
        return Qg.reshape(1,-1)
    
    
    
    
    
#%%    
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
        #intermediateNode.qtotal = np.asarray(intermediateNode.q0) + np.asarray(intermediateNode.q)
        self.nodes = [node1,intermediateNode,node2]
        self.J0, self.invJ0, self.detJ0, self.isJ0constant = self.saveInitialJacobian()
        
  
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
    
    def shapeFunctionDerivative(self,xi_,eta_):
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
        #s1 = -L + 4*xi
        #s2 = L + 4*xi
        s1 = -1 + 2*xi_
        s2 = 1 + 2*xi_
        #dSxxi =  np.array([s1, 0, eta*s1, 0, -8*xi, 0, -8*xi*eta, 0, s2, 0, eta*s2, 0]) * 1/L**2
        dSxxi =  np.array([s1, 0, eta*s1, 0, -4*xi_, 0, -4*xi_*eta, 0, s2, 0, eta*s2, 0]) * 1/L

        dSyxi =  np.array([0, s1, 0, eta*s1, 0, -4*xi_, 0, -4*eta*xi_, 0, s2, 0,  eta*s2]) * 1/L

        #dSxeta = np.array([0, 0, xi*(-L + 2*xi), 0, 0, 0, L**2 - 4*xi**2, 0, 0, 0, xi*(L + 2*xi), 0]) * 1/L**2
        dSxeta = np.array([0, 0, xi_/2*(-1 + xi_), 0, 0, 0, 1-xi_*xi_, 0, 0, 0, xi_/2*(1 + xi_), 0])
        
        dSyeta = np.array([0 ,0 ,0, xi_/2*(-1 + xi_),0 ,0 ,0 ,1-xi_*xi_, 0, 0, 0, xi_/2*(1 + xi_)])
        
        return dSxxi, dSxeta, dSyxi, dSyeta
    
    
    
    def getWeightNodalForces(self,grav):
        '''
        Get nodal forces due to weight
        
        TODO: currently only for initially straight beams

        Parameters
        ----------
        grav : array like
            gravity acceleration.

        Returns
        -------
        Qg : array
            nodal forces.

        '''
        L = self.length
        H = self.height
        W = self.width
        Qg =  L * H * W * 0.25 *  0.3333 * dot(grav, matrix([
            [2, 0, 0, 0, 8, 0, 0, 0, 2, 0, 0, 0],
            [0, 2, 0, 0, 0, 8, 0, 0, 0, 2, 0, 0]]))*eye(len(self.q))*self.parentBody.material.rho
        
        return Qg.reshape(1,-1)





    

            