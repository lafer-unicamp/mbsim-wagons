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

from numpy.matlib import matrix, block, eye, zeros
from numpy.linalg import norm, inv    


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
        self.q0 = matrix([[u1,u2,up1,up2]])
        self.q = 0 * self.q0
        
        
        
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
        return self.q + self.q0
    
    

class material(object):
    '''
    material
    '''
    
    def __init__(self,name):
        
        
        self.name = name
        
class linearElasticMaterial(material):
    '''
    Linear elastic material
    '''
    
    def __init__(self,name,E,nu,rho):
        ''' Material initialization method

        Parameters
        ----------
        name : STRING
            MATERIAL NAME.
        E : DOUBLE
            YOUNG'S MODULUS.
        nu : DOUBLE
            POISSON'S RATIO.
        rho : DOUBLE
            DENSITY.

        Returns
        -------
        None.

        '''
        material.__init__(self,name)
        
        self.E = E
        self.nu = nu
        self.rho = rho
        
        # Lamé constants
        mu = self.E / (1+self.nu) * 0.5
        lam = self.nu * self.E / ( (1+self.nu)*(1-2*self.nu) )
        
        D = eye(6)*0
        
        # following Lai, Introduction to Continuum Mechanics, 2010 notation
        C11 = lam + 2*mu
        D[0,0] = C11
        D[1,1] = C11
        D[2,2] = C11
        
        D[0,1] = lam
        D[1,2] = lam
        D[2,3] = lam
        D[1,0] = lam
        D[2,1] = lam
        D[3,2] = lam
        
        D[3,3] = 2*mu
        D[4,4] = 2*mu
        D[5,5] = 2*mu
        
        self.constitutiveMatrix = D
        
    
    
    
  
class flexibleBody(object):
    '''
    Flexible body class
    '''
    
    def __init__(self,name,material):
        '''
        Flex body initialization method

        Parameters
        ----------
        name : STRING
            NAME OF THE BODY.
        material : MATERIAL
            BODY'S MATERIAL.

        Returns
        -------
        None.

        '''
        
        self.name = name
        self.material = material
        
        
class element(object):
    
    def __init__(self, parentBody):
        self.parentBody = parentBody




class elementLinear(element):
    """
    Planar finite element with linear interpolation
    """
     
    def __init__(self, node1, node2, height, parentBody):
        element.__init__(self,parentBody)
        self.length = norm(node1.qtotal[0:2] - node2.qtotal[0:2])
        self.height = height
        self.nodes = [node1,node2]
        
        L = self.length
        h = self.height
        self.mass = self.parentBody.material.rho * L * h
        self.massMatrix = self.parentBody.material.rho * L * matrix([
            [  h/3,     0,         0,         0,   h/6,     0,         0,         0],
            [    0,   h/3,         0,         0,     0,   h/6,         0,         0],
            [    0,     0,   h**3/36,         0,     0,     0,   h**3/72,         0],
            [    0,     0,         0,   h**3/36,     0,     0,         0,   h**3/72],
            [  h/6,     0,         0,         0,   h/3,     0,         0,         0],
            [    0,   h/6,         0,         0,     0,   h/3,         0,         0],
            [    0,     0,   h**3/72,         0,     0,     0,   h**3/36,         0],
            [    0,     0,         0,   h**3/72,     0,     0,         0,   h**3/36]])
        
        
    def interpolatePosition(self,xi_, eta_):
        """
        Returns the interpolated position given the non-dimensional parameters
        xi_ and eta_. Notice that xi_ and eta_ range from -1 to 1.
        """
        return self.shapeFunctionMatrix(xi_ ,eta_) * block([[x.nodes[0].qtotal.T],[x.nodes[1].qtotal.T]])

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
        
        dSx_dxi =  1/L * matrix([-1,0 ,-eta,0   ,1,0,eta,0])
        dSy_dxi =  1/L * matrix([0 ,-1,0   ,-eta,0,1,0  ,eta])
        dSx_deta = 1/L * matrix([0 ,0 ,S1  ,0   ,0,0,S3 ,0])
        dSy_deta = 1/L * matrix([0 ,0 ,0   ,S1  ,0,0,0  ,S3])
        
        dS = [[dSx_dxi,dSx_deta],[dSy_dxi,dSy_deta]]
        
        return dS[coord][param]
    
    
    
    
    def initialJacobian(self,xi_,eta_):
        q = block([[x.nodes[0].q0.T],[x.nodes[1].q0.T]])
        return self.getJacobian(xi_,eta_,q)
    
    def inverseInitialJacobian(self,xi_,eta_):
        J0 = self.initialJacobian(xi_,eta_)
        return inv(J0)
    
    def currentJacobian(self,xi_,eta_):
        q = block([[x.nodes[0].qtotal.T],[x.nodes[1].qtotal.T]])
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
        invJ0 = self.inverseInitialJacobian(xi_,eta_)
        
        q = block([self.nodes[0].qtotal,self.nodes[1].qtotal]).T
        
        U11 = (self.shapeFunctionDerivative(0,0,xi_,eta_).T*
               self.shapeFunctionDerivative(0,0,xi_,eta_) + 
               self.shapeFunctionDerivative(1,0,xi_,eta_).T*
               self.shapeFunctionDerivative(1,0,xi_,eta_))*q
        U12 = (self.shapeFunctionDerivative(0,0,xi_,eta_).T*
               self.shapeFunctionDerivative(0,1,xi_,eta_) + 
               self.shapeFunctionDerivative(1,0,xi_,eta_).T*
               self.shapeFunctionDerivative(1,1,xi_,eta_))*q
        U21 = (self.shapeFunctionDerivative(0,1,xi_,eta_).T*
               self.shapeFunctionDerivative(0,0,xi_,eta_) + 
               self.shapeFunctionDerivative(1,1,xi_,eta_).T*
               self.shapeFunctionDerivative(1,0,xi_,eta_))*q
        U22 = (self.shapeFunctionDerivative(0,1,xi_,eta_).T*
               self.shapeFunctionDerivative(0,1,xi_,eta_) + 
               self.shapeFunctionDerivative(1,1,xi_,eta_).T*
               self.shapeFunctionDerivative(1,1,xi_,eta_))*q
                     
        U = block([[U11,U12],[U21,U22]])
        print(U)
        return invJ0.T * U * invJ0
    
    
    def getMassMatrix(self):
        
        # Gauss integration points
        gauss = []
        gauss.append(-0.7745967)
        gauss.append(0)
        gauss.append(-gauss[0])
        
        # Gauss weights
        w = []
        w.append(0.555556)
        w.append(0.888889)
        w.append(w[0])
        
        M = 0*eye(8)
        
        for i in range(3):
            for j in range(3):
                S = self.shapeFunctionMatrix(gauss[i],gauss[j])
                M = M + S.T*S * w[i] * w[j]
                
        """we have to multiply by the length and height because
        calculations are carried out on non-dimensional coordinates [-1,1]
        """        
        return self.parentBody.material.rho * M * self.length * self.height / 4
    
    
    def getNodalForces(self):
        # Gauss integration points
        gauss = []
        gauss.append(-0.7745967)
        gauss.append(0)
        gauss.append(-gauss[0])
        
        # Gauss weights
        w = []
        w.append(0.555556)
        w.append(0.888889)
        w.append(w[0])
        
        Qe = zeros([8,1])
        
        # Gaussian quadrature
        for i in range(3):
            for j in range(3):
                print('Not implemented')
        
        
                
           

'''
TEST PROGRAM
'''


steel = linearElasticMaterial('Steel',210e3,0.3,7.85e-6)
body = flexibleBody('Bar',steel)


n1 = node(0,0,0,1)
n2 = node(10,0,0,1)
x = elementLinear(n1, n2, 1,body)
# displace node 2
n2.q = matrix([[0.01,0,0,0]])

#x.strainTensorDerivative(0,0)