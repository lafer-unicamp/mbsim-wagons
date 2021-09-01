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
from scipy.optimize import newton_krylov, anderson, fsolve, broyden1
import matplotlib.pyplot as plt

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
        self._q = dofMatrix.reshape(1,-1)
        
    @property
    def qtotal(self):
        """nodal displacement relative to global frame"""
        return self.q + self.q0
    
    
#################################3
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
        self.mu = self.E / (1+self.nu) * 0.5
        self.lam = self.nu * self.E / ( (1+self.nu)*(1-2*self.nu) )
        
        D = eye(6)*0
        
        # following Lai, Introduction to Continuum Mechanics, 2010 notation
        C11 = self.lam + 2*self.mu
        D[0,0] = C11
        D[1,1] = C11
        D[2,2] = C11
        
        D[0,1] = self.lam
        D[1,2] = self.lam
        D[2,3] = self.lam
        D[1,0] = self.lam
        D[2,1] = self.lam
        D[3,2] = self.lam
        
        D[3,3] = 2*self.mu
        D[4,4] = 2*self.mu
        D[5,5] = 2*self.mu
        
        self.constitutiveMatrix = D
        
    
    
    
########################################
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
        self.elementList = []
        
        
    @property 
    def mass(self):
        mass = 0.0
        for e in self.elementList:
            mass += e.mass
            
        return mass
        
    def assembleFEProblem(self):
        self.assembleMassMatrix()
        
    def assembleMassMatrix(self):
        
        ned = 4             # element degrees of freedom
        nel = len(self.elementList)
        
        M = zeros([ned * (nel + 1), ned * (nel + 1)])
        
        i = 0
        for elem in self.elementList:
            dofstart = ned * i
            dofend = dofstart + 2 * ned
            M[dofstart:dofend, dofstart:dofend] += elem.getMassMatrix()
            i += 1
            
        return M
    
    
    def assembleElasticForceVector(self):
        ned = 4             # element degrees of freedom
        nel = len(self.elementList)
        
        Qe = zeros(ned * (nel + 1))
        
        i = 0
        for elem in self.elementList:
            dofstart = ned * i
            dofend = dofstart + 2 * ned
            Qe[0,dofstart:dofend] += elem.getNodalElasticForces().reshape(1,-1)
            i += 1
            
        return Qe.reshape(-1,1)
    
    def assembleWeightVector(self, g):
        ned = 4             # element degrees of freedom
        nel = len(self.elementList)
        
        Qg = zeros(ned * (nel + 1))
        
        i = 0
        for elem in self.elementList:
            dofstart = ned * i
            dofend = dofstart + 2 * ned
            Qg[0,dofstart:dofend] += elem.getWeightNodalForces(g).reshape(1,-1)
            i += 1
            
        return Qg.reshape(-1,1)
        
    
    def addElement(self, element):
        '''
        

        Parameters
        ----------
        element : ELEMENT
            Element object to be added to elementList.

        Returns
        -------
        None.

        '''
        for el in element:
            el.parentBody = self
        self.elementList.extend(element)
        
    def plot(self, pointsPerElement = 5):
        points = np.linspace(-1.,1.,pointsPerElement)
        
        xy = np.empty([0,2])
        
        for ele in self.elementList:
            for i in range(pointsPerElement):
                xy = np.append(xy,ele.interpolatePosition(points[i],0),axis=0)
                
                
        plt.plot(xy[:,0],xy[:,1],'-o')
        
        return xy
        
        
        
        
        
         

    
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
     
    def __init__(self, node1, node2, height):
        self.length = norm(node1.qtotal[0:2] - node2.qtotal[0:2])
        self.height = height
        self.width = 1.0
        self.nodes = [node1,node2]
        
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
        
        dSx_dxi =  1/L * matrix([-1,0 ,-eta,0   ,1,0,eta,0])
        dSy_dxi =  1/L * matrix([0 ,-1,0   ,-eta,0,1,0  ,eta])
        dSx_deta = 1/L * matrix([0 ,0 ,S1  ,0   ,0,0,S3 ,0])
        dSy_deta = 1/L * matrix([0 ,0 ,0   ,S1  ,0,0,0  ,S3])
        
        dS = [[dSx_dxi,dSx_deta],[dSy_dxi,dSy_deta]]
        
        return dS[coord][param]
    
    
    
    
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
        
        U11 = (self.shapeFunctionDerivative(0,0,xi_,eta_).T*
               self.shapeFunctionDerivative(0,0,xi_,eta_) + 
               self.shapeFunctionDerivative(1,0,xi_,eta_).T*
               self.shapeFunctionDerivative(1,0,xi_,eta_))
        U12 = (self.shapeFunctionDerivative(0,0,xi_,eta_).T*
               self.shapeFunctionDerivative(0,1,xi_,eta_) + 
               self.shapeFunctionDerivative(1,0,xi_,eta_).T*
               self.shapeFunctionDerivative(1,1,xi_,eta_))
        U21 = U12.T
        U22 = (self.shapeFunctionDerivative(0,1,xi_,eta_).T*
               self.shapeFunctionDerivative(0,1,xi_,eta_) + 
               self.shapeFunctionDerivative(1,1,xi_,eta_).T*
               self.shapeFunctionDerivative(1,1,xi_,eta_))
                     
        U = block([[q*U11,q*U12],[q*U21,q*U22]])
        
        deps_dq = np.zeros((2,2,8))
        
        for m in range(8):
            deps_dq[:,:,m] = invJ0.T * U[:,(m,m+8)] * invJ0
        
       
        return deps_dq
    
    
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
    
    def stressTensorByPosition(self,xi_,eta_):
        return self.stressTensor(self.strainTensor(xi_, eta_))
    
    def stressTensor(self,strainTensor):
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
        delta = eye(2)
        
        # gets Lamè constants from material
        mu = self.parentBody.material.mu
        lam = self.parentBody.material.lam
        
        # this is the trace of the strain tensor
        e = strainTensor.diagonal().sum()
        
        # regular stress tensor
        for i in range(2):
            for j in range(2):
                T[i,j] = lam*e*delta[i,j] + 2*mu*strainTensor[i,j]
                
                
        # TODO locking-free stress tensor
                
                
        return T
        
        
    
    def getNodalElasticForces(self):
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
        
        for p in range(3):
            'outer quadrature'
            for b in range(3):
                'inner quadrature'
                detJ0 = det(self.initialJacobian(gauss[p], gauss[b]))
                deps_dq = self.strainTensorDerivative(gauss[p], gauss[b])
                T = self.stressTensor(self.strainTensor(gauss[p], gauss[b]))
                for i in range(2):
                    for j in range(2):
                        'tensor product'
                        for m in range(8):
                            Qe[m,0] += deps_dq[i,j,m]*T[i,j] * detJ0 * w[p] * w[b]
                            
        
        
        return Qe * self.length * self.height / 4
    
    
    def getWeightNodalForces(self,grav):
        L = self.length
        H = self.height
        Qg =  L * H * 0.25 * grav * matrix([
            [2, 0, 0, 0, 2, 0, 0, 0],
            [0, 2, 0, 0, 0, 2, 0, 0]])*eye(8)*self.parentBody.material.rho
        
        return Qg.T

                            
        
       
                



'''
TEST PROGRAM
'''


steel = linearElasticMaterial('Steel',210e3,0.00,7.85e-6)
body = flexibleBody('Bar',steel)


n = []
nel = 1
for i in range(nel+1):
    n.append(node(10 * i/nel,0.0,0.0,1.0))


e = []
for j in range(len(n)-1):
    e.append(elementLinear(n[j],n[j+1],1))


g = matrix([[9.81,0.0]])

body.addElement(e)
body.assembleMassMatrix()
Qe = body.assembleElasticForceVector()
Qg = body.assembleWeightVector(g)
# displace node 2


# def f(z):
    
#     z = matrix(z)
#     z.reshape(1,-1)
    
#     gdl = 0
#     for node in n:
#         node.q = matrix(z[0,gdl:gdl+4])
#         gdl += 4 

    
#     #bcs
#     Qg = body.assembleWeightVector(g)   
#     Qe = body.assembleElasticForceVector()
    
#     conDof = [0,1,2,3]
    
#     Phi = zeros([len(conDof),gdl])
    
#     for d in conDof:
#         Phi[d,d] = 1
    
#     Qbc = Phi.T * z[0,gdl:].T
    
#     f = block([[Qe + Qbc + 1000*Qg],[zeros([4,1])]])
    
#     return f

dt = 5e-5
tend = 0.02

#bcs 
conDof = [0,1,2,3]
gdl = 4*len(n)
Phi = zeros([len(conDof),gdl])

for d in conDof:
    Phi[d,d] = 1

t = []
tnow = 0

# augmented mass matrix

Ma = block([[body.assembleMassMatrix(),-Phi.T],[Phi,zeros([4,4])]])
Mainv = inv(Ma)

z = [np.zeros([4*(nel+1)])]
lam = [np.zeros([4])]
zdot = z.copy()
t = [0.]

while tnow <= tend:
    
    print('t = {}'.format(tnow))
    
    gdl = 0
    for node in n:
        node.q = z[-1][gdl:gdl+4]
        gdl += 4
    
    Qeg = body.assembleElasticForceVector()
    Qgg = body.assembleWeightVector(g)
    
    Qgg = 0.0 * Qg
    Qgg[-3,0] = 1
    
    rhs = np.block([[-Qeg+Qgg],[zeros([4,1])]])
    
    zkdt2 = np.block([z[-1],zeros([1,4])]) / dt / dt
    zdkdt = block([zdot[-1],zeros([1,4])]) / dt
    
    # NOTE: np.dot is faster than regular multiplication
    lhs = np.sum([np.dot(Mainv,rhs),zkdt2.T,zdkdt.T],axis=0)
    
    lhs = matrix(lhs)
    
    z.append(lhs.A1[0:4*(nel+1)] * dt * dt)
    zdot.append((z[-1]-z[-2])/dt)
    lam.append(lhs.A1[4*(nel+1):])
   
    tnow += dt
    t.append(tnow)
    
    

# z0 = newton_krylov(f, zeros([1,(nel+1)*4 + 4]), maxiter=10,verbose=True)

xy = body.plot()      

z = np.array(z)
lam = np.array(lam)
plt.figure()
plt.subplot(2,1,1)
plt.plot(t,z[:,4])
plt.subplot(2,1,2)
plt.plot(t,lam[:,0])



    

            