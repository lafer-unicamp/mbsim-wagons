#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 06:39:30 2021

@author: leonardo
"""

import numpy as np

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
        
        # Lamé constants for 3d stress-strain state
        self.mu = self.E / (1+self.nu) * 0.5
        self.lam = self.nu * self.E / ( (1+self.nu)*(1-2*self.nu) )
        
        D = np.matlib.eye(6)*0
        
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
        
        T = np.zeros([2,2])
        Tc = T.copy()
        delta = np.eye(2)
        
        
        # gets Lamè constants from material
        mu = self.mu
        lam = self.lam
        
        # this is the trace of the strain tensor
        e = strainTensor[0,0] + strainTensor[1,1]
        
        # regular stress tensor
        if not split:
            for i in range(2):
                for j in range(2):
                    #T[i,j] = lam*e*delta[i,j] + 2*mu*strainTensor[i,j]
                    T[0,0] = strainTensor[0,0] + self.nu*strainTensor[1,1]
                    T[1,1] = self.nu*strainTensor[0,0] + strainTensor[1,1] 
                    T[1,0] = (1-self.nu) * strainTensor[1,0]
                    T[0,1] = T[1,0]
                    T = self.E / (1-self.nu**2) * T
                
                
        # locking-free stress tensor
        else:
            E = self.E
            nu = self.nu
            ks = 10*(1+nu)/(12+11*nu)
        
            T[0,0] = E * strainTensor[0,0]
            T[1,1] = E * strainTensor[1,1]
            T[0,1] = 2 * ks * mu * strainTensor[0,1]
            T[1,0] = T[0,1]
            
            Tc[0,0] = nu * ( nu * strainTensor[0,0] + strainTensor[1,1])
            Tc[1,1] = nu * ( strainTensor[0,0] + nu * strainTensor[1,1])
            Tc *=  E / ( 1- nu*nu )
                
                
        return T, Tc