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
        
        # Lamé constants
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