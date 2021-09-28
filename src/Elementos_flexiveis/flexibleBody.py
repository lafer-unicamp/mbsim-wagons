#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 06:35:10 2021

@author: leonardo
"""

import numpy as np
import matplotlib.pyplot as plt

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
        self.totalDof = 0   # total number of degrees of freedom
        
        print('Created body \'{}\' with material \'{}\''.format(name,material.name))
        
        
    @property 
    def mass(self):
        mass = 0.0
        for e in self.elementList:
            mass += e.mass
            
        return mass
        
    def assembleFEProblem(self):
        self.assembleMassMatrix()
        
    def assembleMassMatrix(self):
        print('Assemblying mass matrix')
         
        M = np.matlib.zeros([self.totalDof, self.totalDof])
        
        for elem in self.elementList:
            m = elem.getMassMatrix()
            for dof in elem.globalDof:
                localDof = elem.globalDof.index(dof)
                M[dof, dof] += m[localDof,localDof]
            
        print('Mass matrix assembly done!')
        return M
    
    
    def assembleElasticForceVector(self):
        
        Qe = np.matlib.zeros(self.totalDof)
        
        
        for elem in self.elementList:
            Qelem = elem.getNodalElasticForces().reshape(1,-1)
            for i,dof in enumerate(elem.globalDof):
                Qe[0,dof] += Qelem[0,i]
            
        return Qe.reshape(-1,1)
    
    def assembleWeightVector(self, g=np.array([0,1])):
        ned = 4             # element degrees of freedom
        nel = len(self.elementList)
        
        Qg = np.matlib.zeros(ned * (nel + 1))
        
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
        curGdl = 0
        for el in element:
            el.parentBody = self
            for nd in el.nodes:
                nd.globalDof = list(range(curGdl,curGdl+4))
                curGdl += 4
            curGdl = el.globalDof[-1]-3 
        self.elementList.extend(element)
        self.totalDof = el.globalDof[-1] + 1
        
    def plotPositions(self, pointsPerElement = 5):
        points = np.linspace(-1.,1.,pointsPerElement)
        
        xy = np.empty([0,2])
        
        for ele in self.elementList:
            for i in range(pointsPerElement):
                xy = np.append(xy,ele.interpolatePosition(points[i],0),axis=0)
                
                
        plt.plot(xy[:,0],xy[:,1])
        
        return xy