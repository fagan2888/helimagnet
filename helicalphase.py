# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:38:09 2013

@author: hok1
"""

import math
import numpy as np
from scipy.optimize import minimize

class HelicalPhase:
    def __init__(self, r=-0.0001, a=1.0, c=0.025, u=0.1, H=0, a1=0.01):
        self.recompute(r, a, c, u, H, a1)
    
    def recompute(self, r, a, c, u, H, a1):
        self.r = r
        self.a = a
        self.c = c
        self.u = u
        self.H = H
        self.a1 = a1
        self.compute_all()

    def recompute_config(self, config, ind_param):
        self.recompute(ind_param['r'], config['a'], config['c'], config['u'], ind_param['H'], config['a1'])
            
    def model_betasq(self):
        r = self.r
        a = self.a
        c = self.c
        u = self.u
        H = self.H
        a1 = self.a1
        m0sq = 1/u * (0.5*c*c/a-r-16*a*a*H*H*u/(c*c*c*c))
        mlsq = 16*a*a*H*H/(c*c*c*c)
        betasq = (1-0.5*a/a1*mlsq/m0sq) / 3
        return betasq
        
    def primitiveFreeEnergyDensity(self, q, ml, md, m0, betasq):
        if betasq<0:
            return float('inf')
        r = self.r
        a = self.a
        c = self.c
        u = self.u
        H = self.H
        a1 = self.a1
        fe = 0.5*(r+a*q*q-c*q)*m0*m0+0.25*u*m0*m0*m0*m0
        fe += 0.5*r*ml*ml+0.25*u*ml*ml*ml*ml-H*ml*math.sqrt(1-2*betasq)
        fe += 0.5*r*md*md+0.25*u*md*md*md*md-H*md*math.sqrt(2*betasq)
        fe += 0.5*u*(ml*ml*md*md+ml*ml*m0*m0+2*m0*m0*md*md)
        fe += -a1*q*q*m0*m0*(2-3*betasq)*betasq
        return fe
        
    def primitiveFreeEnergyDensityGradient(self, q, ml, md, m0, betasq):
        if betasq<0:
            return np.zeros(5)
        r = self.r
        a = self.a
        c = self.c
        u = self.u
        H = self.H
        a1 = self.a1

        grad0 = 0.5*(2*a*q-c)
        grad0 += -2*a1*q*m0*m0*(2-3*betasq)*betasq
        
        grad1 = r*ml+u*ml*ml*ml-H*math.sqrt(1-2*betasq)
        grad1 += u*ml*(md*md+m0*m0)
        
        grad2 = r*md+u*md*md*md-H*math.sqrt(2*betasq)
        grad2 += u*(ml*ml+2*m0*m0)*md
        
        grad3 = (r+a*q*q-c*q)*m0+u*m0*m0*m0
        grad3 += u*(ml*ml+2*md*md)*m0
        grad3 += -2*a1*q*q*m0*(2-3*betasq)*betasq

        grad4 = H*ml*betasq/math.sqrt(1-2*betasq)
        grad4 += -H*md*math.sqrt(0.5*betasq)
        grad4 += -2*a1*q*q*m0*m0*(1-3*betasq)

        return np.array([grad0, grad1, grad2, grad3, grad4])               
        
    def packedFreeEnergyDensity(self, xvec):
        return self.primitiveFreeEnergyDensity(xvec[0], xvec[1], xvec[2],
                                               xvec[3], xvec[4])
                                               
    def packedFreeEnergyDensityGradient(self, xvec):
        return self.primitiveFreeEnergyDensityGradient(xvec[0], xvec[1], 
                                                       xvec[2], xvec[3],
                                                       xvec[4])     
        
    def compute_all(self):
        r = self.r
        a = self.a
        c = self.c
        u = self.u
        H = self.H
        #a1 = self.a1
        
        q0 = 0.5*c/a
        ml0 = H/a/q0/q0/1.4142
        md0 = H/a/q0/q0/1.4142
        msp0 = np.sqrt((a*q0*q0-r)/u)
        betasq0 = self.model_betasq()
        init_x0 = np.array([q0, ml0, md0, msp0, betasq0])
        
        res = minimize(self.packedFreeEnergyDensity, init_x0,
                       method='Newton-CG',
                       jac=self.packedFreeEnergyDensityGradient)
        self.q = res.x[0]
        self.ml = res.x[1]
        self.md = res.x[2]
        self.m0 = res.x[3]
        if self.q>0 and self.m0>0 and self.md>0 and res.x[4]>0:
            self.valid = True
            self.beta = np.sqrt(res.x[4])
            self.fden = self.packedFreeEnergyDensity(res.x)
        else:
            self.valid = False
            self.beta = None
            self.fden = float('inf')
            
    def computeQ(self):
        pass
        
    def computeMSP(self):
        pass
        
    def computeML(self):
        pass
        
    def computeMD(self):
        pass
        
    def computeBeta(self):
        pass
        
    def computeFreeEnergyDensity(self):
        pass
    
    def phase_name(self):
        return 'helical phase'