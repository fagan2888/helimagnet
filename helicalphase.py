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
        
    def packedFreeEnergyDensity(self, xvec):
        return self.primitiveFreeEnergyDensity(xvec[0], xvec[1], xvec[2],
                                               xvec[3], xvec[4])
                                               
    def packedFreeEnergyDensityGradient(self, xvec):
        r = self.r
        a = self.a
        c = self.c
        u = self.u
        H = self.H
        a1 = self.a1

        grad0 = 0.5*(2*a*xvec[0]-c)
        grad0 += -2*a1*xvec[0]*xvec[3]*xvec[3]*(2-3*xvec[4])*xvec[4]
        
        grad1 = r*xvec[1]+u*xvec[1]*xvec[1]*xvec[1]-H*math.sqrt(1-2*xvec[4])
        grad1 += u*xvec[1]*(xvec[2]*xvec[2]+xvec[3]*xvec[3])
        
        grad2 = r*xvec[2]+u*xvec[2]*xvec[2]*xvec[2]-H*math.sqrt(2*xvec[4])
        grad2 += u*(xvec[1]*xvec[1]+2*xvec[3]*xvec[3])*xvec[2]
        
        grad3 = (r+a*xvec[0]*xvec[0]-c*xvec[0])*xvec[3]+u*xvec[3]*xvec[3]*xvec[3]
        grad3 += u*(xvec[1]*xvec[1]+2*xvec[2]*xvec[2])*xvec[3]
        grad3 += -2*a1*xvec[0]*xvec[0]*xvec[3]*(2-3*xvec[4])*xvec[4]

        grad4 = H*xvec[1]*xvec[4]/math.sqrt(1-2*xvec[4])
        grad4 += -H*xvec[2]*math.sqrt(0.5*xvec[4])
        grad4 += -2*a1*xvec[0]*xvec[0]*xvec[3]*xvec[3]*(1-3*xvec[4])

        return np.array([grad0, grad1, grad2, grad3, grad4])         
        
    def compute_all(self):
        r = self.r
        a = self.a
        c = self.c
        u = self.u
        H = self.H
        #a1 = self.a1
        
        q0 = 0.5*c/a
        ml0 = H/a/q0/q0
        md0 = 0.0
        msp0 = np.sqrt((a*q0*q0-r)/u)
        betasq0 = self.model_betasq()
        init_x0 = np.array([q0, ml0, md0, msp0, betasq0])
        
        res = minimize(self.packedFreeEnergyDensity, init_x0, method='BFGS',
                       jac=self.packedFreeEnergyDensityGradient)
        self.q = res.x[0]
        self.ml = res.x[1]
        self.md = res.x[2]
        self.m0 = res.x[3]
        self.beta = np.sqrt(res.x[4])
        self.fden = self.packedFreeEnergyDensity(res.x)
            
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
