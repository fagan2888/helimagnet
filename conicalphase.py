# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 17:00:31 2013

@author: hok1
"""

import numpy as np

class ConicalPhase:
    def __init__(self, r=0.0001, a=1.0, c=0.025, u=0.1, H=0):
        self.recompute(r, a, c, u, H)
    
    def recompute(self, r, a, c, u, H):
        self.r = r
        self.a = a
        self.c = c
        self.u = u
        self.H = H
        self.computeQ()
        self.computeMH()
        if self.valid:
            self.computeML()
            self.computeMSP()
            self.computeFreeEnergyDensity()
        else:
            self.fden = float("inf")

    def recompute_config(self, config, ind_param):
        self.recompute(ind_param['r'], config['a'], config['c'], config['u'], ind_param['H'])
        
    def computeQ(self):
        self.q = 0.5 * self.c / self.a
        
    def computeMH(self):
        sqMH = (self.a*self.q*self.q - self.r) / self.u
        if sqMH >= 0:
            self.mH = np.sqrt(sqMH)
            self.valid = True
        else:
            self.valid = False
        
    def computeML(self):
        self.ml = self.H / self.a / self.q / self.q
        
    def computeMSP(self):
        self.msp = np.sqrt(self.mH*self.mH - self.ml*self.ml)
        
    def computeFreeEnergyDensity(self):
        if self.valid:
            self.fden = 0.5*self.r*self.mH*self.mH + 0.5*self.a*self.q*self.q*self.msp*self.msp
            self.fden += -0.5*self.c*self.q*self.msp*self.msp+0.25*self.u*self.mH*self.mH*self.mH*self.mH
            self.fden += -self.H*self.ml
        else:
            self.fden = float("inf")
            
    def phase_name(self):
        return 'conical phase'