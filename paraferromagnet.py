# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 16:22:17 2013

@author: hok1
"""

import numpy as np

class FieldPolarizedMagnet:
    def __init__(self, r=0.001, u=0.0001, H=0.0):
        self.recompute(r, u, H)
        
    def recompute(self, r, u, H):
        self.r = r
        self.u = u
        self.H = H
        self.computeML()
        self.valid = True
        
    def computeML(self):
        complex_mls = np.roots([self.u, 0.0, self.r, -self.H])
        mls = np.real(filter(lambda item: item.imag==0, complex_mls))
        fdens = map(self.primitiveFreeEnergyDensity, mls)
        self.ml, self.fden = min(zip(mls, fdens), key=lambda item: item[1])
        
    def primitiveFreeEnergyDensity(self, ml):
        return (0.5*self.r*ml*ml+0.25*self.u*ml*ml*ml*ml-self.H)
        
    def phase_name(self):
        return 'paramagnet' if self.r >= 0 else 'ferromagnet'
