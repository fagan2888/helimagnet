# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:38:09 2013

@author: hok1
"""

class HelicalPhase:
    def __init__(self, r=0.0001, a=1.0, c=0.025, u=0.1, H=0, a1=0.01, cb=0.0):
        self.recompute(r, a, c, u, H, a1, cb)
    
    def recompute(self, r, a, c, u, H, a1, cb):
        self.r = r
        self.a = a
        self.c = c
        self.u = u
        self.H = H
        self.a1 = a1
        self.cb = cb
        
