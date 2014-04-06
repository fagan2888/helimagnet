# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:38:09 2013

@author: hok1
"""

import sympy
from sympy.functions import im
import math

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
    
    def get_eqns(self):    
        r = self.r
        a = self.a
        c = self.c
        u = self.u
        H = self.H
        a1 = self.a1
        q, betasq, ml, md, m0 = sympy.symbols('q,betasq,ml,md,m0')
        eqn1 = (r+(a-a1*betasq*(2-3*betasq))*q**2-c*q) + u*(m0**2+ml**2+2*md**2)
        eqn2 = (r*ml+u*(ml**2+md**2+m0**2)*ml-H*sympy.sqrt(1-2*betasq))
        eqn3 = r*md+u*(ml**2+md**2)*md+2*u*m0**2*md-H*sympy.sqrt(2*betasq)
        eqn4 = 0.5*(2*a*q-c)-2*a1*q*betasq*(2-3*betasq)
        eqn5 = H*ml/sympy.sqrt(1-2*betasq)-H*md/sympy.sqrt(2*betasq)-2*a1*q**2*m0**2*(1-3*betasq)
        return eqn1, eqn2, eqn3, eqn4, eqn5
        
    def solve_eqn(self):
        q, betasq, ml, md, m0 = sympy.symbols('q,betasq,ml,md,m0')
        nbetasq = self.model_betasq()
        eqn1, eqn2, eqn3, eqn4, eqn5 = self.get_eqns()
        sols = sympy.solve([eqn1.subs(betasq, nbetasq), 
                            eqn2.subs(betasq, nbetasq), 
                            eqn3.subs(betasq, nbetasq),
                            eqn4.subs(betasq, nbetasq)],
                           [q, ml, md, m0])
        isAllReal = lambda item: im(item[0])==0 and im(item[1])==0 and im(item[2])==0 and im(item[3])==0
        sols = filter(isAllReal, sols)
        isValid = lambda item: item[0]>0 and item[1]>=0 and item[2]>=0 and item[3]>0
        sols = filter(isValid, sols)
        return sols
        
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
        
    def compute_all(self):
        betasq = self.model_betasq()
        sols = self.solve_eqn()
        if len(sols) == 0:
            self.valid = False
        else:
            fed = lambda sol: self.primitiveFreeEnergyDensity(sol[0], sol[1],
                                                              sol[2], sol[3],
                                                              betasq)
            fenergies = map(fed, sols)
            sol_fe_pairs = zip(sols, fenergies)
            sol_fe_pairs = sorted(sol_fe_pairs, key=lambda item: item[1])
            sol = sol_fe_pairs[0][0]
            self.q = sol[0]
            self.ml = sol[1]
            self.md = sol[2]
            self.m0 = sol[3]
            self.beta = math.sqrt(betasq)
            self.fe = sol_fe_pairs[0][1]
            self.valid = True
            
    def computeQ(self):
        return float(self.q)
        
    def computeMSP(self):
        return float(self.m0)
        
    def computeML(self):
        return float(self.ml)
        
    def computeMD(self):
        return float(self.md)
        
    def computeBeta(self):
        return float(self.beta)
        
    def computeFreeEnergyDensity(self):
        return float(self.fe)
