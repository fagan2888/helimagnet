# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 16:38:09 2013

@author: hok1
"""

import sympy

class HelicalPhase:
    def __init__(self, r=0.0001, a=1.0, c=0.025, u=0.1, H=0, a1=0.01):
        self.recompute(r, a, c, u, H, a1)
    
    def recompute(self, r, a, c, u, H, a1):
        self.r = r
        self.a = a
        self.c = c
        self.u = u
        self.H = H
        self.a1 = a1
    
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
        eqn1, eqn2, eqn3, eqn4, eqn5 = self.get_eqns()
        sols = sympy.solve([eqn1, eqn2, eqn3, eqn4, eqn5], [q, betasq, ml, md, m0])
        return sols
