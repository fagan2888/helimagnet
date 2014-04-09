# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 16:14:30 2014

@author: hok1
"""

import paraferromagnet as pfm
import conicalphase as cp
import helicalphase as hp

class PhaseDeterminationMachine:
    def __init__(self, configdict, phaseclasses):
        self.configdict = configdict
        self.phaseclasses = phaseclasses
        
    
    
