# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 16:14:30 2014

@author: hok1
"""

import paraferromagnet as pfm
import conicalphase as cp
import helicalphase as hp

class PhaseDeterminationMachine:
    def __init__(self, config, phaseclasses):
        self.config = config
        self.phaseclasses = phaseclasses
        self.points = []

    def compute_phase_energies(self, ind_param):
        phase_energies = {}
        for phaseclass in self.phaseclasses:
            phaseclass.recompute_config(self.config, ind_param)
            phase_energies[phaseclass.phase_name()] = phaseclass.fden
        return phase_energies
        
    def determine_phase(self, ind_param):
        phase_energies = self.compute_phase_energies(ind_param)
        sorted_phase_energies = sorted(phase_energies.items(), key=lambda item: item[1], reverse=True)
        return sorted_phase_energies[0][0]
    
