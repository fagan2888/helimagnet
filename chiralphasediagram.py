# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 16:14:30 2014

@author: hok1
"""

import paraferromagnet as pfm
import conicalphase as cp
import helicalphase as hp
import itertools

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
        sorted_phase_energies = sorted(phase_energies.items(), key=lambda item: item[1])
        phase_name = sorted_phase_energies[0][0]
        self.points.append((ind_param, phase_name))
        return phase_name
    
    def compute_phase_diagram(self, param_ranges):
        self.points = []
        param_iterator = itertools.product(*param_ranges.values())
        param_names = param_ranges.keys()
        for parameters in param_iterator:
            ind_param = {}
            for param_name, param_value in zip(param_names, parameters):
                ind_param[param_name] = param_value
            self.determine_phase(ind_param)