# -*- coding: utf-8 -*-
"""
Created on Tue Apr  8 16:14:30 2014

@author: hok1
"""

import paraferromagnet as pfm
import conicalphase as cp
import helicalphase as hp
import itertools
import numpy as np
import matplotlib.pyplot as plt
import csv
from operator import add

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

    def plot_phase_diagram_2d(self, param_names):
        if len(param_names) < 2:
            return
        horiz_param_val = map(lambda point: point[0][param_names[0]], self.points)
        vert_param_val = map(lambda point: point[0][param_names[1]], self.points)
        phase_names = map(lambda point: point[1], self.points)
        phase_name_set = list(set(phase_names))
        phase_name_code_dict = {}
        for idx, phase_name in zip(range(len(phase_name_set)), phase_name_set):
            phase_name_code_dict[phase_name] = idx
        phase_values = np.array(map(lambda name: phase_name_code_dict[name], phase_names))
        plt.scatter(horiz_param_val, vert_param_val, c=phase_values)

    def save_points(self, filename):
        fields = list(set(reduce(add, map(lambda point: point[0].keys(), self.points))))
        field_to_col_dict = {}
        for idx, field in zip(range(len(fields)), fields):
            field_to_col_dict[field] = idx
        outf = open(filename, 'wb')
        writer = csv.writer(outf)
        writer.writerow(fields+['phase'])
        for point in self.points:
            rowtowrite = [0]*len(fields)
            for param in point[0]:
                rowtowrite[field_to_col_dict[param]] = point[0][param]
            rowtowrite.append(point[1])
            writer.writerow(rowtowrite)
        outf.close()