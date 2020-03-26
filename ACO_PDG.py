# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:32:32 2020

@author: darim
"""

import matplotlib.pyplot as plt
import numpy as np

class ACO_PDG():
    
    def __init__(self, width, height, obstacles, ant_count, step_count, alpha, beta, gamma, evaporation_rate, start, end):
        self.width = width
        self.height = height
        self.obstacles = obstacles
        self.ant_count = ant_count
        self.step_count = step_count
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.evaporation_rate = evaporation_rate
        self.start = start
        self.end = end
    
    def U_att(self, K_a, p, p_e):
        return 1/2. * K_a * (p-p_e)^2
    
    def U_rep(self, K_r, rho, rho_0):
        if rho <= rho_0:
            return 0
        else
            return 1/2. * K_r * (1/rho - 1/rho_0)^2
        
    def U_total(self, u_att, u_rep):
        return u_att + u_rep


if __name__ == "__main__":

    #main


