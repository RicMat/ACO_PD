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
    
    #local pheromone diffusion model
    
    def U_att(self, K_a, p, p_e):
        return 1/2. * K_a * (p-p_e)^2
    
    def U_rep(self, K_r, rho, rho_0):
        if rho <= rho_0:
            return 0
        else:
            return 1/2. * K_r * (1/rho - 1/rho_0)^2
        
    def U_total(self, u_att, u_rep):
        return u_att + u_rep
    
    def diff_pheromone(self, gamma, phi, r, xi):
        return gamma * phi * (r - xi) / r
    
    # geometry optimization path model
    
    def geometry_opt_path(self, old_path):
        new_path = [old_path[0]] #starting point
        i = 1
        
        while i < len(old_path) - 1:
            current_pos=old_path[i] #step 1
            new_path.append(current_pos)
            grid_subset=[]
            for j in range(max(current_pos[0] - 1, 0), min(current_pos[0] + 2, self.width)): 
                for k in range(max(current_pos[1] - 1, 0), min(current_pos[1] + 2, self.height)):
                    if not (j, k) in self.obstacles: #step 2 remove obstacles
                        grid_subset.append((j, k))
            grid_subset.remove(current_pos) #step 2 remove taboo
            grid_subset.remove(old_path[i-1])
            
            signal=0
            for tmp_pos in grid_subset: #step 3
                if tmp_pos in old_path[i:]:
                    tmp_signal = old_path[i:].index(tmp_pos)
                    if tmp_signal > signal:
                        signal = tmp_signal
            
            i = i + signal #step 5
            
        new_path.append(old_path[-1]) #add goal
        return new_path
                
            
        

if __name__ == "__main__":

    #main


