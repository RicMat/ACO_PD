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
        self.pheromone_grid = np.ones(shape=(self.width, self.height))
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
    
    # algorithm
    
    def Brushfire_Algorithm(self, max_iter): #step 2
        '''Find the shortest distance from an obstacle for each grid using Brushfire algorithm.'''
        sdo = float("inf") * np.ones(shape=(self.width, self.height))
        grid_traverse_mask = np.zeros(shape=(self.width, self.height))
        grid_list = []
        for obstacle in self.obstacles:
            sdo[obstacle[0], obstacle[1]] = 0
            grid_list.append(obstacle)
        for iteration in range(max_iter):
            grid_list_tmp = []
            for grid in grid_list:
                for i in range(max(grid[0] - 1, 0), min(grid[0] + 2, 20)):
                    for j in range(max(grid[1] - 1, 0), min(grid[1] + 2, 20)):
                        if grid_traverse_mask[i, j] == 0:
                            sdo[i, j] = min(sdo[i, j], sdo[grid[0], grid[1]] + \
                               np.linalg.norm([grid[0] - i, grid[1] - j]))
                            grid_list_tmp.append((i, j))
            for grid in grid_list:
                grid_traverse_mask[grid[0], grid[1]] = 1
            grid_list = grid_list_tmp
        return sdo
    
    def Potential_Field_Heuristic(self, brushfire_iter): #step 2
        '''Heuristic value over obstacles are not valid and should not be used.'''
        heuristic_grid = np.zeros(shape=(self.width, self.height))
        sdo = self.Brushfire_Algorithm(brushfire_iter)
        for i in range(self.width):
            for j in range(self.height):
                if sdo[i, j] != 0:
                    heuristic_grid[i, j] = self.U_tot(np.array([i, j]), self.end, sdo[i, j], 4)
        return 5 * (1 - heuristic_grid / np.max(heuristic_grid))
    
    def algorithm(self): #core
        d_pheromone = np.zeros(shape=(self.width, self.height)) #step 3
        
    def ant_act(self): #movement - step 4 to 6
        max_steps = 100
        position = self.start
        alpha = 12
        beta = 0.02
        visited_grids = [(position[0],position[1])]
        for i in range(max_steps): #step 6
            available_grids = check_available_grids(visited_grids[i], visited_grids[i-1]) #step 5
            prob = [] #step 4
            for avail in available_grids:
                prob.append(self.pheromone_grid[avail[0], avail[1]] ** alpha * self.heuristic_grid[avail[0], avail[1]] ** beta)
            prob = np.array(prob)
            prob /= np.sum(prob)
            prob = np.cumsum(prob) # sum them because use a random value, e.g. 0.99, so summing u r sure 
            r = np.random.random() # that at least one, the last one, is bigger than any random value
            for best_j in range(len(prob)): # otherwise would take always the same route
                if prob[best_j] >= r:
                    break
            visited_grids.append(available_grids[best_j])
            if available_grids[best_j][0] == self.end[0] and available_grids[best_j][1] == self.end[1]:
                break
        return visited_grids
        
            
        

if __name__ == "__main__":

    #main


