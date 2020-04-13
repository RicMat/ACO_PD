# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:32:32 2020

@author: darim
"""

import matplotlib.pyplot as plt
import numpy as np

class ACO_PDG():
    
    def __init__(self, width, height, obstacles, cycle_num, ant_num, gamma, evaporation_rate, start, end, brushfire_iter):
        self.width = width
        self.height = height
        self.cycle_num = cycle_num
        self.obstacles = obstacles
        self.ant_num = ant_num
        self.gamma = gamma
        self.pheromone_grid = np.ones(shape=(self.width, self.height))
        self.evaporation_rate = evaporation_rate
        self.start = start
        self.end = end
        self.brushfire_iter = brushfire_iter
        self.heuristic_grid = self.Potential_Field_Heuristic(brushfire_iter)
        self.ants_movements = []
    
    #local pheromone diffusion model
    
    def U_att(self, p, p_e):
        K_a = 1.
        return K_a * np.linalg.norm(p - p_e) ** 2 / 2.

    def U_rep(self, rho, rho_o):
        K_r = 25.
        if rho <= rho_o:
            return K_r * (1/rho - 1/rho_o) ** 2 / 2.
        else:
            return 0

    def U_total(self, p, p_e, rho, rho_o):
        return self.U_att(p, p_e) + self.U_rep(rho, rho_o)
    
    def diff_pheromone(self, phi, r, xi):
        return self.gamma * phi * (r - xi) / r
    
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
        heuristic_grid = np.zeros(shape=(self.width, self.height))
        sdo = self.Brushfire_Algorithm(brushfire_iter)
        for i in range(self.width):
            for j in range(self.height):
                if sdo[i, j] != 0:
                    heuristic_grid[i, j] = self.U_total(np.array([i, j]), self.end, sdo[i, j], 4)
        return 5 * (1 - heuristic_grid / np.max(heuristic_grid))
    
    def check_available_grids(self, position, previous):
        x = position[0]
        y = position[1]
        x_p = previous[0]
        y_p = previous[1]
        grid_list = []
        
        #if x == 5:
            #print("position {}   {}".format(x,y))
            #print("before {}   {}".format(x_p,y_p))
        for i in range(max(x - 1, 0), min(x + 2, 20)):
            for j in range(max(y - 1, 0), min(y + 2, 20)):
                if (not (i, j) == (x_p, y_p)) and (not (i, j) == (x, y)) and (not (i, j) in self.obstacles):
                    if ((i == x-1) and (j == y-1) and ((i+1, j) in self.obstacles)) or \
                    ((i == x-1) and (j == y+1) and ((i, j+1) in self.obstacles)) or \
                    ((i == x+1) and (j == y+1) and ((i-1, j) in self.obstacles)) or \
                    ((i == x+1) and (j == y-1) and ((i, j-1) in self.obstacles)) or \
                    ((i == x-1) and (j == y+1) and ((i-1, j) in self.obstacles)) or \
                    ((i == x-1) and (j == y-1) and ((i, j+1) in self.obstacles)) or \
                    ((i == x+1) and (j == y-1) and ((i+1, j) in self.obstacles)) or \
                    ((i == x+1) and (j == y+1) and ((i, j-1) in self.obstacles)):
                        continue
                        #if x == 5:
                            #print("not ok {}   {}".format(i,j))
                    else:
                        grid_list.append((i, j))
                        #if x == 5:
                            #print("ok {}   {}".format(i,j))
                #else:
                    #if x == 5:
                        #if ((i, j) in self.obstacles):
                            #print("scartata wall {}   {}".format(i,j))
        return grid_list
    
    def algorithm(self, Q): #core
        print("\nComputing ... ")
        for t in range(self.cycle_num):
            #plt.figure()
            #plt.imshow(self.pheromone_grid, cmap=plt.cm.gist_heat)
            #plt.show()
            complete_paths = []
            tmp_pheromone = np.zeros(shape=(self.width, self.height))
            for i in range(self.ant_num):
                path = self.ant_act()
                if path[-1] == (self.end[0], self.end[1]):
                    complete_paths.append((self.path_len(path), path)) #step 7
                    for grid in path: # step 8
                        tmp_pheromone[grid[0], grid[1]] += Q / self.path_len(complete_paths[-1][1])
                        for j in range(max(grid[0] - 1, 0), min(grid[0] + 2, 20)):
                            for k in range(max(grid[1] - 1, 0), min(grid[1] + 2, 20)):
                                if not (j, k) in self.obstacles:
                                    self.pheromone_grid[j, k] += self.diff_pheromone(self.pheromone_grid[j, k], 2, np.linalg.norm(np.array(grid) - np.array((j, k))))
                    path = self.geometry_opt_path(path) #step 9
                    complete_paths.append((self.path_len(path), path))
                    for grid in path: # step 10
                        tmp_pheromone[grid[0], grid[1]] += Q / self.path_len(complete_paths[-1][1])
                        for j in range(max(grid[0] - 1, 0), min(grid[0] + 2, 20)):
                            for k in range(max(grid[1] - 1, 0), min(grid[1] + 2, 20)):
                                if not (j, k) in self.obstacles:
                                    self.pheromone_grid[j, k] += self.diff_pheromone(self.pheromone_grid[j, k], 2, np.linalg.norm(np.array(grid) - np.array((j, k))))
                    
            self.pheromone_grid = (1 - self.evaporation_rate) * self.pheromone_grid  + tmp_pheromone
        
        complete_paths = sorted(complete_paths, key=lambda x: x[0])
        print("\nCompleted")
        return complete_paths
                    
       
    def ant_act(self): #movement - step 3 to 6
        max_steps = 100
        position = self.start
        alpha = 1.1
        beta = 12
        visited_grids = [(position[0],position[1])]
        for i in range(max_steps): #step 6
            available_grids = self.check_available_grids(visited_grids[i], visited_grids[i-1]) #step 5
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
    
    def plot(self, path):
        table = np.ones(shape=(self.width, self.height))
        for obstacle in self.obstacles:
            table[obstacle[0], obstacle[1]] = 0
        for grid in path:
            table[grid[0], grid[1]] = 0.5
        plt.figure()
        plt.imshow(table, cmap=plt.cm.gist_heat)
        plt.show()
        return table
    
    def path_len(self, path):
        length = 0
        for i in range(1, len(path)):
            length += np.linalg.norm(np.array(path[i]) - np.array(path[i-1]))
        return length
    
    def __repr__(self):
        o = str(self.width) + "x" + str(self.height) + " grid with obstacles:\n"
        o += str(self.obstacles)
        return o
        
        

if __name__ == "__main__":

    #main
    
    Obstacles = [(1, 16), 
                 (2, 2), (2, 7), (2, 8), (2, 16), 
                 (3, 1), (3, 2), (3, 7), (3, 15), (3, 16), 
                 (4, 2), (4, 3), (4, 6), (4, 7), (4, 15), (4, 16),
                 (5, 2), (5, 6), (5, 7), 
                 (6, 7), (6, 11),
                 (7, 7), (7, 11), (7, 12),
                 (8, 11), (8, 12),
                 (9, 2), (9, 3), (9, 11), (9, 17),
                 (10, 2), (10, 3), (10, 8), (10, 16), (10, 17),
                 (11, 2), (11, 7), (11, 8), (11, 16), (11, 17),
                 (12, 6), (12, 7), (12, 8),
                 (13, 6), (13, 7), (13, 12), (13, 17),
                 (14, 4), (14, 6), (14, 12), (14, 16), (14, 17), (14, 18),
                 (15, 4), (15, 12), (15, 13), (15, 17),
                 (16, 2), (16, 3), (16, 4), (16, 12), (16, 13), (16, 17),
                 (17, 3), (17, 4), (17, 12), (17, 13),
                 (18, 12)]
    
    aco_pdg = ACO_PDG(20, 20, Obstacles, 25, 10, 0.01, 0.5, np.array([0, 0]), np.array([19, 19]), 5)
    heuristic_grid = aco_pdg.heuristic_grid
    plt.figure()
    plt.imshow(heuristic_grid)
    plt.show()
    print("\nFirst run: ")
    path = aco_pdg.ant_act()
    aco_pdg.plot(path)
    print("\nFirst run after geometry optimization: ")
    path = aco_pdg.geometry_opt_path(path)
    aco_pdg.plot(path)
    o = aco_pdg.algorithm(10)
    aco_pdg.plot(o[0][1])
    print("\nShortest path: " + str(o[0][0]))

