#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 19:41:04 2021

@author: dawudabd-alghani
"""

import numpy as np
import pylab as pl
from Ball import Ball
import matplotlib.pyplot as plt

class Simulation():
    
    '''
    Class for simulating how a system of balls in a circular container changes with time (2D).
    The container size and ball postions and velocities are given, or the number of 
    balls is given and that many balls are created that are systematically placed
    inside the container with random velocity that have an average of zero.
    next_collision function finds and performs the next collision. The balls not
    involved in the collision are moved forward by the time jump to the collision.
    run function runs the simulation by performing the next collision repeatedly
    for how many times asked by the user. There is an option in the function to 
    animate the collisions (default is without the animation). If requested specific
    data can be collected and plotted, printed or returned (e.g. pressure, temperature etc).
    '''
    
    # defining time, momentum and distance variables
    # so that they can be used in next_collision and run functions
    t = 0
    p = 0
    distances = []
            
    def __init__(self, container=[Ball(1e20, -100)], balls=[], n=0):
        # defining the container and array of balls as an attribute
        self.__container = container
        self.__balls = balls
    
    #  defining repr function for a simulation
    def __repr__(self):
        return 'Simulation(container=%s,balls=%s)' % (self.__container,self.__balls)
    
    # uses time_to_collision function in Ball script to find the next collision
    # implements the collision using the collide function in the Ball script
    # uses the move function in the Ball script to move all balls that did
    # not collide by the time taken for the collision to occur
    def next_collision(self):
        min_time = np.inf
        
        # double for loop to find the next collison
        for i in range(len(self.__balls)):
            time = self.__balls[i].time_to_collision(self.__container)
            if time < min_time:
                min_time = time
                a = self.__balls[i]
                b = self.__container
            for j in range(len(self.__balls)):
                if i != j:
                    time = self.__balls[i].time_to_collision(self.__balls[j])
                    if time < min_time:
                        min_time = time
                        a = self.__balls[i]
                        b = self.__balls[j]
        
        # moving all balls that did not collide
        for i in range(len(self.__balls)):
            if (self.__balls[i] != a) and (self.__balls[i] != b):
                self.__balls[i].move(min_time)
        
        # returning the time between the previous collision and this collision
        # returning the change in momentum of any balls that collided with the wall
        # these are returned to calculate the pressure
        if b == self.__container:
            p_before = a.p()
            a.collide(b)
            p_after = a.p()
            delta_p_mag = np.linalg.norm(p_after - p_before)
            return(min_time, delta_p_mag)
        # if it did not collide with the container than 0 is returned for the change in momentum
        # this is because only the collisions with the wall contribute to the pressure
        else:
            a.collide(b)
            return(min_time, 0)
        
    # run function runs the simulation for the specified number of frames
    # each frame is a collision
    # has various options on what to do/calculate during the simulation
    # e.g. animate the simulation, find pressure, temperature etc
    def run(self, num_frames, animate=False, pressure=False, distances_from_center=False, \
            distances_between_balls=False, kinetic_energy=False, momentum=False, remit=False, \
                temperature=False, maxwell=False):
        
        #  defining variables to be used when calculating specific things during the simulation
        d_center = []
        d_balls = []
        total_p = 0
        lim = -1*self.__container.rad()
        times = []
        time = 0
        system_kes = []
        system_ps = []
        speeds = []
        
        # if animating then this displays the image of the system each frame
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-lim, lim), ylim=(-lim, lim))
            ax.add_artist(self.__container.get_patch())
            for i  in range(len(self.__balls)):
                ax.add_patch(self.__balls[i].get_patch())
        
        # goes through each frame to run the simulation
        # if neccessary certain information is collected for each frame
        for frame in (range(num_frames)):
            # performs next collision
            t, p = self.next_collision()
            
            # defines veriables to collect data
            system_ke = 0
            system_p = 0
            time += t
            times.append(time)
            
            # calculate total pressure
            if pressure:
                total_p += p
            
            # for loop to collect infomation for every ball at every frame
            for i in range(len(self.__balls)):
                if distances_from_center:
                    distance = np.linalg.norm(self.__balls[i].pos())
                    d_center.append(distance)
                if distances_between_balls:
                    for j in range(len(self.__balls)):
                        if i != j:
                            distance = np.linalg.norm(self.__balls[i].pos()-self.__balls[j].pos())
                            d_balls.append(distance)
                if kinetic_energy:
                    system_ke += self.__balls[i].ke()
                if momentum:
                    system_p += self.__balls[i].p()
                if maxwell:
                    speeds.append(np.linalg.norm(self.__balls[i].vel()))
            
            # collecting data of the whole system every frame
            if kinetic_energy:
                system_kes.append(system_ke)
            if momentum:
                system_ps.append(np.linalg.norm(system_p))
            
            # time per frame if animating
            if animate:
                pl.pause(0.0001)
        
        # showing animation if animating
        if animate:
            pl.show()
        
        # printing pressure if requested
        # option to return pressure instead if requested
        if pressure:
            force = total_p/times[-1]
            pressure = force/(2*np.pi*lim)
            if remit == True:
                return pressure
            else:
                print(pressure, 'Pa')
        
        # plotting histogram of distances from the centre of each ball if requested
        if distances_from_center:
            plt.hist(d_center)
            plt.xlabel('Distances from the centre (m)')
            plt.ylabel('Frequency')
            # plt.title('Histogram of ball distance from container centre')
            print(len(self.__balls), 'balls and', num_frames, 'frames')
        
        # plotting histogram of distances between balls for each ball if requested
        if distances_between_balls:
            plt.hist(d_balls)
            plt.xlabel('Inter-ball separation (m)')
            plt.ylabel('Frequency')
            # plt.title('Histogram of inter-ball separation \n %i balls')
            print(len(self.__balls), 'balls and', num_frames, 'frames')
        
        # plotting kinetic energy each frame if requested
        if kinetic_energy:
            plt.plot(times, system_kes, scaley=False)
            plt.ylim([system_kes[0]-10,system_kes[0]+10])
            plt.xlabel('Time (s)')
            plt.ylabel('System kinetic energy (J)')
            # plt.title('Plot of system kinetic energy with time')
            print(len(self.__balls), 'balls and', num_frames, 'frames')
        
        # plotting momentum each frame if requested
        if momentum:
            plt.plot(times, system_ps, 'x')
            plt.xlabel('Time (s)')
            plt.ylabel('System momentum $(kgms^{-1})$')
            # plt.title('Plot of system momentum with time')
            print(len(self.__balls), 'balls and', num_frames, 'frames')
        
        # printing pressure if requested
        # option to return pressure instead if requested
        if temperature:
            ke = 0
            for i in range(len(self.__balls)):
                self.__balls[i].ke()
                ke += self.__balls[i].ke()
            ke_avg = ke/len(self.__balls)
            T = ke_avg/1.38064852e-23
            if remit == True:
                return T
            else:
                print(T)
         
        # plotting histogram of ball speeds for every ball in every frame
        # plots maxwelll-boltzmann distribution with correct temperature and scaling
        if maxwell:
            # plotting histogram
            bin_heights, bin_edges, patches = plt.hist(speeds)
            plt.xlabel('Speed of balls $(ms^{-1})$')
            plt.ylabel('Frequency')
            # plt.title('Plot of speed distribution')
            print(len(self.__balls), 'balls and', num_frames, 'frames')
            
            # defining maxwell-boltzmann distribution
            def distribution(v,area,T):
                    return area * v * np.exp((-0.5*1*v**2)/(1.38064852e-23*T))
            
            # finding the area of the histogram for scaling the distribution
            area = 0
            for i in range(len(bin_heights)):
                area += bin_heights[i] * (bin_edges[i+1]-bin_edges[i])
                
            # finding temperature so the theoretical distribution matches the histogram
            ke = 0
            for i in range(len(self.__balls)):
                self.__balls[i].ke()
                ke += self.__balls[i].ke()
            ke_avg = ke/len(self.__balls)
            T = ke_avg/1.38064852e-23
            
            # plotting theoretical distribution
            x = np.linspace(0,max(speeds),100)
            plt.plot(x, distribution(x, area, T))
            plt.legend(['Theoretical prediction','Simulated ball speeds'])