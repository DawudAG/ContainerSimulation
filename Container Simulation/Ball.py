#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 19:40:21 2021

@author: dawudabd-alghani
"""

import numpy as np
import pylab as pl

class Ball():
    
    '''
    Class for creating a ball with a given mass, radius, position, velocity and colour (2D).
    Kinetic energy and patch attributes are created for the ball using the given attributes.
    move function moves a ball by a time dt, changing it's position.
    time_to_collision function finds the time for a ball to collide with another ball (or the container).
    collide function collides a ball with another ball (or the container) by moving both balls to the 
    place of the collision and changing their velocites to what they should be after the collision.
    '''
    
    def __init__(self, m=1, rad=0.01, pos = [0,0], vel = [0,0], colour = 'red', n=0, v=1):
        
        # if n is given then the function systematically places n balls by going through 
        # coordinates in a square shape and checking whether they lie in the container
        # the gap between coordinates is varied to achive the correct number of balls
        if n != 0:
            self.balls = []
            g = int(np.sqrt(31397/n))
            for i in range(-100,100,g):
                for j in range(-100,100,g):
                    if i**2 + j**2 < (100-rad)**2:
                        vel = np.random.normal(0,v,[1,2])[0]
                        ball = Ball(m, rad, [i,j], vel)
                        self.balls.append(ball)
        
        # defining a ball if the different attributes are given as parameters
        else:
            self.__m = m
            self.__rad = rad
            self.__pos = np.array(pos).astype(np.float)
            self.__vel = np.array(vel).astype(np.float)
            if rad < 0:
                self.__patch = pl.Circle(self.__pos, rad, fill=False)
            else:
                self.__patch = pl.Circle(self.__pos, rad, fc=colour)
    
    #  defining repr function for a ball
    def __repr__(self):
        return 'Ball(m=%s,rad=%s,pos=%s,vel=%s)' % (self.__m,self.__rad,self.__pos,self.__vel)
    
    # return radius of ball
    def rad(self):
        return self.__rad
    
    # return position of ball
    def pos(self):
        return self.__pos
    
    # return velocity of ball
    def vel(self):
        return self.__vel
    
    # calculate and return kinetic energy of ball
    def ke(self):
        return 0.5 * self.__m * np.dot(self.__vel, self.__vel)
    
    # calculate and return momentum of ball
    def p(self):
        return self.__m * self.__vel
    
    # return patch information of ball
    def get_patch(self):
        return self.__patch
    
    # move the ball by the distance it travels in a time, dt
    def move(self, dt):
        new_pos = self.__pos + self.__vel * dt
        self.__pos = new_pos
        self.__patch.center = new_pos
    
    # finds the time to collide between two balls
    # based on equation 1 in thermodynamics snookered script
    def time_to_collision(self, other):
        pos = self.__pos - other.__pos
        vel = self.__vel - other.__vel
        rad = self.__rad + other.__rad
        pos_mag = np.dot(pos,pos)
        vel_mag = np.dot(vel,vel)
        if vel_mag == 0 or np.dot(pos,vel)**2 - vel_mag\
            * (pos_mag - rad**2) < 0:
            return np.inf #They are parallel. Do not collide
        a = 1/vel_mag
        b = -1 * np.dot(pos,vel)
        c = np.sqrt(np.dot(pos,vel)**2 - vel_mag * (pos_mag - rad**2))
        t1 = a * (b - c)
        t2 = a * (b + c)
        # returns only the smallest positive value as that is the first collision in the future
        if t1 > 1e-5 and t2 > 1e-5:
            return(min(t1, t2))
        if t1 > 1e-5:
            return t1
        if t2 > 1e-5:
            return t2
        return np.inf #Collided in the past. Do not collide
    
    # uses the time_to_collision function and the move function
    # to move the balls to the point of collision
    # collide the balls by changing the velocitites of the balls
    # to what they should be after the collision
    def collide(self, other):
        t = self.time_to_collision(other)
        self.move(t)
        other.move(t)
        v1 = self.__vel
        v2 = other.__vel
        x1 = self.__pos
        x2 = other.__pos
        m1 = self.__m
        m2 = other.__m
        self.__vel = v1 - 2*m2/(m1+m2) * (np.dot(v1-v2,x1-x2)/\
                    np.dot(x1-x2,x1-x2)) * (x1-x2)
        other.__vel = v2 - 2*m1/(m1+m2) * (np.dot(v2-v1,x2-x1)/\
                    np.dot(x2-x1,x2-x1)) * (x2-x1)