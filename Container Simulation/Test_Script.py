#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 13:06:40 2021

@author: dawudabd-alghani
"""

import numpy as np
import matplotlib.pyplot as plt
from Ball import Ball
from Simulation import Simulation

#%%

# Checking that time_to_collision works by moving them by
# that amount of time and seeing if they are colliding
# Radii set to 0 to make it easy to see if they are colliding
# (They will be at the same position if they are)

ball1 = Ball(1, 0, [1,1], [1,1])
ball2 = Ball(1, 0, [-5,0], [7,2])
t = ball1.time_to_collision(ball2)
print(t)
ball1.move(t)
ball2.move(t)
print(ball1.pos())
print(ball2.pos())

# They have the same position so they are colliding
# time_to_collision method works

#%%

# Checking that the collide function moves the balls to
# the point where they collide

ball1 = Ball(1, 1, [1,1], [1,1])
ball2 = Ball(1, 2, [-5,0], [5,0])
ball1.collide(ball2)
print(ball1.pos(), ball1.vel())
print(ball2.pos(), ball2.vel())

# Results were checked manually to find the final 
# positions and velocities of both balls
# Results from paper match the function's results
# collide method works

#%%

# Checking conservation of ke with ball-ball collisions

ball1 = Ball(1, 1, [1,1], [1,1])
ball2 = Ball(4, 2, [-5,0], [7,2])

ke_before_ball1 = ball1.ke()
ke_before_ball2 = ball2.ke()

ball1.collide(ball2)

ke_after_ball1 = ball1.ke()
ke_after_ball2 = ball2.ke()

print(ke_before_ball1 + ke_before_ball2)
print(ke_after_ball1 + ke_after_ball2)

# Both are 107.0 so ke is conserved

#%%

# Checking conservation of ke with ball-container collisions
# Assuming ke of container is 0 throughout

ball1 = Ball(1, 1, [1,1], [1,1])
container = Ball(1e20, -10)

ke_before_ball1 = ball1.ke()

ball1.collide(container)

ke_after_ball1 = ball1.ke()

print(ke_after_ball1)

# The ke before is 1.0 and afterwards is 0.9999999999999996.
# Ke is conserved but the small discrepency is due to the
# container having a finite mass so it gains some ke
# The discrepancy is negligable.

#%%

# Testing the simulation works for many balls

container = Ball(1e20, -10)
ball1 = Ball(1, 1, [-6,0], [1,2])
ball2 = Ball(1, 1, [0,8], [1,-1], 'blue')
ball3 = Ball(1, 1, [-3,0], [1,2], 'green')
ball4 = Ball(1, 1, [-1,2], [5,2], 'orange')
ball5 = Ball(1, 1, [5,0], [1,-3], 'black')
ball6 = Ball(1, 1, [0,5], [1,2], 'yellow')
simulation = Simulation(container, [ball1,ball2,ball3,ball4,ball5,ball6])
simulation.run(100, animate=True)

# Simulation works (checked by eye)

#%%

# Testing simulation for automatically systematically placing many balls 

container = Ball(1e20, -100)
balls = Ball(n=500,rad=2).balls
simulation = Simulation(container,balls)
simulation.run(0, animate=True) # 0 frames to see the initial placement of balls
print(len(balls), 'balls')

# n balls placed without overlap so method works

#%%

# Checking the pressure from many simulations with different numbers of
# frames to check the simulation 

container = Ball(1e20, -100)
balls = Ball(n=10).balls
simulation = Simulation(container,balls)
print(len(balls), 'balls and 1000 frames:')
simulation.run(1000, pressure=True)
print(len(balls), 'balls and 10000 frames:')
simulation.run(10000, pressure=True)
print(len(balls), 'balls and 100000 frames:')
simulation.run(100000, pressure=True)

# Pressure should be independent of time (if the time frame is large enough)
# Program returns 0.008424385501520122, 0.008431384571136142, 0.007840230703778745 Pa
# These are very close to each other. Given higher time frames they would 
# become closer as random error is reduced.
# The numbers will be different each time the program is run due to the random
# velocities generated, but the pressures will always be close to each other

#%%

# Testing validity of the simulation by plotting a histogtram of the
# distances of balls from the center

container = Ball(1e20, -100)
balls = Ball(n=100).balls
simulation = Simulation(container,balls)
simulation.run(1000, distances_from_center=True)

# most balls at the edges as they are colliding with the wall
# more balls further out as larger area

#%%

# Testing validity of the simulation by plotting a histogtram of the
# distances between each pair of balls

container = Ball(1e20, -100)
balls = Ball(n=100).balls
simulation = Simulation(container,balls)
simulation.run(1000, distances_between_balls=True)

# 

#%%

# Checking kinetiic energy is conserved

container = Ball(1e20, -100)
balls = Ball(n=10).balls
simulation = Simulation(container,balls)
simulation.run(1000, kinetic_energy=True)

# System kinetic energy is constant so kinetic energy is conserved

#%%

# Checking momentum is conserved

container = Ball(1e20, -100)
balls = Ball(n=100).balls
simulation = Simulation(container,balls)
simulation.run(1000, momentum=True)

# The momentum is fluctuating as the container does not have infinite mass
# Some momentum is transferred between the container and balls
# There is overall increase or decrease in momentum so it is conserved

#%%

# Pressure against temperature graph

temperatures = []
pressures = []
container = Ball(1e20, -100)
for v in range(1,20):
    balls = Ball(n=20,v=v).balls
    simulation = Simulation(container,balls)
    pressures.append(simulation.run(1000, pressure=True, remit=True))
    temperatures.append(simulation.run(1000, temperature=True, remit=True))
plt.plot(temperatures, pressures, 'x')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
# plt.title('Plot of pressure against temperature')
coefs, cov = np.polyfit(temperatures,pressures,1, cov=True)
ffit = np.poly1d(coefs)
plt.plot(np.array(temperatures),ffit(np.array(temperatures)))
print('Gradient:', coefs[0], '±', np.sqrt(cov[0,0]))
print('Value of k from gradient =', coefs[0]*((np.pi*100**2)/len(balls)), '±', np.sqrt(cov[0,0])*((np.pi*100**2)/len(balls)))
print(len(balls), 'balls and 1000 frames')

# The gradient was used to find the value of k, and it is accurate to 3sf

#%%

# Pressure against temperature graph for different areas

for r in range(100,201,50):
    temperatures = []
    pressures = []
    container = Ball(1e20, -r)
    for v in range(1,20):
        balls = Ball(n=10,rad=j,v=v).balls
        simulation = Simulation(container,balls)
        pressures.append(simulation.run(1000, pressure=True, remit=True))
        temperatures.append(simulation.run(1000, temperature=True, remit=True))
    coefs, cov = np.polyfit(temperatures,pressures,1, cov=True)
    ffit = np.poly1d(coefs)
    plt.plot(np.array([0,4e25]),ffit(np.array([0,4e25])), label=j)
    print('gradient for container with radius', r, ':', coefs[0], '±', np.sqrt(cov[0,0]))
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
# plt.title('Plot of pressure against temperature for different container areas')
plt.legend(['area = 31416 $m^{2}$','area = 70686 $m^{2}$','area = 125664 $m^{2}$'])
print(len(balls), 'balls and 1000 frames')

# Gradient becomes less steep as the area increases
# This is expected as gradient*area must be constant so they are inversely proportional

#%%

# Pressure against temperature graph for different numbers of balls

for n in range(10,31,10):
    temperatures = []
    pressures = []
    container = Ball(1e20, -100)
    for v in range(1,20):
        balls = Ball(n=n,rad=j,v=v).balls
        simulation = Simulation(container,balls)
        pressures.append(simulation.run(1000, pressure=True, remit=True))
        temperatures.append(simulation.run(1000, temperature=True, remit=True))
    coefs, cov = np.polyfit(temperatures,pressures,1, cov=True)
    ffit = np.poly1d(coefs)
    plt.plot(np.array([0,4e25]),ffit(np.array([0,4e25])), label=j)
    print('gradient for container with', n, 'balls :', coefs[0], '±', np.sqrt(cov[0,0]))
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
# plt.title('Plot of pressure against temperature for varying numbers of balls')
plt.legend(['number = 10','number = 20','number = 30'])
print('1000 frames')

# Gradient becomes steeper as the number of balls increases
# This is expected as gradient/balls must be constant so they are directly proportional

#%%

# Pressure against temperature graph for different radii

for j in [1,2,4,10]:
    temperatures = []
    pressures = []
    container = Ball(1e20, -100)
    for v in range(1,20):
        balls = Ball(n=10,rad=j,v=v).balls
        simulation = Simulation(container,balls)
        pressures.append(simulation.run(1000, pressure=True, remit=True))
        temperatures.append(simulation.run(1000, temperature=True, remit=True))
    coefs, cov = np.polyfit(temperatures,pressures,1, cov=True)
    ffit = np.poly1d(coefs)
    plt.plot(np.array([0,4e25]),ffit(np.array([0,4e25])))
    print('gradient for balls with radius', j, ':', coefs[0], '±', np.sqrt(cov[0,0]))
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
# plt.title('Plot of pressure against temperature for different ball radii')
plt.legend(['radius = 1','radius = 2','radius = 5','radius = 10'])
print(len(balls), 'balls and 1000 frames')

# as ball radius increases the gradient increases slightly.
# this is expected as a larger ball radius means a higher affect due to the balls
# having a larger volume, which reduces the effective volume of the container.

#%%

# Plotting the speed distribution of the balls

container = Ball(1e20, -100)
balls = Ball(n=250).balls
simulation = Simulation(container,balls)
simulation.run(2000, maxwell=True)

# The speed distribution follows a maxwell-boltzmann curve for a large enough system, as expected

#%%

# Finding a nd b constant in Van der Waals' law
# a relates to intermolecular forces and there are none so a = 0
# Can plot PT and find gradient to find b

temperatures = []
pressures = []
container = Ball(1e20, -100)
for v in range(1,20):
    balls = Ball(n=20,rad=2,v=v).balls
    simulation = Simulation(container,balls)
    pressures.append(simulation.run(1000, pressure=True, remit=True))
    temperatures.append(simulation.run(1000, temperature=True, remit=True))
plt.plot(temperatures, pressures, 'x')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
# plt.title('Plot of pressure against temperature')
coefs, cov = np.polyfit(temperatures,pressures,1, cov=True)
ffit = np.poly1d(coefs)
plt.plot(np.array(temperatures),ffit(np.array(temperatures)))
print('Gradient:', coefs[0], '±', np.sqrt(cov[0,0]))
print('a = 0')
print('b =', ((coefs[0]*np.pi*100**2)-(len(balls)*1.38064852e-23))/(coefs[0]*len(balls)), '±',\
      ((len(balls)*1.38064852e-23)/(coefs[0]**2*len(balls)))*cov[0,0], 'm^2')
print('This is for balls of radius 2. The value of b is dependant on ball radius.')
print(len(balls), 'balls and 1000 frames')

# The values of a and b were determined
# a = 0 as there are no simulated inter-ball forces
# b is dependant on the radius of balls
# a larger radius means a higher value of b, as the balls occupy a larger volume