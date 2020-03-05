# -*- coding: utf-8 -*-
"""
Simple Orbit Integrator graphing script

This program takes the data file 'variables.dat' that was produced by the
orbit_integrator.py script and produces graphs of energy, radius and velocity
against time.
"""

import scipy as sci
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation

data = sci.genfromtxt('variables.dat', skip_header=True)

fig = plt.figure(1)
plt.title('Gravitational Potential Energy over Time')
plt.xlabel('Time [Gyr]')
plt.ylabel('Energy []')
plt.plot(data[:,0], data[:,1])
plt.savefig('energy_vs_time.png')

fig = plt.figure(2)
plt.title('Orbital Radius over Time')
plt.xlabel('Time [Gyr]')
plt.ylabel('Radius [kpc]')
plt.plot(data[:,0], data[:,2])
plt.savefig('radius_vs_time.png')

fig = plt.figure(3)
plt.title('Orbital Velocity over Time')
plt.xlabel('Time [Gyr]')
plt.ylabel('Velocity [km/s]')
plt.plot(data[:,0], data[:,3])
plt.savefig('velocity_vs_time.png')

fig = plt.figure(4, dpi=500)
plt.title('X vs Y of Orbit')
plt.xlabel('x [kpc]')
plt.ylabel('y [kpc]')
plt.scatter(data[:,4],data[:,5],s=0.1,)
plt.savefig('x_vs_y.png')