# -*- coding: utf-8 -*-
"""
Simple Orbit Integrator

This program simulates a one-body solar system; a single point mass orbiting around a larger point mass.
These represent a planet orbiting around its parent star. This is done by using an orbit
 integrator algorithm.
"""

import scipy as sci
import astropy.constants as const

def v(r):

    vel = sci.sqrt(G*m1/(orbit_rad*2.0))

    return vel

def f(r):

    r=sci.linalg.norm(r2-r1) #Calculate magnitude or norm of vector

    force = -1 * ((const.G.value * m1 * m2) / sci.mod(r)**3) * r

    return force

def vn_next_half(vn, xn, dt):

    next_half_vn = vn + f(xn) * dt/2.0

    return next_half_vn

def xn_next(xn, next_half_vn, dt):

    next_xn = xn + next_half_vn * dt

    return next_xn

def vn_next(next_half_vn, next_xn, dt):

    next_vn = next_half_vn + f(next_xn) * dt/2.0

    return next_vn


# class Star(star_mass):

#     mass = star_mass
#     pos = (0.0,0.0,0.0)

# class Planet(planet_mass, x_pos, y_pos, z_pos):

#     mass = planet_mass
#     pos = (x_pos, y_pos, z_pos)

iterations = 1000

sol_mass = const.M_sun.value
kpc = const.kpc.value
G = 43007.105731706317

orbit_rad = 8.0 # kpc
apo = 8.0 # kpc
peri = 1.0 # kpc
dt = 1E-3 # ~ E-3 Gyr

coords = sci.matrix(3)

print(planet_mass)
print(orbit_rad)

# Masses of 
m1 = 9.0 # central mass in E10 solar masses
m2 = 1.0 # object mass in E10 solar masses
# ^ random value

# Starting positions
r1 = sci.array([0.0,0.0,0.0], dtype='float64')
r2 = sci.array([8.0,0.0,0.0], dtype='float64')

# Initial velocities
v1 = sci.array([0.0,0.0,0.0])
v2 = sci.array([0.0,10.0,0.0])

v_com = (m1*v1+m2*v2)/(m1+m2)

for i in range(iterations):

