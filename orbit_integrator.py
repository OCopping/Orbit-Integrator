# -*- coding: utf-8 -*-
"""
Simple Orbit Integrator

This program simulates a one-body solar system; a single point mass orbiting around a larger point mass.
These represent a planet orbiting around its parent star. This is done by using an orbit
 integrator algorithm.
"""

import scipy as sci
import scipy.linalg
import astropy.constants as const

def v(r):

    vel = sci.sqrt(G*m1/(orbit_rad*2.0))

    return vel

def a(xn):

    r = xn-r1
    print('r=',r)
    r_norm = sci.linalg.norm(xn-r1) #Calculate magnitude or norm of vector
    print('r_norm=',r_norm)

    accel = -1 * ((const.G.value * m1) / r**3) * r_norm
    print(accel)

    return accel


def vn_next_half(vn, xn, dt):

    next_half_vn = vn + a(xn) * dt/2.0

    return next_half_vn

def xn_next(xn, next_half_vn, dt):

    next_xn = xn + next_half_vn * dt

    return next_xn

def vn_next(next_half_vn, next_xn, dt):

    next_vn = next_half_vn + a(next_xn) * dt/2.0

    return next_vn


# class Star(star_mass):

#     mass = star_mass
#     pos = (0.0,0.0,0.0)

# class Planet(planet_mass, x_pos, y_pos, z_pos):

#     mass = planet_mass
#     pos = (x_pos, y_pos, z_pos)

iterations = 10

sol_mass = const.M_sun.value
kpc = const.kpc.value
G = 43007.105731706317

orbit_rad = 8.0 # kpc
apo = 8.0 # kpc
peri = 1.0 # kpc
dt = 1 # Gyr

coords = sci.matrix(3)

# Masses
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

with open('orbit_data', 'w') as orbit_file:

    orbit_file.write('# r1 v1 r2 v2\n')

    data = scipy.concatenate((r1,v1,r2,v2))
    orbit_file.write(str(data)+'\n')

    for i in range(iterations):

        #for index, dimension in enumerate(v2):
        half_vn = vn_next_half(v2, r2, dt)
        next_xn = xn_next(r2, half_vn, dt)
        r2 = next_xn
        v2 = vn_next(half_vn, next_xn, dt)
            
        #sci.savetxt(orbit_file, (r1,v1,r2,v2))
        print(str(sci.concatenate(r1,v1,r2,v2)))
        
