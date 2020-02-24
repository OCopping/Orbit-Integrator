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
from pprint import pprint
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation

# Function to define the equation of circulae velocity
def v_circ(r):

    vel = sci.sqrt(G*m1/orbit_rad)

    return vel

# Function to define the equation of elliptical velocity
def v_elip(r):

    vel = sci.sqrt(G*m1*(2.0/r * 1/a))

    return vel

def semi_maj(r, e):

    semi_maj_axis = (r*(1 + e))/(1 - e**2)

    return semi_maj_axis

# Function to define the equation of relative motion, which calculates the 
# acceleration of the orbiting body
def accel(xn):

    print('r1',r1)

    r = xn-r1
    r_mag = sci.linalg.norm(r) #Calculate magnitude or norm of vector

    accel = -1 * ((G * m1) / r_mag**3) * r

    return accel


def xn_next(xn, next_half_vn, dt):

    next_xn = xn + next_half_vn * dt

    return next_xn

def vn_next(vn, accel, dt):

    next_vn = vn + accel * dt/2.0

    return next_vn

def calc_time_step(central_pos,obj_pos,accel):

    eta = 0.02

    dt = eta * sci.sqrt(abs(sci.linalg.norm(obj_pos)- \
        sci.linalg.norm(central_pos))/sci.linalg.norm(accel))

    #print(dt)

    return dt



#iterations = 20000

sol_mass = const.M_sun.value
kpc = const.kpc.value
G = 43007.105731706317

orbit_rad = 8.0 # kpc
apo = 8.0 # kpc
peri = 1.0 # kpc
e = 0.1

a = semi_maj(peri,e)
print(a)

# Masses
m1 = 9.0 # central mass in E10 solar masses
#m2 = 1.0 # object mass in E10 solar masses
# ^ random value

t_max = 1500
t_current = 0.0
dt = (2.0*sci.pi*orbit_rad)/v_elip(peri) * 0.01 # Gyr

# Starting positions
r1 = sci.array([a*e,0.0,0.0], dtype='float64')
r2 = sci.array([peri+a*e,0.0,0.0], dtype='float64')

#points = sci.repeat(r2[:,0],10)
points = sci.array([r2,r2,r2,r2,r2,r2,r2,r2,r2,r2,r2,r2,r2,r2,r2]) #15 points

# Initial velocities
print(v_elip(peri))
print(v_elip(apo))
v1 = sci.array([0.0,0.0,0.0], dtype='float64')
v2 = sci.array([0.0,v_elip(peri),0.0], dtype='float64')

#v_com = (m1*v1+m2*v2)/(m1+m2)

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)
line, = ax.plot([], [], [])

# Setting the axes properties
ax.set_xlim3d([-10.0, 10.0])
ax.set_xlabel('X')

ax.set_ylim3d([-10.0, 10.0])
ax.set_ylabel('Y')

ax.set_zlim3d([-10.0, 10.0])
ax.set_zlabel('Z')

ax.set_title('Orbital integrator test')


def init():
    line.set_data([], [])
    line.set_3d_properties([])
    return line,


## ANIMATED ORBITAL TRAJECTORY 
def animate_orbit(t):

    global r1, r2, v2, t_current, points

    if t_current > t:
        exit()

    accel1 = accel(r2)
    dt1 = calc_time_step(r1,r2,accel1)

    #for index, dimension in enumerate(v2):
    half_vn = vn_next(v2, accel1, dt1)
    next_xn = xn_next(r2, half_vn, dt1)
    r2 = next_xn

    pprint(r2)
    points = sci.concatenate(([r2],points[0:-1]), axis=0)
    
    accel2 = accel(next_xn)
    dt2 = calc_time_step(r1,r2,accel2)

    v2 = vn_next(half_vn, accel2, dt2)
        
    t_current += dt2

    line.set_data(points[:,0], points[:,1])
    line.set_3d_properties(points[:,2])

    return line,

ax.scatter3D(r1[0],r1[1],r1[2],color='r',s=10)
#plt.plot(data[:,0], data[:,1], data[:,2])

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, animate_orbit, frames=int(t_max), 
                                   init_func=init, interval=5, blit=False)

plt.show()

#plt.savefig('3dtest.png')