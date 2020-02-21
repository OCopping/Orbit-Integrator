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
import pprint
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

def v_circ(r):

    vel = sci.sqrt(G*m1/orbit_rad)

    return vel

def v_elip(r):

    vel = sci.sqrt(G*m1*(2.0/r * 1/a))

    return vel

def semi_maj(r, e):

    semi_maj_axis = (r*(1 + e))/(1 - e**2)

    return semi_maj_axis

def accel(xn):

    r = xn-r1
    r_mag = sci.linalg.norm(r) #Calculate magnitude or norm of vector
    #print('r_mag =', r_mag)

    accel = -1 * ((G * m1) / r_mag**3) * r
    #print('a=',accel)

    return accel


def xn_next(xn, next_half_vn, dt):

    next_xn = xn + next_half_vn * dt

    return next_xn

def vn_next(vn, xn, dt):

    next_vn = vn + accel(xn) * dt/2.0

    return next_vn


iterations = 1000

sol_mass = const.M_sun.value
kpc = const.kpc.value
G = 43007.105731706317

orbit_rad = 8.0 # kpc
apo = 8.0 # kpc
peri = 1.0 # kpc
e = 0.12

a = semi_maj(peri,e)
print(a)

coords = sci.matrix(3)

# Masses
m1 = 9.0 # central mass in E10 solar masses
#m2 = 1.0 # object mass in E10 solar masses
# ^ random value

t_max = 1
t_current = 0.0
dt = (2.0*sci.pi*orbit_rad)/v_elip(peri) * 0.01 # Gyr


# Starting positions
r1 = sci.matrix([a*e,0.0,0.0], dtype='float64')
r1.shape = (3,1)
r2 = sci.matrix([peri+a*e,0.0,0.0], dtype='float64')
r2.shape = (3,1)

#pprint.pprint(r1)

# Initial velocities
# v1 = sci.array([0.0,0.0,0.0])
# v2 = sci.array([0.0,10.0,0.0])

#v1 = sci.asmatrix(v1, dtype='float64')
#v2 = sci.asmatrix(v2, dtype='float64')

v1 = sci.matrix([0.0,0.0,-0.2], dtype='float64')
v1.shape = (3,1)
v2 = sci.matrix([0.0,v_elip(peri),0.0], dtype='float64')
v2.shape = (3,1)

#print(v(8.0))

#v_com = (m1*v1+m2*v2)/(m1+m2)

with open('orbit_data.dat', 'w') as orbit_file:

    orbit_file.write('# t    x    y    z\n')
    orbit_file.write('{0} {1} {2}\n'.format(r2[0,0], r2[1,0], r2[2,0]))

    for i in range(iterations):

        #for index, dimension in enumerate(v2):
        half_vn = vn_next(v2, r2, dt)
        next_xn = xn_next(r2, half_vn, dt)
        r2 = next_xn
        v2 = vn_next(half_vn, next_xn, dt)
            
        t_current += dt

        orbit_file.write('{0} {1} {2}\n'.format(r2[0,0], r2[1,0], r2[2,0]))

data = sci.genfromtxt('orbit_data.dat')

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# Fifty lines of random 3-D lines
#data = [Gen_RandLine(25, 3) for index in range(50)]

# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()
#lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

# Setting the axes properties
ax.set_xlim3d([-10.0, 10.0])
ax.set_xlabel('X')

ax.set_ylim3d([-10.0, 10.0])
ax.set_ylabel('Y')

ax.set_zlim3d([-10.0, 10.0])
ax.set_zlabel('Z')

ax.set_title('3D Test')

# Creating the Animation object
#line_ani = animation.FuncAnimation(fig, update_lines, 25, fargs=(data, lines),
#                                   interval=50, blit=False)

#p3.scatter([0],[0],[0],color='r',s=100)
plt.plot(data[:,0], data[:,1], data[:,2])

plt.show()

plt.savefig('3dtest.png')