# -*- coding: utf-8 -*-
"""
Simple Orbit Integrator

This program simulates a one-body solar system; a single point mass orbiting around a larger point mass.
These represent a planet orbiting around its parent star. This is done by using an orbit
 integrator algorithm.
"""

import scipy as sci
import scipy.linalg
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
import argparse


#------------------------------------------------------------------------------


# Function to define the equation of circulae velocity
def v_circ(r):

    vel = sci.sqrt(G*m1/orbit_rad)

    return vel

def semi_maj(r, e):

    semi_maj_axis = (r*(1 + e))/(1 - e**2)

    return semi_maj_axis

# Function to define the equation of relative motion, which calculates the 
# acceleration of the orbiting body
def accel(xn):

    r = xn-r1
    r_mag = sci.linalg.norm(r) #Calculate magnitude or norm of vector

    accel = -1 * ((G * m1) / r_mag**3) * r

    return accel

def NFW_pot_accel(xn):

    M_nfw = 80.0 # *E10 sol masses
    rs = 16.0 # kpc
    c = 15.3

    r = xn-r1
    r_mag = sci.linalg.norm(r) #Calculate magnitude or norm of vector


    #accel = -r* (( c+1)*G*M_nfw / ((r_mag**3)*((sci.log(1+c)*(1+c))-c))) * (((r_mag+rs)*sci.log((r_mag+rs)/rs)-r) / (r_mag+rs))

    accel = -r*G*M_nfw*(1+c) / (r_mag**3*((sci.log(1+c)*(1+c))-c) )  * ( sci.log(1+r_mag/rs) - r_mag/(r_mag+rs))

    return accel


def xn_next(xn, next_half_vn, dt):

    next_xn = xn + (next_half_vn * dt)

    return next_xn

def vn_next(vn, accel, dt):

    next_vn = vn + (accel * dt/2.0)

    return next_vn

def calc_time_step(central_pos,obj_pos,accel):

    eta = 0.08

    dt = eta * sci.sqrt(abs(sci.linalg.norm(obj_pos)- \
        sci.linalg.norm(central_pos))/sci.linalg.norm(accel))

    return dt

def calc_energy(v, m, r):

    E = 0.5*sci.power(v, 2) - (G*m)/r

    return E

def calc_nfw_energy(v, m, xn):

    M_nfw = 80.0 # *E10 sol masses
    rs = 16.0 # kpc
    c = 15.3

    r = xn-r1
    r_mag = sci.linalg.norm(r) #Calculate magnitude or norm of vector


    E = -G*M_nfw/r_mag *(sci.log(1+r_mag/rs)/(sci.log(1+c)-c/(1+c)))

    return E


#------------------------------------------------------------------------------


# Take an argument for the type of potential for the simuation
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--pot_type', type=str, required=True, \
    help="Gravitational potential model to use [Newton, NFW]")
parser.add_argument('-t', '--time_step', type=float, required=False, \
    help="Time step length, can be fixed (e.g. 0.0001) or dynamic (0.0)", \
    default=0.0001)
args = parser.parse_args()

# G vakue for simplified units
G = 43007.105731706317

orbit_rad = 8.0 # kpc
apo = 8.0 # kpc
peri = 1.0 # kpc
e = 0.1

# Masses
m1 = 9.0 # central mass in E10 solar masses
#m2 = 1.0 # object mass in E10 solar masses
# ^ random value

t_max = 10 # ~Gyrs
t_current = 0.0
# const_dt = (2.0*sci.pi*orbit_rad)/v_circ(peri) * 0.01 # Gyr

# Starting positions
r1 = sci.array([0.0,0.0,0.0], dtype='float64')
r2 = sci.array([8.0,0.0,0.0], dtype='float64')

points = sci.array([r2]*5)# 5 points on the line

# Initial velocities
# v1 = sci.array([0.0,0.0,0.0], dtype='float64')
# v2 = sci.array([0.0,v_elip(peri),0.0], dtype='float64')
v2 =  sci.array([0.0,71.0,0.0], dtype='float64')

#v_com = (m1*v1+m2*v2)/(m1+m2)

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)
line, = ax.plot([], [], [])
time_text = ax.text(0.02, 0.95,0.01, s='', transform=ax.transAxes)

# Setting the axes properties
ax.set_xlim3d([-10.0, 10.0])
ax.set_xlabel('X')

ax.set_ylim3d([-10.0, 10.0])
ax.set_ylabel('Y')

ax.set_zlim3d([-10.0, 10.0])
ax.set_zlabel('Z')

ax.set_title('Orbital integrator test')


#------------------------------------------------------------------------------


def init():
    line.set_data([], [])
    line.set_3d_properties([])
    time_text.set_text('')
    return line,

## ANIMATED ORBITAL TRAJECTORY 
def animate_orbit(t):

    global r1, r2, v2, t_current, points

    if t_current > t_max:
        exit()


    # KDK algorithm
    if args.pot_type == 'Newton':
        accel1 = accel(r2)
    elif args.pot_type == 'NFW':
        accel1 = NFW_pot_accel(r2)

    if args.time_step == 0.0:
        dt = calc_time_step(r1,r2,accel1)
    elif args.time_step > 0.0:
        dt = float(args.time_step)

    half_vn = vn_next(v2, accel1, dt)
    next_xn = xn_next(r2, half_vn, dt)
    r2 = next_xn

    points = sci.concatenate(([r2],points[0:-1]), axis=0)
    
    if args.pot_type == 'Newton':
        accel2 = accel(next_xn)
    elif args.pot_type == 'NFW':
        accel2 = NFW_pot_accel(next_xn)

    v2 = vn_next(half_vn, accel2, dt)

    v = sci.linalg.norm(v2)
    r = sci.linalg.norm(r2)
    energy = calc_nfw_energy(v, m1, r)
        
    t_current += dt

    # Write values to file
    fileout.write('{0} {1} {2} {3} {4} {5} {6}\n'.format(t_current,energy,r,v,r2[0],r2[1],r2[2]))

    line.set_data(points[:,0], points[:,1])
    line.set_3d_properties(points[:,2])

    time_text.set_text('t = {0:8.5f} Gyr'.format(t_current))

    return line,


#------------------------------------------------------------------------------


ax.scatter3D(r1[0],r1[1],r1[2],color='r',s=10)

with open('variables.dat','w') as fileout:

    fileout.write('# Time (Gyr) | Energy () | Radius (kpc) | Velocity (km/s) | x | y | z\n')

    # Creating the Animation object
    line_ani = animation.FuncAnimation(fig, animate_orbit, frames=int(t_max), 
                                    init_func=init, interval=5, blit=False)

    plt.show()

# Set up formatting for the movie files
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

#line_ani.save('orbit_test.mp4', writer=writer)
#plt.savefig('3dtest.png')