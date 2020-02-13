# -*- coding: utf-8 -*-
"""
Simple Orbit Integrator

This program simulates a one-body solar system; a single point mass orbiting around a larger point mass.
These represent a planet orbiting around its parent star. This is done by using an orbit
 integrator algorithm.
"""

import numpy as np
import astropy.constants as const

def f(x):

    val = 1

    return val

def vn_next_half(vn, xn, dt):

    next_half_vn = vn + f(xn) * dt/2.0

    return next_half_vn

def xn_next(xn, next_half_vn, dt):

    next_xn = xn + next_half_vn * dt/2.0

    return next_xn

def vn_next(next_half_vn, next_xn, dt):

    next_vn = next_half_vn + f(next_xn) * dt/2.0

    return next_vn


sol_mass = const.M_sun.value
mass = 9E10

print(sol_mass)