# -*- coding: utf-8 -*-

"""
Chapter 3 - Meteor Dynamics

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-07)
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log

def datm(y):
  """Computes the atmospheric density at height y(cm)
  """
  return exp(-6.65125 - 1.39813e-6 * y)

def effes(x, y, u, v, m, k1, k2, tau):
  """Computes fx, fy, fu, fv, fm and s (the speed) for a position (x,y)
  velocity (u,v) and mass (m), and with parameters k1,k2 and tau
  """
  fx = u
  fy = v
  s = sqrt(u**2 + v**2)
  rho = datm(y)
  print('u={:f}, v={:f}, m={:f}, s={:f}, rho={:f}'.format(u, v, m, s, rho))
  fu = -k1 * rho * s * u * exp(-1/3 * log(m))
  fv = -k1 * rho * s * v * exp(-1/3 * log(m)) - 980
  fm = -k2 * rho * s**3 * exp(2/3 * log(m))
  return [fx, fy, fu, fv, fm, s]

print('Astrophysics with a PC : METEOR')
print('------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of initial conditions and parameters : ')
y   = start_parameter('Initial height (km)             : ', 1)
u   = start_parameter('Initial horizontal speed (km/s) : ', 2)
v   = start_parameter('Initial vertical speed (km/s)   : ', 3)
m   = start_parameter('Initial mass (gram)             : ', 4)
k1  = start_parameter('Parameter K1  : ', 5)
k2  = start_parameter('Parameter K2  : ', 6)
tau = start_parameter('Parameter tau : ', 7)

# transform input data x, y, u and v from km to cm
# and make sure that y is negative
y = y + 100000
x = .0
u = u * 100000
v = -1 * fabs(v * 100000)
minit = m
t = .0
i = 1.0
print('')
print('   i   t    x    y    u    v    m    mag')
print (' ')

while m >= minit: # This is the main loop of the program

  # next 3 nested if-blocks select the time the time step dt
  if m > .8 * minit:
    dt = .1
  elif m > .5 * minit:
    dt = .05
  elif m > .35 * minit:
    dt = .02
  else:
    dt = .01

  t = t + dt

  # next line computes right hand sides of the 5 differential equations (state i)
  fx, fy, fu, fv, fm, s = effes(x, y, u, v, m, k1, k2, tau)

  # next 5 instructions compute the predicted state 'i+1
  x1 = x + dt * fx
  y1 = y + dt * fy
  u1 = u + dt * fu
  v1 = v + dt * fv
  m1 = m + dt * fm

  # next line computes rigt hand sides of the 5 differential equations (state 'i+1)
  fx1, fy1, fu1, fv1, fm1, s1 = effes(x1, y1, u1, v1, m1, k1, k2, tau)

  # next 5 instructions compute the correct state at 'i+1
  x = x + .5 * dt * (fx + fx1)
  y = y + .5 * dt * (fy + fy1)
  u = u + .5 * dt * (fu + fu1)
  v = v + .5 * dt * (fv + fv1)
  m = m + .5 * dt * (fm + fm1)

  # Computation of the apparent magnitude :
  print('fm1={:f}, s1={:f}'.format(fm1, s1))
  e = -.5 * tau * fm1 * s1**2
  mag = 5 * log(y, 10) - 2.5 * log(e, 10) - 8.795

  # Results are shwon on screen :
  if y > 0:
    print('{: 2d}/', i, t, x / 100000, y / 100000, u / 100000, v / 100000, m, mag)
  else:
    print('')
    print('Meteoroid has reached the ground')
    m = 0

  i = i + 1

  if i % 15 == 0:
    # this if makkes the program pause after every 15 iterations
    print('')
    raw_input('press enter to continue')
    print('')

print('')
print('Program terminated.')

