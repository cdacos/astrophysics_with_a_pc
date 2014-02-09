# -*- coding: utf-8 -*-

"""
Chapter 3 - Meteor Dynamics

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-07)

Example:
$ python ch03_meteor_dynamics.py 160 20 40 0.01 1 1.1e-11 0.02
According to Table 3-1 should produce:
0.10  2.0000 156.0000 19.99999 -40.00095 0.0099986 10.77
1.00 19.9992 119.9968 19.99488 -39.99956 0.0095010  4.77
1.45 28.9885 102.0126 19.93035 -39.87487 0.0049890  2.17
1.60 31.9708  96.0459 19.80846 -39.63245 0.0014891  2.05
1.70 33.9402  92.1055 19.50003 -39.01633 0.0000793  3.69
This is not what this program outputs. No explanation yet.
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow

def datm(y):
  """Computes the atmospheric density at height y(cm)
  """
  return exp(-6.65125 - 1.39813e-6 * y)

def effes(u, v, m, k1, k2):
  """Computes fx, fy, fu, fv, fm and s (the speed) for a position (x,y)
  velocity (u,v) and mass (m), and with parameters k1,k2 and tau
  """
  fx = u
  fy = v
  s = sqrt(u**2 + v**2)
  rho = datm(y)
  fu = -k1 * rho * s * u * pow(m, -1.0/3)
  fv = -k1 * rho * s * v * pow(m, -1.0/3) - 980
  fm = -k2 * rho * s*s*s * pow(m, 2.0/3)
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
y = y * 100000
x = .0
u = u * 100000
v = -1 * fabs(v * 100000)
minit = m
t = .0
i = 1
print('')
print(' i    t      x        y        u        v        m        mag')
print (' ')

while m >= minit * .01: # This is the main loop of the program

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
  fx, fy, fu, fv, fm, s = effes(u, v, m, k1, k2)

  # next 5 instructions compute the predicted state 'i+1
  x1 = x + dt * fx
  y1 = y + dt * fy
  u1 = u + dt * fu
  v1 = v + dt * fv
  m1 = m + dt * fm

  # next line computes rigt hand sides of the 5 differential equations (state 'i+1)
  fx1, fy1, fu1, fv1, fm1, s1 = effes(u1, v1, m1, k1, k2)

  # next 5 instructions compute the correct state at 'i+1
  x = x + .5 * dt * (fx + fx1)
  y = y + .5 * dt * (fy + fy1)
  u = u + .5 * dt * (fu + fu1)
  v = v + .5 * dt * (fv + fv1)
  m = m + .5 * dt * (fm + fm1)

  # Computation of the apparent magnitude :
  e = -.5 * tau * fm1 * s1**2
  mag = 5. * log(y, 10) - 2.5 * log(e, 10) - 8.795

  # Results are shwon on screen :
  if y > 0:
    print('{:2d} {:5.2f} {:8.4f} {:8.4f} {:9.5f} {:9.5f} {:9.7f} {:5.2f}'.format(i, t, x / 100000, y / 100000, u / 100000, v / 100000, m, mag))
  else:
    print('')
    print('Meteoroid has reached the ground')
    m = 0

  if i % 10 == 0:
    # this if makkes the program pause after every 15 iterations
    print('')
    raw_input('press enter to continue')
    print('')

  i = i + 1

print('')
print('Program terminated.')

