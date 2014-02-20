# -*- coding: utf-8 -*-

"""
Chapter 4 - The Restricted Three-Body Problem

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-07)

Examples:
$ python ch04_restricted_3_bodies.py 0.000953875 -0.509046125 0.883345912 0.0258975212 0.0149272418 0.4
for trojan asteroids in Sun-Jupiter system (orbit A in text)

$ python ch04_restricted_3_bodies.py 0.000953875 -0.524046125 0.909326674 0.0646761399 0.0367068277 0.4
for trojan asteroids in Sun-Jupiter system (orbit B in text)

$ python ch04_restricted_3_bodies.py 0.0121396054 -0.4978603946 0.8833459119 0.0265752203 0.0146709149 0.4
for Earth-Moon system (orbit C in text)

$ python ch04_restricted_3_bodies.py 0.0121396054 -0.5128603946 0.9093266740 0.0682722747 0.0334034039 0.4
for Earth-Moon system (orbit D in text)

$ python ch04_restricted_3_bodies.py 0.000953875 -0.647717531 0.0 0.0 -0.6828143998 0.02
for Hilda in Sun-Jupiter system (ideal)

$ python ch04_restricted_3_bodies.py 0.000953875 -0.4952265404 -0.4163448036 0.4389046359 -0.5230661767 0.05
for Hilda in Sun-Jupiter system (more realistic)

$ python ch04_restricted_3_bodies.py 0.0000525 -0.6073955952 -0.7774968265 0.1083342234 -0.08463997159 0.05
for Pluto-Neptune
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow

def effes(x, y, u, v, mu):
  """Computes the right hand sides (fu and fv) of the differential equations
  for the components of the velocity for input values of x,y,u and v, and mu
  """
  r1 = sqrt((x - mu) * (x - mu) + y * y)
  r2 = sqrt((x + 1 - mu) * (x + 1 - mu) + y * y)
  fu = -(1 - mu) * (x - mu) / r1**3 - mu * (x + 1 - mu) / r2**3 + x + 2 * v
  fv = -(1 - mu) * y / r1**3 - mu * y / r2**3 + y - 2 * u
  return [fu, fv]

print('Astrophysics with a PC : RESTRICTED THREEBODY PROBLEM')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of initial conditions and parameters : ')
print('')
mu = start_parameter('Mass parameter mu         : ', 1)
x  = start_parameter('Initial conditions : x(0) : ', 2)
y  = start_parameter('                     y(0) : ', 3)
u  = start_parameter('                     u(0) : ', 4)
v  = start_parameter('                     v(0) : ', 5)
dt = start_parameter('Time step                 : ', 6)

t = dt
n = 20
ni = 1
ch = ''

while ch == '':
  print('   i       t          x            y             u            v')
  print('')

  # next for-cycle computes 20 new iterations
  for i in range(1, n+1):

    # next block computes half step values of Cauchy method
    x1 = x + .5 * dt * u
    y1 = y + .5 * dt * v
    fu, fv = effes(x, y, u, v, mu)
    u1 = u + .5 * dt * fu
    v1 = v + .5 * dt * fv

    # next blovk computes new state i+1 of Cauchy method
    x = x + dt * u1
    y = y + dt * v1
    fu, fv = effes(x1, y1, u1, v1, mu)
    u = u + dt * fu
    v = v + dt * fv

    # new state is shown on screen
    print('{:4d} {: 10.4f} {: 11.9f} {: 11.9f} {: 12.10f} {: 12.10f}'.format(ni, t, x, y, u, v))
    t = t + dt
    ni = ni + 1

  ch = raw_input('Press enter to continue, any other value to stop ')



