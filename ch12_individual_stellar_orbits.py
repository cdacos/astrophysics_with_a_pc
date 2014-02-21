# -*- coding: utf-8 -*-

"""
Chapter 12 - Individual Stellar Orbits in the Galaxy

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-21)

Examples:
$ python ch12_individual_stellar_orbits.py 10 0 0 150 180
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow, log10, pi, sin, tan, cos

def effes(r, z, e, fe, a0, dc, kr, kz):
  """ evaluates the right hand sides of the differentia; equation of motion
      for the coordinates r and z. Both right hand sides are computed in the
      same cycle applying Simpson's rule for the numerical calculation of a
      defined integral
  """
  k1 = .0
  k2 = .0
  db = fe / 20.0
  b1 = .0
  b2 = db
  b3 = 2 * db
  fr1 = .0
  fz1 = .0

  # next for cycle computes ten pairs of subintervals for numerical 
  # integration with Simpson's rule
  for ii in range(1, 11):
    a2 = 1.0 / e * sqrt(pow(r * sin(b2), 2) + pow(z * tan(b2), 2))
    a3 = 1.0 / e * sqrt(pow(r * sin(b3), 2) + pow(z * tan(b3), 2))
    fr2 = exp(-a2 / a0) * pow(sin(b2), 2)
    fz2 = exp(-a2 / a0) * pow(tan(b2), 2)
    fr3 = exp(-a3 / a0) * pow(sin(b3), 2)
    fz3 = exp(-a3 / a0) * pow(tan(b3), 2)
    k1 = k1 + db * (fr1 + 4.0 * fr2 + fr3) / 3.0
    k2 = k2 + db * (fz1 + 4.0 * fz2 + fz3) / 3.0
    b1 = b1 + 2 * db
    b2 = b2 + 2 * db
    b3 = b3 + 2 * db
    fr1 = fr3
    fz1 = fz3

  # multiple the two integrals with the other factors of the equations
  # of motion
  kr = -4 * pi * sqrt(1 - e * e) / pow(e, 3) * dc * r * k1
  kz = -4 * pi * sqrt(1 - e * e) / pow(e, 3) * dc * z * k2

  return [r, z, e, fe, a0, dc, kr, kz]

print('Astrophysics with a PC : INDIVIDUAL STELLAR ORBITS')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of initial conditions and parameters : ')
r   = start_parameter('Initial conditions :  r(0) : ', 1)
z   = start_parameter('Initial conditions :  z(0) : ', 2)
u   = start_parameter('Initial conditions :  u(0) : ', 3)
v   = start_parameter('Initial conditions :  v(0) : ', 4)
vt0 = start_parameter('Initial conditions : vt(0) : ', 5)

# initialize galaxy model parameters
dt = .001
e = .99
dc = 11613.5
a0 = 2.8
fe = 1.4292567

# initialize main cycle
h = r * vt0
past = 0
n = 15
print('  i      t         x         z         u         v')
t = .0
ni = 1

# here start main cycle computing blocks of 15 iterations
ch = ''
while ch != 's' and ch != 'S':
  for i in range(1, n+1):

    # results at half the step (i+1/2)
    r1 = r + .5 * dt * u
    z1 = z + .5 * dt * v
    [r, z, e, fe, a0, dc, kr, kz] = effes(r, z, e, fe, a0, dc, .0, .0)
    u1 = u + .5 * dt * (kr + h * h / r / r / r)
    v1 = v + .5 * dt * kz

    # results at the full step (i+1)
    r = r + dt * u1
    z = z + dt * v1
    [r1, z1, e, fe, a0, dc, kr1, kz1] = effes(r1, z1, e, fe, a0, dc, .0, .0)
    u = u + dt * (kr1 + h * h / r1 / r1 / r1)
    v = v + dt * kz1
    
    t = t + dt

    # show newly computed iteration on screen
    print('{: 3d} {: 9.4f} {: 9.4f} {: 9.4f} {: 9.4f} {: 9.4f}'.format(ni, t, r, z, u, v))
    ni = ni + 1

  ch = raw_input('Enter S to stop or any key to continue')



