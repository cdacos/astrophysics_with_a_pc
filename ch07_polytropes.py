# -*- coding: utf-8 -*-

"""
Chapter 7 - Polytropes

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-20)

Examples:
$ python ch07_polytropes.py 1.5 .05 2 3
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow, log10

print('Astrophysics with a PC : POLYTROPES')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of initial conditions and parameters : ')
n    = start_parameter('polytrope index : ', 1)
dr   = start_parameter('stepsize        : ', 2)
mass = start_parameter('mass            : ', 3)
rad  = start_parameter('radius          : ', 4)

# show heading of main table of results on screen
print('  i   x = r/rn       f           h        log(P/Pc)   log(d/dC)     l*mr')

# compute polytrope results in center of the star
x = 0
f = 1
h = 0
i = 0
p = pow(f, n + 1)
d = pow(f, n)
m = -x * x * h
print('{: 3d} {: 9.4f} {: 11.5f} {: 11.5f} {: 11.4f} {: 11.4f} {: 11.4f}'.format(i, x, f, h, log10(p), log10(d), m))

# compute first step
x = dr
f = 1.0 - pow(x, 2) / 6.0 + pow(x, 4) * n / 120.0
h = -x / 3.0 + pow(x, 3) * n / 30.0

i = 1
p = pow(f, n + 1)
d = pow(f, n)
m = -x * x * h
print('{: 3d} {: 9.4f} {: 11.5f} {: 11.5f} {: 11.4f} {: 11.4f} {: 11.4f}'.format(i, x, f, h, log10(p), log10(d), m))

# initialize main cycle
i = 2
verder = 1

while verder > 0:
  x12 = x + .5 * dr
  f12 = f + .5 * dr * h

  # check whether f12 (= F at half step) is till positive
  if f12 > 0:
    h12 = h + .5 * dr * (-pow(f, n) - 2 * h / x)
    x1 = x + dr
    f1 = f + dr * h12
    # check whether f1 ( = F at new state) is still positive
    if f1 > 0:
      h1 = h + dr * (-pow(f12, n) - 2 * h12 / x12)

      # compute pressure, density and mass
      p = pow(f1, n + 1)
      d = pow(f1, n)
      m = -x1 * x1 * h1

      x = x1
      f = f1
      h = h1

      # show results on screen
      print('{: 3d} {: 9.4f} {: 11.5f} {: 11.5f} {: 11.4f} {: 11.4f} {: 11.4f}'.format(i, x, f, h, log10(p), log10(d), m))

      # pause every 10 iterations
      if i > 0 and i % 10 == 0:
        raw_input('Press Enter to continue')

      # prepare for next iteration
      i = i + 1
    else:
      # this else is reached if f1 was negative
      verder = 0
  else:
    verder = 0

print('')
raw_input('Press Enter to proceed to the general characteristics')

# compute general characteristics and surface data
xm = x - f / h
hm = h + (xm - x) * (-pow(f, n) - 2 * h / x)
fm = 0
pc = 9.048e14 * mass * mass / (n + 1) / hm / hm / pow(rad, 4)
dm = 1.42 * mass / pow(rad, 3)
dc = -dm * xm / 3.0 / hm
lam = -xm * xm * hm / mass
rn = rad / xm

# show general charcateristics and surface data on screen
print('')
print('central pressure (Pc) : ', pc)
print('average density (dm)  : ', dm)
print('central density (dc)  : ', dc)
print('mass parameter (L)    : ', lam)
print('distance unit (rn)    : ', rn)
print('x-final               : {: 10.4f} '.format(xm))
print('-x2*h (final)         : {: 10.4f} '.format(-xm * xm * hm))

