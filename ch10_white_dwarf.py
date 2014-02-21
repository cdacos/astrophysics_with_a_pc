# -*- coding: utf-8 -*-

"""
Chapter 10 - The Structure of White Dwarfs

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-21)

Examples:
$ python ch10_white_dwarf.py 8 80
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow, log10#, pi

def density(p, xold, rho, xnew):
  """ computes the density rho and the correspondingvalsue xnew from the
      input pressure P and starting guess xold
  """
  k = p / a
  stopcrit = 0
  while stopcrit != 1:
    xnew = xold - (f(xold) - k) / dfdx(xold)
    if abs(xold - xnew) < .000001:
      stopcrit = 1
    else:
      xold = xnew
  rho = 1964000.0 * pow(xnew, 3)
  return [p, xold, rho, xnew]

def dfdx(x):
  """ evaluates the derviatives of f(x)
  """
  return 8 * pow(x, 4) / sqrt(1 + pow(x, 2))

def f(x):
  """ evaluates the function f(x)
  """
  return x * (2 * pow(x, 2) - 3) * sqrt(1 + pow(x, 2)) + 3 * log(x + sqrt(1 + pow(x, 2)))

def mass(d, r):
  """ computes the right hand side of the differential equation of mass continuity
  """
  return 4 * pi * d  * r * r

def press(m, r, d):
  """ computes the right hand side of the equation of hydrostatic equilibrium
  """
  return -g * m * d / pow(r, 2)

# some physical constants used in the program :
g = 6.673e-8
m0 = 2e33
r0 = 6.96e10
a = 6.01e22
pi = 3.1415926536

print('Astrophysics with a PC : WHITE DWARF')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of initial conditions and model parameters : ')
rhc  = start_parameter('Log10 of the central density : ', 1)
drkm = start_parameter('Stepsize in km               : ', 2)

# central values of the physical variables
rhoc = pow(10, rhc)
dr = 100000.0 * drkm
print('')
m = .0
r = .0
xc = pow(rhoc / 1964000.0, 1 / 3.0)
pc = a * f(xc)
i = 0

# display heading of the table with the results
print(' i      r          Mr        log(P)     log(rho)       x')

# display central values (layer zero)
print('{: 3d} {: 10.7f} {: 10.7f} {: 10.7f} {: 10.7f} {: 10.7f}'.format(i, r / r0, m / m0, log10(pc), log10(rhoc), xc))

# compute first step from layer zero to layer one, and show them on screen
p = pc - 2 / 3.0 * g * pi * pow(rhoc * dr, 2)
m = 4 / 3.0 * pi * rhoc * pow(dr, 3)
[p, xc, d, x] = density(p, xc, .0, .0)
r = dr
i = 1
print('{: 3d} {: 10.7f} {: 10.7f} {: 10.7f} {: 10.7f} {: 10.7f}'.format(i, r / r0, m / m0, log10(p), log10(d), x))

# prepare for the other layers (from layer two to the surface)
surface = 0
i = 2

# here starts the main cycle for the other layers
while surface != 1:

  r1 = r + .5 * dr
  p1 = p + .5 * dr * press(m, r, d)

  # check if pressure (hald step i+1/2) is still positive
  # otherwise surface is reached
  if p1 > 0:
    m1 = m + .5 * dr * mass(d, r)
    [p1, x, d1, x1] = density(p1, x, .0, .0)
    r = r + dr
    p = p + dr * press(m1, r1, d1)

    # check if pressure (full step i+1) is still positive,
    # otherwise surface is reached
    if p > 0:
      m = m + dr * mass(d1, r1)
      [p, x1, d, x] = density(p, x1, d, x)

      # show results of newly computed layer on screen
      print('{: 3d} {: 10.7f} {: 10.7f} {: 10.7f} {: 10.7f} {: 10.7f}'.format(i, r / r0, m / m0, log10(p), log10(d), x))
      i = i + 1

      # pause every 10 layers
      if i % 10 == 1:
        raw_input('Press Enter to continue')
    else:
      surface = 1
  else: # this else is reached if pressure (i+1/2) was negative
    surface = 1

# show general results on screen
print('')
print('Total mass : {: 6.2f}'.format(m / m0))







