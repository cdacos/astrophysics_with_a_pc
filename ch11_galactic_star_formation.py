# -*- coding: utf-8 -*-

"""
Chapter 11 - Star Formation in the Galaxy

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-21)

Examples:
$ python ch11_galactic_star_formation.py .15 .1 1 10 2.5
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow, log10, pi

def effes(a, m, k1, k2, n):
  """ computes the right hand sides of the differential equations for the
      variables a (atomic fraction) and m (molecular fraction)
  """
  fa = 1 - a - m - k1 * pow(m, 2) * a
  if m > .0:
    fm = k1 * pow(m, 2) * a + k2 * (a - 1 + m) * pow(m, n)
  else:
    fm = 0
  return [fa, fm]

print('Astrophysics with a PC : GALACTIC STAR FORMATION')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of initial conditions and parameters : ')
m  = start_parameter('Initial fraction of molecular clouds : ', 1)
a  = start_parameter('Initial fraction of atomic gas       : ', 2)
n  = start_parameter('Parameter n  : ', 3)
k1 = start_parameter('Parameter k1 : ', 4)
k2 = start_parameter('Parameter k2 : ', 5)
print('')

s = 1 - a - m
dx = .02

# heading of main cycle
print('  i      x         a         m           s')

# prepare for main cycle
x = 0
ni = 1
print('{: 3d} {: 6.2f} {: 10.4f} {: 10.4f} {: 10.4f}'.format(0, x, a, m, s))

# here start main cycle, computing blocks of 20 iterations
ch = ''
while ch != 's' and ch != 'S':
  for i in range(1, 20):

    # values at half the step (i+1/2)
    [fa, fm] = effes(a, m, k1, k2, n)
    m1 = m + .5 * dx * fm
    a1 = a + .5 * dx * fa
    s1 = 1 - m1 - a1

    # values at the full step (i+1)
    [fa, fm] = effes(a1, m1, k1, k2, n)
    m = m + dx * fm
    a = a + dx * fa
    s = 1 - m - a
    x = x + dx

    # show results of newly computed layer on screen
    print('{: 3d} {: 6.2f} {: 10.4f} {: 10.4f} {: 10.4f}'.format(ni, x, a, m, s))

    ni = ni + 1

  # Ask user whether to compute 20 iterations more or not
  ch = raw_input('Enter S to stop or any other value to continue')



