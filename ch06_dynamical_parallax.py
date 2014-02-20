# -*- coding: utf-8 -*-

"""
Chapter 6 - The Dynamical Parallax

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-20)

Examples:
$ python ch06_dynamical_parallax.py 78.8 17.6 .3 1.7 .06 .3
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow

print('Astrophysics with a PC : DYNAMICAL PARALLAX')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of the observed data : ')
p   = start_parameter('Orbital period (years)                     : ', 1)
a   = start_parameter('Apparent distance (arc seconds)            : ', 2)
mv1 = start_parameter('Apparent magnitude of first component      : ', 3)
mv2 = start_parameter('Apparent magnitude of second component     : ', 4)
bc1 = start_parameter('Bolometric magnitude of first component    : ', 5)
bc2 = start_parameter('Bolometric magnitude of second component   : ', 6)
print('')
print('  i       m1        m2      dist       par       Mb1       Mb2')

eps = .01 # pre-defined constant in text

# select starting values for the two masses
m1 = 1.0
m2 = 1.0
stopcrit = 0
i = 1

m11 = 0.0
m22 = 0.0
dis = 0.0

# main cycle that stops when the two masses have converged
while stopcrit != 1 and i < 15:
  # compute new approximations for the two masses
  par = a / pow(p, 2 / 3.0) / pow(m1 + m2, 1 / 3.0)
  mabs1 = mv1 + 5 + 5 * log(par, 10.0)
  mabs2 = mv2 + 5 + 5 * log(par, 10.0)
  mb1 = mabs1 - bc1
  mb2 = mabs2 - bc2
  m11 = pow(10, .58 - .112 * mb1) # = new approximation of first mass
  m22 = pow(10, .58 - .112 * mb2) # = new approximation of second mass
  dis = 1 / par * 3.26            # = new approximation of the distance

  # show iteration on screen
  print('{: 3d} {: 9.2f} {: 9.2f} {: 9.3f} {: 9.2f} {: 9.2f} {: 9.2f}'.format(i, m1, m2, par, dis, mb1, mb2))

  # check convergence (if yes, stopcrit becomes 1)
  if abs(m1 - m11) < eps and abs(m2 - m22) < eps:
    stopcrit = 1
  else:
    m1 = m11
    m2 = m22
    i = i + 1

# show final results on screen
if i > 15:
  print('Method does not converge for your input')
else:
  print('Final results : mass of 1st component : {: 6.2f}'.format(m11))
  print('                mass of 2nd component : {: 6.2f}'.format(m22))
  print('              distance in light years : {: 6.2f}'.format(dis))



