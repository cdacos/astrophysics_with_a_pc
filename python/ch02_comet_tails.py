# -*- coding: utf-8 -*-

"""
Chapter 2 - Comet Tails

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-02)

Example:
$ python ch02_comet_tails.py .5 .95 1 .03
"""

from __future__ import print_function
from helpers import start_parameter
from math import cos, sin, sqrt, pi

print('Astrophysics with a PC : COMET TAILS')
print('------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input parameters and orbital elements : ')
ap  = start_parameter('Perihelion distance (A.U.)      : ', 1)
ecc = start_parameter('Eccentricity of the comet orbit : ', 2)
mu  = start_parameter('Parameter 1 - mu                : ', 3)
g   = start_parameter('Outflow velocity                : ', 4)
p = ap * (1 + ecc) # (4)

for i in range(-4, 5):
  # this for-loop considers 9 positions of the comet in its orbit
  # CLS

  # Compute position of the nucleus and parameters A1, A2 and A3
  nu = 0.5 * i
  r = p / (1 + ecc * cos(nu)) # (6)
  x = r * cos(nu) # (7)
  y = r * sin(nu) # (8)
  a1 = (sqrt(2) / sqrt(mu)) * r # (13a)
  a2 = (4 * ecc * r * sin(nu)) / (3 * mu * sqrt(p)) # (13b)
  a3 = (2 * sqrt(2 * p)) / (3 * r * sqrt(mu)) # (13c)

  print('\nPosition {: 2f} : true anomaly = {: 5.3f}    r = {: 9.6f}   x = {: 10.7f}   y = {: 10.7f}'.format(i + 5, nu, r, x, y))
  print('                     a1 = {: 9.6f}        a2 = {: 9.6f}        a3 = {: 9.6f}'.format(a1, a2, a3))
  print('')
  print('        s            t             x\'              y\'')

  for j in range(-1, 2):
    # this for-loop considers the 3 syndynames for each of the 9 positions
    gg = j * pi / 2
    ggdegree = j * 90
    print('Syndyname for G = {: 3.1f}'.format(ggdegree))

    for k in range(1, 10):
      # this for-loop computes 9 points (s,t) and their (x',y') transformation
      s = 0.05 * k
      t = g * sin(gg) * (a1 * sqrt(s) - a2 * s) + a3 * s * sqrt(s) # (14)
      xx = (s * x + t * y + r * x) / r # (11a)
      yy = (s * y - t * x + r * y) / r # (11b)
      #0if k % 2 == 1:
      print('     {: 5.2f}     {: 11.8f}     {: 11.7f}     {: 11.7f}'.format(s, t, xx, yy))

  # LOCATE 23, 1
  raw_input('Press enter key to continue')

