# -*- coding: utf-8 -*-

"""
Chapter 13 - Cosmological Models ofr the Universe

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-21)

Examples:
$ python ch13_universe_model.py .35 .35
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow, log10, pi, sin, tan, cos

print('Astrophysics with a PC : UNIVERSE MODEL')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of of the parameters : ')
s = start_parameter('sigma(0) : ', 1)
q = start_parameter('q(0)     : ', 2)

# this for omputes the past (era=0) and the future (era=1)
for era in [0, 1]:
  x = .0
  y = 1.0
  z = 1.0
  ctn = 1

  print('')
  # select the time step
  if era == 0:
    dx = -.02
    print('Computations for the past')
  else:
    dx = .02
    print('Computations for the future')

  # display heading of table of results
  print('')
  print('    x          y           z')
  count=0
  print('{: 9.2f} {: 11.7f} {: 11.7f}'.format(x, y, z))

  # main cycle of the iterative procedure
  while ctn != 0:

    # decrease time step if function y is too steep
    if abs(z) > 2:
      dx = .01 * dx / abs(dx)

    # results at half step
    x12 = x + .5 * dx
    y12 = y + .5 * dx * z

    # check if scale factor y12 is still positive
    if y12 > 0:
      z12 = z + .5 * dx * (-s / y / y + (s - q) * y)

      # results for the full step
      x = x + dx
      y = y + dx * z12

      # check if scale factor y1 is still positive
      if y > 0:
        z = z + dx * (-s / y12 / y12 + (s - q) * y12)
        count = count + 1

        # show results of newly computed iteration on screen
        print('{: 9.2f} {: 11.7f} {: 11.7f}'.format(x, y, z))

        # offer the user the possibility to stop the actual era
        # every 15 iterations
        if (count + 1) % 15 == 0:
          answer = ''
          while answer != 'y' and answer != 'n':
            answer = raw_input('Continue y/n?')
          if answer == 'y':
            ctn = 1
          else:
            ctn = 0

      else: # This else is reached when Y1 was negative
        if era == 0:
          print('Model starts from Big Bang')
        else:
          print('Model ends in Collapse')

        raw_input('Press Enter to continue')


