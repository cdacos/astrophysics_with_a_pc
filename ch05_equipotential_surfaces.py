# -*- coding: utf-8 -*-

"""
Chapter 5 - Equipotential Surfaces of the Two-Body Problem

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-20)

Examples:
$ python ch05_equipotential_surfaces.py .4 1.05 0 1 .05
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow

def dvdx(x, y, mu):
  """evaluates partial derivative dV/dx in (x,y)
  """
  return -(1 - mu) * (x - mu) / pow(r1(x, y, mu), 3) - mu * (x + 1 - mu) / pow(r2(x, y, mu), 3) + x

def dvdy(x, y, mu):
  """evaluates partial derivative dV/dy in (x,y)
  """
  return -(1 - mu) * y / pow(r1(x, y, mu), 3) - mu * y / pow(r2(x, y, mu), 3) + y

def r1(x, y, mu):
  """computes distance from point (x,y) to first primary
  """
  return sqrt(pow(x - mu, 2) + pow(y, 2))

def r2(x, y, mu):
  """computes distance from point (x,y) to second primary
  """
  return sqrt(pow(x + 1 - mu, 2) + pow(y, 2))

def vxy(x, y, mu):
  """evaluates potential energy fucntion V(x,y) at point (x,y)
  """
  return (1 - mu) / r1(x, y, mu) + mu / r2(x, y, mu) + (pow(x, 2) + pow(y, 2)) / 2

print('Astrophysics with a PC : EQUIPOTENTIAL CURVES')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of initial conditions and parameters : ')
print('')
mu = start_parameter('Mass parameter mu         : ', 1)
x  = start_parameter('Initial conditions : x(0) : ', 2)
y  = start_parameter('                     y(0) : ', 3)

k = vxy(x, y, mu)
print('Potential constant K = {: 11.7f}'.format(k))

# Ask user in which direction to move
drc = start_parameter('Enter 1 to compute Eq.curve as y(x), 2 to compute Eq.curve as x(y)', 4)

if drc == '1':
  dx = start_parameter('stepsize dx = ', 5)
else:
  dy = start_parameter('stepsize dy = ', 5)

# here starts main cycle
stp = 0

while stp != 1:
  # use midpoint method to compute y as function of x if drc is 1
  if drc == '1':
    x12 = x + .5 * dx
    y12 = y - .5 * dx * dvdx(x, y, mu) / dvdy(x, y, mu)
    x = x + dx
    y = y - dx * dvdx(x12, y12, mu) / dvdy(x12, y12, mu)

    # add one correction of the computed value of y
    y = y - (vxy(x, y, mu) - k) / dvdy(x, y, mu)

    # use midpoint method to compute x as function of y if drc is 2
  else:
    y12 = y + .5 * dy
    x12 = x - .5 * dy * dvdy(x, y, mu) / dvdx(x, y, mu)
    y = y + dy
    x = x - dy * dvdy(x12, y12, mu) / dvdx(x12, y12, mu)

    # add one correction of the computed value of x
    x = x - (vxy(x, y, mu) - k) / dvdy(x, y, mu)

  # show results on screen
  print('')
  print('x = {: 11.7f}       y = {: 11.7f}      K = {: 11.7f}'.format(x, y, k))
  print('')

  # ask the user what to do next
  print('Enter c to change between x(y) and y(x)')
  print(' s to stop, or any key to continue with actual y(x) or y(x)')
  chopt = raw_input('')

  # if user has answered with S : prepare to stop the program
  if chopt == 's' or chopt == 'S':
    stp = 1

  # if user answered with C : change the direction and select new stepsize
  if chopt == 'c' or chopt == 'C':
    if drc == '1':
      drc = '2'
      dy = float(raw_input('stepsize dy = '))
    else:
      drc = '1'
      dx = float(raw_input('stepsize dx = '))

