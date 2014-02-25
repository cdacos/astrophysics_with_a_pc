# -*- coding: utf-8 -*-

"""
Chapter 9 - Stellar Atmospheres

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-20)

Examples:
$ python ch09_stellar_atmosphere.py 10000 4 1.048
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow, log10, pi

def absorp(t, d):
  """ computes the absorption coefficient for given temperature T and density d
  """
  return 1.984e24 * d / pow(t, 3.5)

def dens(t, p, mu):
  """ computes the density for input values of temperature T, pressure P and mu
  """
  return p * mu / t / r

def disp(i, tau, t, pg, pr, rho, kk, ge, z):
  """ this is called each layer to display the results for that layer
  """
  print('{: 3d} {: 8.5f} {: 8.1f} {: 9.2f} {: 9.2f} {: 9.3f} {: 9.5f} {: 9.2f} {: 9.1f}'.format(i, tau, t, pg, pr, rho * 1e11, kk, ge, z / 100000.0))

def geffect(gs, teff, k, tau):
  """ computes the effective gravitational acceleration
  """
  return gs - k * a * pow(teff, 4) / 4.0 * (1 + .459 * exp(-3.4488 * tau))

def radpress(t):
  """ computes the radiation pressure for given temperature T and density d
  """
  return a / 3.0 * pow(t, 4)

def stap(teff, gs, mu, dtau, tau, t, pg, d, kk, ge, kk1, d1):
  """ computes the physical quantities at layer i+1 starting from layer i
  """
  tau1 = tau + .5 * dtau
  t1 = temp(tau1, teff)
  pg1 = pg + .5 * dtau * ge / kk
  d1 = dens(t1, pg1, mu)
  kk1 = absorp(t1, d1)
  ge1 = geffect(gs, teff, kk1, tau1)
  tau = tau + dtau
  t = temp(tau, teff)
  pg = pg + dtau * ge1 / kk1
  d = dens(t, pg, mu)
  kk = absorp(t, d)
  ge = geffect(gs, teff, kk, tau)
  return [teff, gs, mu, dtau, tau, t, pg, d, kk, ge, kk1, d1]

def temp(tau, teff):
  """ computes the temperature for given optical depth tau in effective temperature teff
  """
  q = .7104 - .1331 * exp(-3.4488 * tau)
  return teff * pow(.75 * (tau + q), .25)

# some physical constants used in the program :
g = 6.673e-8
a = 7.56464e-15
r = 8.314e7

print('Astrophysics with a PC : STELLAR MODEL')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input of initial conditions and model parameters : ')
teff = start_parameter('Effective temperature         : ', 1)
gs   = start_parameter('Surface gravit.accel.(log10)  : ', 2)
mu   = start_parameter('Average mean molecular weight : ', 3)
gs = pow(10, gs)

# displat heading of main table
print('')
print(' i     tau        T        Pg         Pr         d        k        ge       z(km)')

# compute central value of physical quantities
tau = .0
d = 1e-13
t = temp(tau, teff)
pg = r * d * t / mu
pr = radpress(t)
kk = absorp(t, d)
ge = geffect(gs, teff, kk, tau)
z = .0
i = 0
disp(i, tau, t, pg, pr, d, kk, ge, z)

# compute physical quantities in the layer 1
i = 1
dtau = .001
[teff, gs, mu, dtau, tau, t, pg, d, kk, ge, kk1, d1] = stap(teff, gs, mu, dtau, tau, t, pg, d, kk, ge, .0, .0)
pr = radpress(t) # carlos - this was missing?
z = dtau / kk1 / d1
disp(i, tau, t, pg, pr, d, kk, ge, z)

# next for-cycle computes the other layers
for i in range(2, 33):
  dtau = .25 * tau
  [teff, gs, mu, dtau, tau, t, pg, d, kk, ge, kk1, d1] = stap(teff, gs, mu, dtau, tau, t, pg, d, kk, ge, kk1, d1)
  dz = dtau / kk1 / d1
  pr = radpress(t)
  z = z + dz
  disp(i, tau, t, pg, pr, d, kk, ge, z)

  # pause after having shown layer 15
  if i == 15:
    raw_input('Press Enter to continue')

print('')
print('Model complete.')
