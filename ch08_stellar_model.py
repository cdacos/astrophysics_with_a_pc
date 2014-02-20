# -*- coding: utf-8 -*-

"""
Chapter 8 - Homogeneous Stellar Models

'Astrophysics with a PC' by Paul Hellings, ISBN 943396-43-3
Copyright (c) 1994 Paul Hellings. All rights reserved.

Python version of the book's QuickBasic source
by Carlos da Costa https://github.com/cdacos (2014-02-20)

Examples:
$ python ch07_polytropes.py 1.5 .05 2 3
"""

from __future__ import print_function, division
from helpers import start_parameter
from math import fabs, exp, sqrt, log, pow, log10, pi

def e(d, t):
  """ computes the energy production for given density d and temperature t
  """
  print(t)
  tt = exp(1 / 3.0 * log(t / 1e9))
  p1 = 1 + tt * ( .133 + tt * (1.09 + tt * .938))
  p2 = 1 + tt * (.027 + tt * (-.788 + tt * (-.149 + tt * (.261 + tt * .127))))
  e1 = 23760.0 / pow(tt, 2) * p1 * exp(-3.38 / tt)
  e2 = 8.6665e25 / pow(tt, 2) * p2 * exp(-15.228 / tt - pow(tt, 6) / 9.5481)
  return d * (pow(xx, 2) * e1 + .02 * xx * e2)

def poly(n, x, f, h, pc, dc, rn):
  """ computes pressure, density, mass and distance to the centre starting
      from the polytrope results (x, f and h) and parameters n, pc and dc
  """
  p = pc * pow(f, n + 1)
  d = dc * pow(f, n)
  mr = -4 * pi * dc * pow(rn, 3) * pow(x, 2) * h
  r = rn * x
  return [p, d, mr, r]

def temp(mu, p, d):
  """ solves the equation of state to compute the temperature from the 
      pressure, density and mean molecular weight
  """
  tt = mu * p / rgas / d
  for i in range(1, 11):
    tt = mu / rgas / d * (p - 1 / 3.0 * a * pow(tt, 4))
  return tt

# declare a number of constant used in the program
g = 6.673e-8
a = 7.56464e-15
rgas = 8.314e7
xx = .7
yy = .27
zz = .03
mu = .618238
m0 = 2e33
r0 = 6.96e10
l0 = 3.83e33

print('Astrophysics with a PC : STELLAR MODEL')
print('--------------------------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
mtot = start_parameter('Input parameter : Approximation; of; total; mass(2 - 15) : ', 1)
print('')

# show heading of main table
print('i    Mr/Mo   log(p)   log(T)    log(d)   r/r0   log(E)    log(L)   x   f   h')

w = log10(mtot)

# central values of temperature and density
tc = 7.23937 + .2724354 * w - .0401771 * w * w
tc = pow(10, tc)
dc = 2.27899 - 1.658707 * w + .29329095 * w * w

# compute the value of the fitmass :
if mtot < 4:
  ffit = 9.0
else:
  if mtot < 10: # TODO - check if "m" was a bug in original text
    ffit = 19.58794 - 17.58794 * w
  else:
    ffit = 2.0

# compute other central quantities
pgc = rgas * dc * tc / mu
prc = 1 / 3.0 * a * pow(tc, 4)
ptc = pgc + prc
betac = pgc / ptc
beta = 1 - 2 / 3.0 * (1 - betac)
fb = (8 - 6 * beta) / (32 - 24 * beta - 3 * pow(beta, 2))
n = (1 - fb) / fb
rn = sqrt((n + 1) * ptc / 4.0 / pi / g / dc / dc)
eec = e(dc, tc)
i = 0
x = .0
f = 1.0
h = .0
m = .0
r = .0
l = 10.0
[p, d, m, r] = poly(n, x, f, h, ptc, dc, rn)

# print first line of table (central value)
print('{:2d} {: 8.5f} {: 7.4f} {: 7.4f} {: 8.4f} {: 7.4f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f}'.format(i, m / m0, log10(p), log10(tc), log10(d), r / r0, log10(eec), log10(1 / 10.0), x, f, h))

# compute and show first step
i = 1
l = 0
dx = .1
x = dx
f = 1 - 1 / 6.0 * pow(dx, 2) + n / 120.0 * pow(dx, 4)
h = -1 / 3.0 * dx + n / 30.0 * pow(dx, 3)
[p, d, m, r] = poly(n, x, f, h, ptc, dc, rn)
t = temp(mu, p, d)
ee = e(d, t)
dr = dx * rn
l = 4 / 3.0 * pi * dc * eec * pow(dr, 3)
print('{:2d} {: 8.5f} {: 7.4f} {: 7.4f} {: 8.4f} {: 7.4f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f}'.format(i, m / m0, log10(p), log10(tc), log10(d), r / r0, log10(eec), log10(1 / 10.0), x, f, h))

# Start of main cycle

i = 2
for zone in range(1, 3):  # zone = 1 during convective region
                          #      = 2 during radiative region
  stp = 0
  while stp != 1:
    flast = f # Save previous value of polytrope variable in flast
    x12 = x + .5 * dx
    f12 = f + .5 * dx * h

    # check whether f12 f12 is still positive
    if f12 > 0:
      h12 = h + .5 * dx * (-pow(f, n) - 2 * h / x)
      p12 = ptc * pow(f12, n + 1)
      d12 = dc * pow(f12, n)
      t12 = temp(mu, p12, d12)

      # compute energy production when in convective zone
      if zone == 1:
        ee12 = e(d12, t12)
      x = x + dx
      f = f + dx * h12

      # check whether f is still positive
      if f > 0:
        h = h + dx * (-pow(f12, n) - 2 * h12 / x12)
        p = ptc * pow(f, n + 1)
        d = dc * pow(f, n)
        m = -4 * pi * dc * pow(rn, 3) * pow(x, 2) * h
        r = rn * x
        t = temp(mu, p, d)

        # compute new value of luminosity when in convective zone
        if zone == 1:
          l = l + 4 * pi * d12 * dx * ee12 * pow(rn, 3) * pow(x12, 2)

        # compute energy production for new state when in convective
        # zone. In radiative zone, ee is put to 1, so that log(ee) becomes 0
        if zone == 1:
          ee = e(d, t)
        else:
          ee = 1

        # show results of layer i+1 on screen
        print('{:2d} {: 8.5f} {: 7.4f} {: 7.4f} {: 8.4f} {: 7.4f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f} {: 7.3f}'.format(i, m / m0, log10(p), log10(tc), log10(d), r / r0, log10(eec), log10(1 / 10.0), x, f, h))
      else:
        surface = 1 # surface has been reached
    else:
      surface = 1 # surface has been reached

    # check if in convective zone
    if zone == 1:
      test = 1.339944e9 * p / m * l / pow(t, 4) / fb

      # check if boundary of convective zone is reached
      if test < 1:
        raw_input('Boundary of convective core is reached. Press Enter')

        # compute fitting parameters
        f = ffit
        ptc = p / pow(ffit, 4)
        dc = d / pow(ffit, 3)
        rn = sqrt(ptc / pi / g / pow(dc, 2))
        x = r / rn
        h = -m / 4.0 / pi / dc / pow(rn, 3) / pow(x, 2)
        n = 3
        dx = .04
        stp = 1
        surface = 0

        # show fitting parameters
        print('')
        print('---------- Fitting parameters -----------')
        print('New centr.press.    : {:8.5f}'.format(ptc))
        print('New centr.density   : {:8.5f}'.format(dc))
        print('xfit                : {:8.5f}'.format(x))
        print('ffit                : {:8.5f}'.format(ffit))
        print('hfit                : {:8.5f}'.format(h))
        print('new rn              : {:8.5f}'.format(rn))
        print('')

        raw_input('Press Enter to continue')
        print('')

      # next else is entered when in radiative zone
      else:
        dx = 1.1 * dx

        # check if surface is reached
        if surface == 1:

          # compute exact location of surface and surface data
          stp = 1
          xs = x - flast / h
          mt = m + .5 * pi * d * pow(rn, 3) * pow(x + xs, 2)
          rad = r + rn * (xs - x)
          logteff = 3.7613 + .25 * log10(1 / 10.0) - .5 * log10(rad / r0)

          raw_input('Press Enter to continue')
          print('')

          # print surface data
          print('----------------- surface data -----------------')
          print('Mass (in Mo)           : {: 7.2f}'.format(mt / m0))
          print('radius (Ro)            : {: 7.2f}'.format(rad / ro))
          print('Luminosity (log(L/Lo)) : {: 7.2f}'.format(log10(1 / 10.0)))
          print('Effect.temp.(log)      : {: 7.2f}'.format(logteff))





  



