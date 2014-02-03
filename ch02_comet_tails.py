from __future__ import print_function
import math

print('Astrophysics with a PC : COMET TAILS')
print('------------------------------------')
print('')
print('-------Minimal solution program--------')
print('')
print('Input parameters and orbital elements : ')
ap  = float(raw_input('Perihelion distance (A.U.)      : '))
ecc = float(raw_input('Eccentricity of the comet orbit : '))
mu  = float(raw_input('Parameter 1 - mu                : '))
g   = float(raw_input('Outflow velocity                : '))
p = ap * (1 + ecc) # (4)

for i in range(-4, 5):
  # this for-loop considers 9 positions of the comet in its orbit
  # CLS

  # Compute position of the nucleus and parameters A1, A2 and A3
  nu = 0.5 * i
  r = p / (1 + ecc * math.cos(nu)) # (6)
  x = r * math.cos(nu) # (7)
  y = r * math.sin(nu) # (8)
  a1 = (math.sqrt(2) / math.sqrt(mu)) * r # (13a)
  a2 = (4 * ecc * r * math.sin(nu)) / (3 * mu * math.sqrt(p)) # (13b)
  a3 = (2 * math.sqrt(2 * p)) / (3 * r * math.sqrt(mu)) # (13c)

  print('\nPosition {: 2f} : true anomaly = {: 2.3f}    r = {: 2.3f}   x = {: 2.3f}   y = {: 2.3f}'.format(i + 5, nu, r, x, y))
  print('                     a1 = {: 3.5f}        a2 = {: 3.5f}        a3 = {: 3.5f}'.format(a1, a2, a3))
  print('')
  print('         s            t            x\'           y\'')

  for j in range(-1, 2):
    # this for-loop considers the 3 syndynames for each of the 9 positions
    gg = j * math.pi / 2
    ggdegree = j * 90
    print('Syndyname for G = {: 3.1f}'.format(ggdegree))

    for k in range(1, 10):
      # this for-loop computes 9 points (s,t) and their (x',y') transformation
      s = 0.05 * k
      t = g * math.sin(gg) * (a1 * math.sqrt(s) - a2 * s) + a3 * s * math.sqrt(s) # (14)
      xx = (s * x + t * y + r * x) / r # (11a)
      yy = (s * y - t * x + r * y) / r # (11b)
      #0if k % 2 == 1:
      print('     {: 3.5f}     {: 3.5f}     {: 3.5f}     {: 3.5f}'.format(s, t, xx, yy))

  # LOCATE 23, 1
  raw_input('Press enter key to continue')

