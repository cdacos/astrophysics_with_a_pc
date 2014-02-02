from __future__ import print_function
import math

print('Astrophysics with a PC : COMET TAILS')
print('------------------------------------')
print('\n')
print('-------Minimal solution program--------')
print('\n')
print('Input parameters and orbital elements : ')
ap  = float(raw_input('Perihelion distance (A.U.)      : '))
ecc = float(raw_input('Eccentricity of the comet orbit : '))
mu  = float(raw_input('Parameter 1 - mu                : '))
g   = float(raw_input('Outflow velocity                : '))
p = ap * (1 + ecc)

for i in range(-4, 4):
  # this for-loop considers 9 positions of the comet in its orbit
  # CLS

  # Compute position of the nucleus and parameters A1, A2 and A3
  nu = 0.5 * i
  r = p / (1 + ecc + math.cos(nu))
  x = r * math.cos(nu)
  y = r * math.sin(nu)
  a1 = math.sqrt(2 / mu) * r
  a2 = 4 * ecc * r * math.sin(nu) / 3 / mu / math.sqrt(p)
  a3 = math.sqrt(8 * p / mu) / 3 / r

  print('Position {:2f}   : true anomaly = {:2.3f}    r = {:2.3f}   x = {:2.3f}   t = {:2.3f}'.format(i + 5, nu, r, x, y))
  print('                   a1 = {:3.5f}        a2 = {:3.5f}        a3 = {:3.5f}'.format(a1, a2, a3))
  print('\n')
  print('                     s                    t                  x                y')

  for j in range(-1, 1):
    # this for-loop considers the 3 syndynames for each of the 9 positions
    gg = j * math.pi / 2
    ggdegree = j * 90
    print('Syndyname for G = {:3.1f}'.format(ggdegree))

    for k in range(1, 9):
      # this for-loop computes 9 points (s,t) and their (x',y') transformation
      s = 0.05 * k
      t = g * math.sin(gg) * (a1 * math.sqrt(s) - a2 * s) + a3 * s * math.sqrt(s)
      xx = (s * x + t * y + r * x) / r
      yy = (s * y - t * x + r * y) / r
      if k % 2 == 1:
        print('     {:3.5f}     {:3.5f}     {:3.5f}     {:3.5f}'.format(s, t, xx, yy))

  # LOCATE 23, 1
  raw_input('Press enter key to continue')

