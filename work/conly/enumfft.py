
import sys
import math

"""
>python enumfft.py |grep Variable.w1\  |sort -n -k 6
"""

def factor_single(n, f):
  c = 0
  while n % f == 0:
    c += 1
    n /= f
  return c, n

def factor235(n):
  l = list()
  for f in (2, 3, 5):
    c, n = factor_single(n, f)
    l.append(c)
  return tuple(l)

def next235(t):
  first_time = 1

  i2 = t[0]
  i3 = t[1]
  i5 = t[2]
  while i5 <= 13:
    while i3 <= 19:
      while i2 <= 30:
        if first_time:
          first_time = 0

        elif i2 == 0 and i3 == 0 and i5 == 0:
          pass

        else:
          if 2**i2 * 3**i3 * 5**i5 > 1<<31:
            break

          return i2, i3, i5

        i2 += 1
      i2 = 0
      i3 += 1
    i3 = 0
    i5 += 1

  return (30, 19, 13)

def find_minsize(results, allbounds):
  allnames = ("w1", "w2", "ww", "ww2", "ww3", "ww4", "c", "d")

  powers235 = list()
  for p in range(1, 30):
    powers235.append(1<<p)
    if p <= 19: powers235.append(3**p)
    if p <= 13: powers235.append(5**p)
  powers235.sort()

  mm = dict()
  for v in powers235:
    pass

  for n in results:
    for name in allnames:
      if not results[n].has_key(name):
        raise ValueError, "%d %s" % (n, name)

      bound = results[n][name]
      fcn, mul, dflt = allbounds[name]
      alloc_size = max(int(mul*fcn(n)), 1)
      if bound >= alloc_size:
        print "Exceeded", name, n, bound, alloc_size, 1.0 / (fcn(n) / bound), mul, dflt, "Problem" * int(dflt < bound)
      if bound >= dflt:
        print "Variable", name, n, bound, alloc_size, 1.0 / (fcn(n) / bound), mul, dflt

NDA2 = 65536
NBLK = 16
NP=8

def enumfft(fname):
  results = dict()
  ctpl = (1, 0, 0)

  bounds = {
      "INPUT":  [float, 1.0, 0],
      "OUTPUT": [float, 1.0, 0],
      "w1":  [math.sqrt,  1.1,    NDA2/2],
      "w2":  [math.sqrt,  0.375,  NDA2/2],
      "ww":  [math.sqrt,  1.0,    NDA2],
      "ww2": [math.sqrt,  3.9,    NDA2],
      "ww3": [math.sqrt,  5.4773, NDA2],
      "ww4": [float,     1.0/256, NDA2+2**13], # there are instances where the space is O(N/81), O(N/162), O(N/243): they all need 72899
      "c":   [math.sqrt, 16.75,    NDA2*NBLK/2],
      "d":   [math.sqrt,  1.0,    NDA2],
      }

  maxmult = dict()
  maxmultall = dict()

  for fline in open(fname):
    if fline.find("[") < 0: # skip lines without brackets
      continue

    name, rest = fline.split("[")
    n, size = map(int, rest.split("]"))

    if not results.has_key(n):
      results[n] = dict()

    results[n][name] = size # structure of 'results' results[32768]["ww4"] = 127

    if "INPUT" == name:
      t = factor235(n)
      if t != ctpl:
        print "Error:", n, t, ctpl
      ctpl = next235(t)

    fnctn, mult, sizemin = bounds[name]
    if fnctn(n)*mult < size and size > sizemin:
      print name, n, size, fnctn(n)*mult, mult, sizemin, size / fnctn(n)

    maxmultall[name] = max(size / fnctn(n), maxmultall.get(name, 0.0))

    if size > sizemin:
      maxmult[name] = max(size / fnctn(n), maxmult.get(name, 0.0))

  ffmt = "%7.4f"
  for k in maxmultall:
    print "Multipliers:", "%6s" % k, ffmt % maxmult.get(k, -1.0), ffmt % bounds[k][1], ffmt % maxmultall[k]

  find_minsize(results, bounds)

def main(argv):
  enumfft("enumerate_all.txt")

  return 0

if "__main__" == __name__:
  sys.exit(main(sys.argv))
