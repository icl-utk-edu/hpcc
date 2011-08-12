
"""
Preprocess C files.

Change c_re() c_im() under some circumstances.
"""

import sys

def errlog(msg):
  sys.stderr.write(str(msg))
  sys.stderr.write("\n")

def proc_c_xx(oldline):
  l = list()
  idx = -1
  while 1:
    idx = oldline.find("c_", idx+1)
    if idx < 0:
      break

    if oldline[idx+2:].startswith("re") or oldline[idx+2:].startswith("im"):
      idx = oldline.find("(", idx)+1
      while oldline[idx].isspace():
        idx += 1

    l.append(idx)

  newline = oldline

  for idx in l:
    oparen = oldline.find("(", idx)
    cparen = oldline.find(")", idx)
    sqobrkt = oldline.find("[", idx)
    sqcbrkt = oldline.find("]", idx)

    if sqobrkt < 0 or sqobrkt > cparen: # '[' is not there or is beyond ')'
      continue

    if oparen >= 0 and oparen < sqobrkt: # if '(' is there and it's before '['
      continue

    newline = newline[:sqobrkt] + "->sqbracket(" + newline[sqobrkt+1:sqcbrkt] + ")" + newline[sqcbrkt+1:]

  return newline


def cpp(fname):
  for fline in open(fname):
    if fline.find("c_re") >= 0 or fline.find("c_im"):
      newline = proc_c_xx(fline)
    else:
      newline = fline
    print newline,

def main(argv):
  for a in argv[1:]:
    print "/****", a, "****/"
    cpp(a)

if "__main__" == __name__:
  sys.exit(main(sys.argv))
