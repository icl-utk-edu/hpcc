# -*- coding: utf-8 -*-
import sys

from xml.dom.minidom import parse

class XMLNODE:
  prefix_el = "HPCC:"

  def __init__(self, node):
    self.node = node

  def __getattr__(self, name):
    prfx = self.prefix_el
    if not name.startswith(prfx):
      name = prfx + name

    name = name.replace("_", "-")

    for n in self.node.childNodes:
      #print "N", name, n.nodeName, n.attributes, n.nodeValue, len(n.childNodes)
      if len(n.childNodes) > 0 and n.nodeName == name:
        return n.childNodes[0].nodeValue

class XML:
  site_el = "HPCC:Site"
  id_el = "HPCC:ID"

  def __init__(self, filename_or_file):
    self.dom = parse(filename_or_file)

  def __getitem__(self, idx):
    sidx = "%d" % idx
    if idx < 1 or idx > 500:
      raise ValueError, sidx
    for n in self.dom.childNodes[0].childNodes:
      if n.ELEMENT_NODE == n.nodeType and n.nodeName == self.site_el:
        name = self.id_el
        for nn in n.childNodes:
          if nn.nodeName == name:
            if nn.childNodes[0].nodeValue == sidx:
              return XMLNODE(n)

  def min_id(self): return self.minmax_id()[0]
  def max_id(self): return self.minmax_id()[1]
  def minmax_id(self):
    min_idx = 391
    max_idx = 1
    for n in self.dom.childNodes[0].childNodes:
      if n.ELEMENT_NODE == n.nodeType and n.nodeName == self.site_el:
        name = self.id_el
        for nn in n.childNodes:
          if nn.nodeName == name:
            idx = int(nn.childNodes[0].nodeValue)
            if idx < min_idx:
              min_idx = idx
            if idx > max_idx:
              max_idx = idx
    return min_idx, max_idx

def main(argv):
  fname = argv[1]
  d = XML(fname)
  for idx in range(d.min_id(), d.max_id()+1):
    nde = d[idx]
    if nde is None: continue
    print idx, nde.HPL, nde.SingleMPIProcessDGEMM, nde.HPLNodes

if "__main__" == __name__:
  sys.exit(main(sys.argv))
