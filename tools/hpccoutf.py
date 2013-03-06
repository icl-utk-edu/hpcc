#! /usr/bin/env python

import re
import sys

HPCC_out = dict(re.findall(r"^(\w+)=(\d.*)$", sys.stdin.read(), re.MULTILINE))

Walk_Order = (
  "HPL_Tflops",                      "PTRANS_GBs",
  "MPIRandomAccess_GUPs",            "MPIFFT_Gflops",
  "StarSTREAM_Triad*CommWorldProcs", "StarSTREAM_Triad",
  "StarDGEMM_Gflops",                "RandomlyOrderedRingBandwidth_GBytes",
  "RandomlyOrderedRingLatency_usec" )

Walk_Units = (
  "Tera Flops per Second",   "Tera Bytes per Second",
  "Giga Updates per Second", "Tera Flops per Second",
  "Tera Bytes per Second",   "Giga Bytes per Second",
  "Giga Flops per Second",   "Giga Bytes per second",
  "micro-seconds");

Cross_Walk = {
    "HPL_Tflops"           : "G-HPL",
    "PTRANS_GBs"           : "G-PTRANS",
    "MPIRandomAccess_GUPs" :  "G-RandomAccess",
    "MPIFFT_Gflops"        :  "G-FFT",
    "StarSTREAM_Triad*CommWorldProcs"     :  "EP-STREAM Sys",
    "CommWorldProcs"                      :  "MPI Processes",
    # StarSTREAM_Triad * CommWorldProcs   :   EP-STREAM Sys
    "StarSTREAM_Triad"                    :  "EP-STREAM Triad",
    "StarDGEMM_Gflops"                    :  "EP-DGEMM",
    "RandomlyOrderedRingBandwidth_GBytes" :  "RandomRing Bandwidth",
    "RandomlyOrderedRingLatency_usec"     :  "RandomRing Latency",
}

def show_all():
  for key in sorted(HPCC_out.keys()):
    print key +":", HPCC_out[key]

def show_web():
  count = 0
  for key in Walk_Order:
    if key == "StarSTREAM_Triad*CommWorldProcs":
      print key, Cross_Walk[key], float(HPCC_out["StarSTREAM_Triad"]) * int(HPCC_out["CommWorldProcs"]), Walk_Units[count]
    else:
      print key, Cross_Walk[key], HPCC_out[key], Walk_Units[count]
    count += 1

Show_all = 1
Show_web = 0

if Show_all:
    show_all()

if Show_web:
    show_web()
