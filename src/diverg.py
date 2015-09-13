#!/usr/bin/python
# diverg.py
#
# This program compares the output of two runs of the lego program and summarizes
# the difference between them using the Kullback-Leibler divergence.
from math import log
import sys

def usage():
    print "usage: diverg.py inputfile1 inputfile2"
    print "       Input file name \"-\" means standard input."
    print "       At most one input may be \"-\"."
    exit(1)

def openInput(fname):
    if fname == "-":
        return sys.stdin
    else:
        return open(fname, "r")

if len(sys.argv) != 3:
    usage()
fname1 = sys.argv[1]    
fname2 = sys.argv[2]
if fname1 == fname2 == "-":
    usage()
f1 = openInput(fname1)
f2 = openInput(fname2)

pat1 = []
prob1 = []
pat2 = []
prob2 = []

for line in f1:
    line = line.strip()
    if len(line)==0 or line[0] == '#':
        continue
    line = line.split()
    pat1.append(line[0])
    prob1.append(float(line[1]))

for line in f2:
    line = line.strip()
    if len(line)==0 or line[0] == '#':
        continue
    line = line.split()
    pat2.append(line[0])
    prob2.append(float(line[1]))

if len(pat1) != len(pat2):
    print "Input files must have the same number of patterns,"
    print "but file1 has %d and file2 has %d" % (len(pat1), len(pat2))
    exit(1)

# Find field widths
wid1 = max(9, len(fname1))
wid2 = max(9, len(fname2))
widpat=7
for i in range(len(pat1)):
    if pat1[i] != pat2[i]:
        print "Mismatch in %d'th pattern." % i
        print "  pat1=%s pat2=%s" % (pat1[i], pat2[i])
        exit(1)
    widpat = max(widpat, len(pat1[i]), len(pat2[i]))

fmt = "%%%ds %%%ds %%%ds %%6s" % (widpat, wid1, wid2)
print fmt % ("SitePat", fname1, fname2, "KL")
fmt = "%%%ds %%%d.5f %%%d.5f %%6.3f" % (widpat, wid1, wid2)

KLsum = 0.0
for i in range(len(prob1)):
    kl = prob2[i]*log(prob2[i]/prob1[i])  # Kullback-Leibler divergence
    KLsum += kl
    print fmt % (pat1[i], prob1[i], prob2[i], kl),
    if abs(kl) >= 0.01:
        print "*",
    print

fmt = "%%%ds %%%ds %%%ds %%9.6f" % (widpat, wid1, wid2)
print fmt % ("", "", "", KLsum),
    

