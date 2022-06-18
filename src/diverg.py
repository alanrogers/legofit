#!/usr/bin/env python

###
#@file diverg.py
#@page diverg
#@brief Compare two sets of site-pattern counts or frequencies
#
#diverg.py: compare two sets of site-pattern counts or frequencies
#=================================================================
#
#This program compares two sets of site-pattern counts or frequencies
#and summarizes the difference between them using the Kullback-Leibler
#divergence.
#
#    usage: diverg.py inputfile1 [inputfile2]
#           Input file name "-" means standard input.
#           At most one input may be "-".
#
#Each input file should look like:
#
#    # SitePat           Count
#          x:y          227852
#          x:n             663
#          y:n           17410
#
#Here, the first line is an optional comment, which is used to label
#the columns. The first column contains site-pattern labels,
#which consist of population labels separated by colons. The second
#column contains the contribution of each site pattern. These need not
#be integers. They will be normalized in the output so that they sum
#to unity.
#
#This input format is also the output format of the program @ref legosim
#"legosim". This makes it possible to pipe legosim output into
#diverg.py. For example, suppose we have a file called "input.lgo" in
#@ref lgo "lgo" format, and a file called "sitepat.txt" containing
#observed site-pattern counts such as those above. Then the command
#
#    legosim -i 1000 input.lgo | diverg.py sitepat.txt -
#
#would compare the data in sitepat.txt to the output of the legosim
#command. This would produce something like the following:
#
#    SitePat sitepat.txt         -     KL
#        x:y     0.92651   0.87006 -0.055 *
#        x:n     0.00270   0.00671  0.006
#        y:n     0.07079   0.12323  0.068 *
#                                   0.019732
#
#The "sitepat.txt" column contains the data from file sitepat.txt,
#re-expressed as relative frequencies. The "-" column summarizes the
#legosim output in the same fashion. The KL column lists contributions to
#the Kulback-Leiebler (KL) divergence. When these contributions are
#small, they are approximately equal to the difference between columns
#3 and 2. Large KL contributions are marked with an asterisk. The
#bottom entry in this column is the KL divergence, a measure of the
#discrepancy between the two frequency distributions.

from math import log
import sys

def usage():
    print "usage: diverg.py inputfile1 [inputfile2]"
    print "       Input file name \"-\" means standard input."
    print "       At most one input may be \"-\"."
    print "       If only 1 file is given, print normalized frequencies"
    exit(1)

def openInput(fname):
    if fname == "-":
        return sys.stdin
    else:
        return open(fname, "r")

dokl = False
fname2 = None

if len(sys.argv) == 3:
    dokl = True
    fname2 = sys.argv[2]
elif len(sys.argv) != 2:
    usage()
fname1 = sys.argv[1]
if fname1 == fname2 == "-":
    usage()
f1 = openInput(fname1)
if dokl:
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

if dokl:    
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
if dokl:
    wid2 = max(9, len(fname2))
else:
    wid2=0
widpat=7
if dokl:
    for i in range(len(pat1)):
        if pat1[i] != pat2[i]:
            print "Mismatch in %d'th pattern." % i
            print "  pat1=%s pat2=%s" % (pat1[i], pat2[i])
            exit(1)
        widpat = max(widpat, len(pat1[i]), len(pat2[i]))
else:
    for i in range(len(pat1)):
        widpat = max(widpat, len(pat1[i]))

# Normalize probabilities
s1 = sum(prob1)
prob1 = [z/s1 for z in prob1]
if dokl:
    s2 = sum(prob2)
    prob2 = [z/s2 for z in prob2]

if not dokl:
    fmt = "%%%ds %%%ds" % (widpat, wid1)
    print fmt % ("SitePat", fname1)
    fmt = "%%%ds %%%d.5f" % (widpat, wid1)
    for i in range(len(prob1)):
        print fmt % (pat1[i], prob1[i])
    exit(0)

fmt = "%%%ds %%%ds %%%ds %%7s" % (widpat, wid1, wid2)
print fmt % ("SitePat", fname1, fname2, "KL")
fmt = "%%%ds %%%d.5f %%%d.5f %%7.4f" % (widpat, wid1, wid2)

KLsum = 0.0
for i in range(len(prob1)):
    kl = prob2[i]*log(prob2[i]/prob1[i])  # Kullback-Leibler divergence
    KLsum += kl
    print fmt % (pat1[i], prob1[i], prob2[i], kl),
    if abs(kl) >= 0.01:
        print "*",
    print

fmt = "%%%ds %%%ds %%%ds %%9.6g" % (widpat, wid1, wid2)
print fmt % ("", "", "", KLsum),


