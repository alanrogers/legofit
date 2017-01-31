#!/usr/bin/python
###
#@file flatfile.py
#@page flatfile
#@brief Read legofit output files; write flat file of estimates
#
#flatfile.py: read legofit output files; write flat file of estimates
#====================================================================
#
#    usage: flatfile.py [options] <file1> <file2> ...
#
#    where <file*> files are legofit output files, which must each
#    estimate the same parameters
#    Options may include:
#
#      -h     or --help         Print this message.
#
#    The program writes to standard output.
#
#For example, the command
#
#    flatfile.py *.legofit
#
#would summarize all .legofit files in the curent directory.
#
# @copyright Copyright (c) 2017, Alan R. Rogers
# <rogers@anthro.utah.edu>. This file is released under the Internet
# Systems Consortium License, which can be found in file "LICENSE".

import sys
from math import floor, ceil
import datetime

# Print usage message and abort
def usage(msg1):
    if len(msg1) > 0:
        print >> sys.stderr, msg1
    msg = \
        """
usage: flatfile.py [options] <file1> <file2> ...

where the "file" arguments files are legofit output files
Options may include:

  -h     or --help         Print this message.

The program writes to standard output.
"""
    print >> sys.stderr, msg
    exit(1)

# Parse legofit output file.  Return a tuple containing two lists:
# first a list of parameter names; second a list of parameter
# estimates.
def parselegofit(fname):
    ifile = open(fname, "r")
    parmap = {}
    estmap = {}

    for line in ifile:
        if line[0] == "#":
            continue
        line = line.split("=")
        if len(line) != 2:
            continue

        key = line[0].strip()
        key = key.replace("2","two")
        if "Gaussian" in line[1]:
            value = 1.0
        else:
            value = float(line[1].strip())

        if key in parmap:
            estmap[key] = value
        else:
            parmap[key] = value

    ifile.close()

    parnames = sorted(estmap.keys())
    estimates = len(parnames)*[0.0]
    for i in range(len(parnames)):
        estimates[i] = estmap[parnames[i]]
    return (parnames, estimates)

fnames = []

# Loop over command line arguments, ignoring the 0th.
# (The 0th is just the name of the program.)
i = 1
while(True):
    if i >= len(sys.argv):
        break
    elif sys.argv[i]=="-h" or sys.argv[i]=="--help":
        usage("")
    elif sys.argv[i][0] == "-":
        usage("Unknown argument: %s" % sys.argv[i])
    else:
        fnames.append(sys.argv[i])
    i += 1

if len(fnames) < 1:
    usage("Command line must list at least 1 input file")

print "# flatfile.py run at: %s" % datetime.datetime.now()
print "# input files:",
for i in range(len(fnames)):
    print fnames[i],
print

mat = []
parnames = []
npar = 0

for name in fnames:
    parnames2, estimates = parselegofit(name)

    if npar > 0 and parnames != parnames2:
        print >> sys.stderr, "Input files estimate different parameters"
        print >> sys.stderr, "  1:", parnames
        print >> sys.stderr, "  2:", parnames2
        exit(1)

    parnames = parnames2
    npar = len(parnames)
    mat.append(estimates)

for name in parnames:
    print " %s" % name,
print

for row in mat:
    for val in row:
        print " %s" % val,
    print

