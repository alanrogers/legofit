#!/usr/bin/env python3
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
#    where <file*> files are legofit output files.
#    Options may include:
#
#      -t     or --transpose    Rows are parameters rather than data sets.
#      -h     or --help         Print this message.
#
#    The program writes to standard output.
#
#For example, suppose that s2.legofit contains estimates from the real
#data, and estimates from bootstrap replicates are in files with names
# s2boot0.legofit, s2boot1.legofit, and so on. Then the following command
#
#    flatfile.py s2.legofit s2boot*.legofit
#
#Parameters that are missing from a .legofit file will print as "None".
#
#The output begins with two lines of comment, which begin with a sharp
#character in column 1 and give (1) the date and time at which the
#program was run, and (2) the names of all the input files. The output
#then continues with (3) a row for each input file and a column for
#each parameter. Columns are separated by a single space character. To
#make the rows and columns pretty, try:
#
#    flatfile.py s2.legofit s2boot*.legofit | grep -v ^# | column -t
#
#The "grep" and "column" utilities are standard on Linux and OS X but
#not on Windows.

# @copyright Copyright (c) 2017, Alan R. Rogers
# <rogers@anthro.utah.edu>. This file is released under the Internet
# Systems Consortium License, which can be found in file "LICENSE".

import sys
from math import floor, ceil
import datetime
from os.path import basename

# Print usage message and abort
def usage(msg1):
    if len(msg1) > 0:
        print(msg1, file=sys.stderr)
    msg = \
        """
usage: flatfile.py [options] <file1> <file2> ...

where the "file" arguments files are legofit output files
Options may include:

  -t     or --transpose    Rows are parameters rather than data sets.
  -h     or --help         Print this message.

The program writes to standard output.
"""
    print(msg, file=sys.stderr)
    exit(1)

# Parse legofit output file.  Return a map relating parameter names to
# estimated parameter values.
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
        if "Gaussian" in line[1]:
            value = 1.0
        else:
            value = line[1].strip()

        # In legofit output, pairs are printed twice. First as initial
        # values, and second as estimates. This code adds a pair to
        # parmap the first time it is seen and to estmap the second
        # time. Thus, parmap will contain the initial values and
        # estmap the estimates.
        if key in parmap:
            estmap[key] = value
        else:
            parmap[key] = value

    ifile.close()
    return estmap

fnames = []
transpose = False

# Loop over command line arguments, ignoring the 0th.
i = 1
while(True):
    if i >= len(sys.argv):
        break
    elif sys.argv[i]=="-h" or sys.argv[i]=="--help":
        usage("")
    elif sys.argv[i]=="-t" or sys.argv[i]=="--transpose":
        transpose = True
    elif sys.argv[i][0] == "-":
        usage("Unknown argument: %s" % sys.argv[i])
    else:
        fnames.append(sys.argv[i])
    i += 1

if len(fnames) < 1:
    usage("Command line must list at least 1 input file")

print("# flatfile.py run at: %s" % datetime.datetime.now())
print("# input files:", end=' ')
for i in range(len(fnames)):
    print(basename(fnames[i]), end=' ')
print()

allmaps = []  # allmaps[i] is the dictionary for file i
allnames = set([]) # set of all parameter names

# Get allmaps and allnames.
for name in fnames:
    estmap = parselegofit(name)
    if len(estmap) == 0:
        print("ERR: file %s has no parameters" % name, file=sys.stderr)
        sys.exit(1)
    allmaps.append(estmap)
    allnames |= set(estmap.keys())

allnames = sorted(list(allnames))
npar = len(allnames)
nfile = len(fnames)

# mat is a rectangular matrix, with missing parameters set to None.
mat = nfile*[None]
for i in range(nfile):
    mat[i] = [None for j in range(npar)]
    k = 0
    for j in range(npar):
        # If parameter name is in map i, then put the corresponding
        # value into mat. Otherwise, it retains its initial value of
        # None.
        if allnames[j] in allmaps[i]:
            mat[i][j] = allmaps[i][allnames[j]]

if transpose:
    nrows = npar
    ncols = len(fnames)
    print("param", end=' ')
    for j in range(ncols):
        print(fnames[j], end=' ')
    print()
    for i in range(nrows):
        print(allnames[i], end=' ')
        for j in range(ncols):
            if mat[j][i] == None:
                print("NA", end=' ')
            else:
                print(mat[j][i], end=' ')
        print()
else:
    for name in allnames:
        print("%s" % name, end=' ')
    print()

    for row in mat:
        for val in row:
            if val == None:
                print("NA", end=' ')
            else:
                print("%s" % val, end=' ')
        print()

