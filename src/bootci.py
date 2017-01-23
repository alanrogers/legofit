#!/usr/bin/python
###
#@file bootci.py
#@page bootci
#@brief Calculate confidence interval from legofit output files.
#
#bootci.py: calculate confidence interval from legofit output
#============================================================
#
#    usage: bootci.py [options] <central_tendency> <boot1> <boot2> ...
#
#    where <central_tendency> is the legofit output for the real data
#    and the <boot*> files are legofit output for bootstrap replicates.
#    Options may include:
#
#      -c <x> or --conf <x>     Set confidence level.
#      -h     or --help         Print this message.
#      -l <x> or --label <x>    Set label.
#
#    The program writes to standard output.
#
#For example, suppose that we have run legofit against the real data
#to produce a file called lancre.legofit, and we've run it against 50
#bootstrap replicates to produce files boot1.legofit, boot2.legofit,
#and so on. Then the command
#
#    bootci.py lancre.legofit boot*.legofit
#
#would read these files and print output like this
#
#    # bootci.py run at: 2017-01-06 22:32:58.916983
#    # real data: lancre.legofit
#    # bootstrap replicates: boot0.legofit boot1.legofit boot10.legofit ...
#    # confidence: 0.950
#          par             est             low            high
#          Tnd  25882.18097500  25783.98581875  25918.98318675
#          Txy   1910.32458600   1898.67020375   2020.23384350
#           mD      0.00030900      0.00001650      0.00133550
#           mN      0.01769600      0.01396675      0.01715075
#       twoNnd   1204.62497100   1168.42503700   1327.14358450
#       twoNxy  59669.81112500  59120.53155525  61060.16095975
#     twoNxynd  34355.15657100  34093.55639600  34552.91827875
#
#Here, "par" lists the parameter names, "est" gives the estimates, and
#"low" and "high" give the lower and upper confidence bounds.
#
#By default, the confidence level is 95%, but this can be changed with
#the --conf command-line option.
#
#The --label <x> option will add a column in which every row contains
#the character string <x>. This is useful in data analysis, if one
#concatenates the output from several runs of bootci.py. The column of
#labels makes it possible to distinguish the estimates from different
#runs of bootci.py.
#
#For example,
#
#    bootci.py --label lancre lancre.legofit boot*.legofit
#
#would generate
#
#    # bootci.py run at: 2017-01-06 22:32:58.916983
#    # real data: lancre.legofit
#    # bootstrap replicates: boot0.legofit boot1.legofit boot10.legofit ...
#    # confidence: 0.950
#          par             est             low            high lbl
#          Tnd  25882.18097500  25783.98581875  25918.98318675 lancre
#          Txy   1910.32458600   1898.67020375   2020.23384350 lancre
#           mD      0.00030900      0.00001650      0.00133550 lancre
#           mN      0.01769600      0.01396675      0.01715075 lancre
#       twoNnd   1204.62497100   1168.42503700   1327.14358450 lancre
#       twoNxy  59669.81112500  59120.53155525  61060.16095975 lancre
#     twoNxynd  34355.15657100  34093.55639600  34552.91827875 lancre
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
usage: bootci.py [options] <central_tendency> <boot1> <boot2> ...

where <central_tendency> is the legofit output for the real data
and the <boot*> files are legofit output for bootstrap replicates.
Options may include:

  -c <x> or --conf <x>     Set confidence level.
  -h     or --help         Print this message.
  -l <x> or --label <x>    Set label.

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

# Interpolate in order to approximate the value v[p*(len-1)].  Raise
# exception if len(v)==0.
def interpolate(p, v):
    if len(v) == 0:
        raise Exception("interpolate can't deal with 0-length vector")

    goal = p*(len(v) - 1)
    i = int(floor(goal))
    j = int(ceil(goal))
    if i == j:
        return v[i]
    w = goal - i
    return (1.0-w)*v[i] + w*v[j]

conf = 0.95
realdata = None
bootnames = []
lbl = None

# Loop over command line arguments, ignoring the 0th.
# (The 0th is just the name of the program.)
i = 1
while(True):
    if i >= len(sys.argv):
        break
    elif sys.argv[i]=="-h" or sys.argv[i]=="--help":
        usage("")
    elif sys.argv[i]=="-c" or sys.argv[i]=="--confidence":
        i += 1
        if i >= len(sys.argv):
            usage("Missing arg to -c or --confidence")
        conf = float(sys.argv[i])
    elif sys.argv[i]=="-l" or sys.argv[i]=="--label":
        i += 1
        if i >= len(sys.argv):
            usage("Missing arg to -l or --label")
        lbl = sys.argv[i]
    elif sys.argv[i][0] == "-":
        usage("Unknown argument: %s" % sys.argv[i])
    else:
        if realdata == None:
            realdata = sys.argv[i]
        else:
            bootnames.append(sys.argv[i])
    i += 1

if len(bootnames) < 2:
    usage("Command line must list at least 3 input files")

print "# bootci.py run at: %s" % datetime.datetime.now()
print "# real data:", realdata
print "# bootstrap replicates:",
for i in range(len(bootnames)):
    print bootnames[i],
print
print "# %s: %0.3f" % ("confidence", conf)

parnames, realEst = parselegofit(realdata)
mat = []
mat.append(realEst)

for name in bootnames:
    parnames2, estimates = parselegofit(name)

    if parnames != parnames2:
        print >> sys.stderr, "Input files estimate different parameters"
        print >> sys.stderr, "  1:", parnames
        print >> sys.stderr, "  2:", parnames2
        exit(1)

    mat.append(estimates)

# transpose matrix, so that mat[i] is a list of estimates for a single
# parameter.
mat = zip(*mat)
assert len(mat) == len(parnames)

# number of parameters
npar = len(mat)

tailProb = (1.0 - conf)/2.0

if lbl:
    print "%10s %15s %15s %15s %s" % ("par", "est", "low", "high", "lbl")
else:
    print "%10s %15s %15s %15s" % ("par", "est", "low", "high")
for i in range(npar):
    v = sorted(mat[i])
    lowBnd = interpolate(tailProb, v)
    highBnd = interpolate(1.0-tailProb, v)
    if lbl:
        print "%10s %15.8f %15.8f %15.8f %s" % \
            (parnames[i], realEst[i], lowBnd, highBnd, lbl)
    else:
        print "%10s %15.8f %15.8f %15.8f" % \
            (parnames[i], realEst[i], lowBnd, highBnd)
