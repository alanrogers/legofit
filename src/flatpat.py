#!/usr/bin/python

###
#@file flatpat.py
#@page flatpat
#@brief Read files of site pattern counts or frequencies and write a
#flat file with a column for each site pattern and a row for each
#input file.
#

from math import log
import sys

def usage():
    print "usage: flatpat.py inputfile1 [inputfile2 ...]"
    exit(1)

# Return a tuple containing a list of site pattern labels for the
# current file and a list of site pattern frequencies.
def process_file(fname):
    f = open(fname)

    lbl = []
    frq = []

    # Fill arrays of labels and frequencies
    reading_patterns = False
    for line in f:
        line = line.strip()
        if len(line)==0:
            continue
        line = line.strip().split()
        if not reading_patterns:
            if len(line) >= 2 and line[1] == "SitePat":
                reading_patterns = True
            continue
        if reading_patterns:
            lbl.append(line[0])
            currfrq = line[1]
            frq.append(currfrq)
            s += currfrq

    f.close()

    return (lbl, frq)
    
if len(sys.argv) < 2:
    usage()

# abort with usage message if there are any flag arguments    
for arg in sys.argv[1:]:
    if arg[0] == '-':
        usage()

# Print a line of output for each input file. Each line has a column
# for each site pattern. Make sure all files have the same site
# pattern labels.
first = True
for fname in sys.argv[1:]:
    lblvec, frqvec = process_file(fname)
    assert len(lblvec) == len(frqvec)
    if first:
        for lbl in lblvec:
            print lbl,
        print
        first = False
        lblvec0 = [lbl for lbl in lblvec]
    else:
        if lblvec != lblvec0:
            print >> sys.stderr, \
                "Pattern label mismatch in file \"%s\"" % fname
            sys.exit(1)
    for frq in frqvec:
        print frq,
    print


    


