/**
 * @file clic.c
 * @author Daniel R. Tabin
 * @brief Functions for Composite Likelihood Information Criterion.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "strdblstck.h"

// Function prototypes
void usage(void);

//vars

const char *usageMsg =
     "usage: bepe <file.pts> <boot1.legofit> <boot2.legofit> ...\n"
     "  where file.pts is the .pts file produced by legofit with the real\n"
     "  data, and each \"boot\" file is the legofit output from one bootstrap\n"
     "  replicate. Must include .pts file and at least 2 boostrap files.\n";

int main(){
  // Command line arguments specify file names
  if(argc < 4)
      usage();

  const char *ptsfname = argv[1];
  int nfiles = argc-2;  // number of bootstrap files
  const char *bootfname[nfiles];
  for(i=0; i < nfiles; ++i)
      bootfname[i] = argv[i+2];

  // Read bootstrap files into an array of FIFO stacks
  StrDblStack *stack[nfiles];
  for(i=0; i < nfiles; ++i) {
      stack[i] = parseLegofit_BEPE(bootfname[i]);
      if(i>0) {
          if(StrDblStack_compare(stack[0], stack[i])) {
              fprintf(stderr, "%s:%d: inconsistent parameters in"
                      " files %s and %s\n", __FILE__,__LINE__,
                      bootfname[0], bootfname[i]);
              exit(EXIT_FAILURE);
          }
      }
  }
}
