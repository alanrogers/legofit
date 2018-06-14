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
    "usage: bepe -d <real_data> -D <bootdata_1> <bootdata_2> ..."
    " -L <boot1.legofit> <boot2.legofit> ...\n"
    "  where real_data is the real data,\n"
    "  each \"bootdata_\" file is the legofit output from one bootstrap"
    " replicate,\n"
    "  and each \"boot#.legofit\" file is the legofit output from one"
    " bootstrap replicate\n"
    "  Must include real_data file and at least 2 boostrap replicates.\n";

 void usage(void) {
     fputs(usageMsg, stderr);
     exit(EXIT_FAILURE);
 }

int main(int argc, char **argv){
  // Command line arguments specify file names
  if(argc < 7)
      usage();

  const char *realfName = argv[2];
  int nfiles = ((argc-5)/2);  // number of bootstrap files
  const char *legofname[nfiles];
  const char *datafname[nfiles];
  for(int i = 0; i < nfiles; ++i)
      datafname[i] = argv[i+4];
  for(int i = 0; i < nfiles; ++i)
      legofname[i] = argv[i+5+nfiles];

  // Read bootstrap files into an array of FIFO stacks
  StrDblStack* data_stack[nfiles];
  StrDblStack* lego_stack[nfiles];

  StrDblStack* real_stack = parseLegofit_BEPE(realfName);

  for(int i = 0; i < nfiles; ++i) {
      if(StrDblStack_compare(real_stack, lego_stack[i]) ||
         StrDblStack_compare(real_stack, data_stack[i])) {
          fprintf(stderr, "%s:%d: inconsistent parameters in"
                  " files\n", __FILE__,__LINE__);
          exit(EXIT_FAILURE);
      }
      lego_stack[i] = parseLegofit_BEPE(legofname[i]);
      data_stack[i] = parseLegofit_BEPE(datafname[i]);
  }
}
