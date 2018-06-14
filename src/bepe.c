/**
 * @file clic.c
 * @author Daniel R. Tabin
 * @brief Functions for Composite Likelihood Information Criterion.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "strdblstck.h"
#include <math.h>

// Function prototypes
void usage(void);

//vars

const char *usageMsg =
    "usage: bepe <real_data> <bootdata_1> <bootdata_2> ..."
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
  if(argc < 5)
      usage();

  const char *realfName = argv[1];
  int nfiles = ((argc-2)/2);  // number of bootstrap files
  const char *legofname[nfiles];
  const char *datafname[nfiles];

  for(int i = 0; i < nfiles; ++i)
      datafname[i] = argv[i+2];
  for(int i = 0; i < nfiles; ++i)
      legofname[i] = argv[i+3+nfiles];

  // Read bootstrap files into an arrays of FIFO queues
  StrDblStack* data_stack[nfiles];
  StrDblStack* lego_stack[nfiles];

  StrDblStack* real_stack = parseLegofit_BEPE(realfName);


  for(int i = 0; i < nfiles; ++i) {
      lego_stack[i] = parseLegofit_BEPE(legofname[i]);
      data_stack[i] = parseLegofit_BEPE(datafname[i]);

      if(StrDblStack_compare(real_stack, lego_stack[i]) ||
         StrDblStack_compare(real_stack, data_stack[i])) {
          fprintf(stderr, "%s:%d: inconsistent parameters in"
                  " files\n", __FILE__,__LINE__);
          exit(EXIT_FAILURE);
      }
  }

  //normalize the queues
  real_stack = normalize(real_stack);
  for(int i = 0; i < nfiles; ++i) {
      lego_stack[i] = normalize(lego_stack[i]);
      data_stack[i] = normalize(data_stack[i]);
  }

  //find MSD
  double real_msd = 0;
  double boot_msd = 0;
  double msd;

  StrDblStack* temp_L;
  StrDblStack* temp_D;
  StrDblStack* temp_d = real_stack;

  for (int i = 0; i < nfiles; ++i){
    temp_L = lego_stack[i];
    for (int j = 0; j < nfiles; ++j){
      temp_D = data_stack[i];
      for (int k = 0; k < nfiles; ++k){
        if(j == 0){
          real_msd += pow((temp_L->strdbl.val
                          - temp_d->strdbl.val),2);
          temp_d = temp_d->next;
          printf("%lf\n", real_msd);
        }
        boot_msd += pow((temp_L->strdbl.val
                        - temp_D->strdbl.val),2);
        temp_D = temp_D->next;
        printf("%lf\n", boot_msd);
      }
    }
  }

  real_msd = (real_msd / nfiles);
  boot_msd = (boot_msd / (nfiles * nfiles));

  msd = real_msd + boot_msd;

  printf("MSD = %lf\n", msd);
}
