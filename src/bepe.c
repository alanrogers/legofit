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
    "usage: bepe <real_data> <bootdata_1> <bootdata_2> ..."
    "  -L <boot1.legofit> <boot2.legofit> ...\n"
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
  int nfiles = 0;

  for(int i = 2; i < argc; i++){
    if (argv[i][0] == '-'){
      if(strcmp(argv[i], "-L") == 0){
        break;
      }
      else{
        usage();
      }
    }
    else{
      nfiles++;
    }
  }

  int nfiles_temp = 0;
  for(int i = (3+nfiles); i < argc; i++){
    if (argv[i][0] == '-'){
      usage();
    }
    else{
      nfiles_temp++;
    }
  }

  if(nfiles_temp != nfiles){
      fprintf(stderr, "%s:%d\n"
              " Inconsistent number of files!"
              " %d data files and %d legofit files\n", __FILE__,__LINE__, nfiles, nfiles_temp);
      usage();
  }

  const char *legofname[nfiles];
  const char *datafname[nfiles];

  for(int i = 0; i < nfiles; ++i)
      datafname[i] = argv[i+1];
  for(int i = 0; i < nfiles; ++i)
      legofname[i] = argv[i+3+nfiles];

  // Read bootstrap files into an arrays of FIFO queues
  StrDblStack* data_stack[nfiles];
  StrDblStack* lego_stack[nfiles];

  StrDblStack* real_stack = parseSitPat(realfName);


  for(int i = 0; i < nfiles; ++i) {
      lego_stack[i] = parseSitPat(legofname[i]);
      data_stack[i] = parseSitPat(datafname[i]);

      if(StrDblStack_compare(real_stack, lego_stack[i])) {
          fprintf(stderr, "%s:%d: inconsistent parameters in"
                  " files%s and %s\n", __FILE__,__LINE__, realfName, legofname[i]);
          exit(EXIT_FAILURE);
      }
      if(StrDblStack_compare(real_stack, data_stack[i])) {
          fprintf(stderr, "%s:%d: inconsistent parameters in"
                  " files%s and %s\n", __FILE__,__LINE__, realfName, datafname[i]);
          exit(EXIT_FAILURE);
      }
  }

  //normalize the queues
  StrDblStack_normalize(real_stack);
  for(int i = 0; i < nfiles; ++i) {
      StrDblStack_normalize(lego_stack[i]);
      StrDblStack_normalize(data_stack[i]);
  }

  //find MSD
  double real_msd = 0;
  double boot_msd = 0;
  double bepe;
  double x;

  int npat;

  StrDblStack* temp_L;
  StrDblStack* temp_D;
  StrDblStack* temp_d;

  for (int i = 0; i < nfiles; ++i){
    temp_L = lego_stack[i];
    temp_D = data_stack[i];
    temp_d = real_stack;
    npat = StrDblStack_length(temp_D);
    for (int j = 0; j < npat; ++j){
      x = (temp_d->strdbl.val - temp_L->strdbl.val);
      real_msd += (x*x);

      x = (temp_D->strdbl.val - temp_L->strdbl.val);
      boot_msd += (x*x);

      temp_D = temp_D->next;
      temp_L = temp_L->next;
      temp_d = temp_d->next;

      fprintf(stderr, "Real: %lf\nBoot: %lf", real_msd, boot_msd);
    }
  }

  real_msd /= (nfiles*npat);
  boot_msd /= (nfiles*npat);

  bepe = real_msd + boot_msd;

  printf("BEPE = %lf\n", bepe);
}
