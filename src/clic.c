/**
 * @file clic.c
 * @author Daniel R. Tabin and Alan R. Rogers
 * @brief Functions for Composite Likelihood Information Criterion.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "hessian.h"
#include "strdblstck.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void matrix_mult(double** m1, int m1_rows, int m1_cols, double** m2, int m2_rows, int m2_cols, double** results);

void usage(void);

double KL_to_lnL(double KL, double* p_matrix, int p_matrix_size, double sum);

/*
 Matrix multiplier.  Takes two matricies of doubles m1 and m2 and their sizes,
 as well as a results matrix, where the results go.
*/
 void matrix_mult(double** m1, int m1_rows, int m1_cols,
                  double** m2, int m2_rows, int m2_cols,
                  double** results){
   int i, j, k;

   for (i = 0; i < m1_rows; i++){         //go through each row
     for(j = 0; j < m2_cols; j++){        //go through each coloumn
       int x = 0;
       for (k = 0; k < m2_rows; k++){
         x += (m1[i][k] * m2[k][j]);      //calculate what goes there
       }
       results[i][j] = x;                 //put it in results
     }
   }
 }

#if 0
/*
  Parse legofit output.
  Creates a double array with the first dimension being number of files,
  and the second being number number of paramaters.  This is used as a
  file parser
*/
 double** get_fit_param_array(int num_files, char *file_name[num_files]){
   char file_base[100];

   FILE* f;

   double** array = (double**) malloc(num_files * sizeof(double*));

   for (int i = 0; i < num_files; i++){             //go through each file
     f = fopen(file_name[i], "r");
     if(f == NULL){
         fprintf(stderr, "Error, invalid file name: %s\n", file_name);
         exit(EXIT_FAILURE);
     }

     char input[100];
     int param_num = 0;

     array[i] = (double*) malloc(num_params * sizeof(double));

     do {                                       //fscanf until past DiffEv
       fscanf(f, "%s", input);
     } while(strcmp(input, "DiffEv") != 0);

     while(fscanf(f, "%s", input)){             //find and place data into array
       if(strcmp(input, "=") == 0){
         fscanf(f, "%lf", &array[i][param_num]);
         param_num++;
         if(param_num >= num_params){
           break;
         }
       }
     }
   }

   return array;
 }

/*
  Takes the start of the string and finds how many of that file you have, then
  sends that information to get_fit_param_array and returns its result
*/
 double** get_fit_param_array_num_unkown(char* title, int num_params){
   int num_files = 0;

   char file_name[200];
   FILE* f;

   snprintf(file_name, sizeof(file_name), "%sboot%d.legofit",title, num_files);

   while ((f = fopen(file_name, "r"))){
      snprintf(file_name, sizeof(file_name), "%sboot%d.legofit",title, num_files);
     num_files++;
   }

   return get_fit_param_array(title, num_files, num_params);
 }
#endif

/*
  Takes KL and converts it to the natural log of likelihood
*/

 double KL_to_lnL(double KL, double* p_matrix, int p_matrix_size, double sum){
   double lnL;
   double entropy = 0.0;

   for (int i = 0; i < p_matrix_size; i++){
     entropy -= (p_matrix[i] * log(p_matrix[i]));
   }

   lnL = -sum*(entropy + KL);

   return lnL;
 }

const char *usageMsg =
     "usage: clic <file.pts> <boot1.legofit> <boot2.legofit> ...\n"
     "  where file.pts is the .pts file produced by legofit with the real\n"
     "  data, and each \"boot\" file is the legofit output from one bootstrap\n"
     "  replicate. Must include .pts file and at least 2 boostrap files.\n";

void usage(void) {
    fputs(usageMsg, stderr);
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv){
    int i, j;

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
        stack[i] = parseLegofit_CLIC(bootfname[i]);
        if(i>0) {
            if(StrDblStack_compare(stack[0], stack[i])) {
                fprintf(stderr, "%s:%d: inconsistent parameters in"
                        " files %s and %s\n", __FILE__,__LINE__,
                        bootfname[0], bootfname[i]);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Use the stacks to populate an array of parameter names
    // and a matrix of parameter values. Rows are bootstrap
    // replicates; columns are parameters.
    int npar = StrDblStack_length(stack[0]);
    char *parname[npar];
    double datmat[nfiles][npar];
    for(i=0; i < nfiles; ++i) {
        for(j=0; j < npar; ++j) {
            StrDbl strdbl;
            stack[i] = StrDblStack_pop(stack[i], &strdbl);
            datmat[i][j] = strdbl.val;
            if(i==0)
                parname[j] = strdup(strdbl.str);
        }
        assert(stack[i] == NULL); // check that stacks are freed
    }

    // Print data matrix with column header
    for(j=0; j < npar; ++j)
        printf(" %s", parname[j]);
    putchar('\n');
    for(i=0; i<nfiles; ++i) {
        for(j=0; j < npar; ++j)
            printf(" %lg", datmat[i][j]);
        putchar('\n');
    }

    // Make covariance matrix
    gsl_matrix *c_matrix = gsl_matrix_alloc(npar,npar);
    make_covar_matrix(nfiles, npar, datmat, c_matrix);

    // Print it
    for (j = 0; j < npar; j++)
        printf(" %8s", parname[j]);
    putchar('\n');
    for (i = 0; i < npar; i++){
        for (j = 0; j < npar; j++){
            printf(" %8.2lg", gsl_matrix_get(c_matrix, i, j));
        }
        printf("\n");
    }

    Hessian hesobj = hessian(ptsfname);
    gsl_matrix *H = hesobj.hessian;
    char **Hparname = hesobj.parname;
    double lnL = hesobj.lnL;

    // Do the two sets of parameters match? (One from bootstrap files,
    // the other from pts file.)
    int mismatch = 0;
    if(hesobj.npar != npar)
        mismatch=1;
    for(i=0; mismatch==0 && i<npar; ++i) {
        if(0 != strcmp(parname[i], Hparname[i]))
            mismatch=1;
    }
    if(mismatch) {
        fprintf(stderr,"%s:%d: mismatch between parameters"
                " in files %s and %s\n", __FILE__,__LINE__,
                ptsfname, bootfname[0]);
        exit(EXIT_FAILURE);
    }

    // allocate matrix to hold product H*c_matrix
    gsl_matrix *HC = gsl_matrix_alloc(npar,npar);

    // form matrix product HC = H*c_matrix
    gsl_blas_dgemm(CblasNoTrans, // don't transpose H
                   CblasNoTrans, // don't transpose C
                   1.0,          // scalar
                   H,            // Hermitian matrix
                   c_matrix,     // covariance matrix
                   0.0,          // scalar
                   HC            // product
                   );

    double trace=0.0;
    for(i=0; i<npar; ++i)
        trace += gsl_matrix_get(HC, i, i);

    /*
      This is the information criterion of Varin, Cristiano and
      Vidoni, Paolo. 2005. A note on composite likelihood inference
      and model selection. Biometrika 92(3):519-528. Eqn 5, p 523.

      Note that "trace" should be negative at a local maximum, so clic
      will be smaller than lnL.
     */
    double clic = lnL + trace;

    printf("CLIC: %0.8lg\n", clic);

    gsl_matrix_free(H);
    gsl_matrix_free(HC);
    for(i=0; i<npar; ++i)
        free(Hparname[i]);
    free(Hparname);
    return 0;
 }
