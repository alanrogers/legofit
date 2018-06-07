#define NDEBUG
/**
 * @file clic.c
 * @author Daniel R. Tabin
 * @brief Functions for Composite Likelihood Information Criterion.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <math.h>
 #include "clic.h"


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

/*
  Parse legofit output.
  Creates a double array with the first dimension being number of files,
  and the second being number number of paramaters.  This is used as a
  file parser
*/
 double** get_fit_param_array(char* title, int num_files, int num_params){
   char file_base[100];
   char file_name[200];

   snprintf(file_base, size_of(file_base), "%sboot",title);

   FILE* f;

   double** array = (double**) malloc(num_files * sizeof(double*));

   for (int i = 0; i < num_files; i++){             //go through each file
     snprintf(file_name, size_of(file_name), "%sboot%d.legofit", file_base, i);

     f = fopen(file_name, "r");
     if(f == NULL){
       printf("Error, invalid file name\n");
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

   snprintf(file_name, size_of(file_name), "%sboot%d.legofit",title, num_files);

   while ((f = fopen(file_name, "r"))){
      snprintf(file_name, size_of(file_name), "%sboot%d.legofit",title, num_files);
     num_files++;
   }

   return get_fit_param_array(title, num_files, num_params);
 }


 /*
    Takes an array of paramaters and returns a matrix of covariances
    NOTE: This flips the array and the matrix goes from having the first
    index be the file to the second (and third) index being the file.
 */
 double** make_covar_matrix(double** array, int files, int params){
   double param_averages[params];
   double covar_sum;

   for(int i = 0; i < params; i++){
     param_averages[i] = 0;
     for (int j = 0; j < files; j++){
       param_averages[i] += array[j][i];
     }
     param_averages[i] /= files;
   }

   double** covar_matrix = (double**) malloc(params * sizeof(double*));

   for(int i = 0; i < params; i++){
     covar_matrix[i] = (double*) malloc(params * sizeof(double));
     for (int j = 0; j < params; j++){
       covar_sum = 0;
       for (int k = 0; k < files; k++){
         covar_sum += (array[k][j] - param_averages[j])
           * (array[k][i] - param_averages[i]);
       }
       covar_matrix[i][j] = (covar_sum/files);
     }
   }

   return covar_matrix;
 }

/*
  Takes KL and converts it to the natural log of likelihood
*/

 double KL_to_lnL(double KL, double* p_matrix, int p_matrix_size, double sum){
   double likelihood;
   double entropy = 0;

   for (int i = 0; i < p_matrix_size; i++){
     entropy += (p_matrix[i] * log(p_matrix[i]));
   }

   likelihood = entropy - KL;
   likelihood = likelihood * sum;

   return likelihood;
 }

#ifdef TEST

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

 int main(){
   double** array;
   double** c_matrix;
   int a = 1;
   int b = 9;

   array = get_fit_param_array("s1", a, b);

   for (int i = 0; i < a; i++){
     for (int j = 0; j < b; j++){
       printf("%lf\n", array[i][j]);
     }
   }

   c_matrix = make_covar_matrix(array, 1, 9);

   for (int i = 0; i < b; i++){
     for (int j = 0; j < a; j++){
       printf("%f ", c_matrix[i][j]);
     }
     printf("\n");
   }

   printf("done\n");
 }

#endif
