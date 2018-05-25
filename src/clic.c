/**
 * @file clic.c
 * @author Daniel R. Tabin
 * @brief Functions for Composite Likelihood Information Criterion.
 * @copyright Copyright (c) 2016, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

 #include "clic.h"


/*
 Matrix multiplier.  Takes two matricies of floats m1 and m2 and their sizes,
 as well as a results matrix, where the results go.
*/
 void mm_mult(float **m1, int m1_rows, int m1_cols,
              float **m2, int m2_rows, int m2_cols,
              float **results){
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
