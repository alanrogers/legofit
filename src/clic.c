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
 #include "clic.h"


/*
 Matrix multiplier.  Takes two matricies of floats m1 and m2 and their sizes,
 as well as a results matrix, where the results go.
*/
 void matrix_mult(float** m1, int m1_rows, int m1_cols,
                  float** m2, int m2_rows, int m2_cols,
                  float** results){
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
  Puts an int (data) into a string (buffer)
*/
 char* itoa(int data, char* buffer){
   int temp = data;      //copy data
   int size = 1;         //base string size
   while ((int)(temp/10) != 0){      //while there is still stuff in temp
     temp /= 10;                     //divide temp by 10
     size++;                         //add one to size
   }
   for(int i = size; i > 0; i--){
     *(buffer+size-1) = (char)((data%10)+'0');   //place data in array
     data /= 10;
   }
   return buffer;        //return
 }

  //very basic method, turns tree strings into one
 char* tri_cat(char* a, char* b, char* c, char* str){
   strcpy(str,a);
   strcat(str,b);
   strcat(str,c);
   return str;
 }

 float** get_fit_param_array(char* title, int num_files){

 }

/*
  Takes the start of the string and finds how many of that file you have, then
  sends that information to get_fit_param_array and returns its result
*/
 float** get_fit_param_array_num_unkown(char* title){
   int num_files = 0;
   char* file_base = title;
   strcat(file_base, "boot");
   char file_num[5];
   char file_name[20];
   FILE* f;

   while (f = fopen(tri_cat(file_base, itoa(num_files, file_num), ".legofit", file_name), "r")){
     num_files++;
   }

   return get_fit_param_array(title, num_files);
 }

 int main(){
   printf("test\n");
 }
