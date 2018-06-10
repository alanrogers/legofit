/**
 * @file clic.c
 * @author Daniel R. Tabin and Alan R. Rogers
 * @brief Functions for Composite Likelihood Information Criterion.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

 #include "clic.h"
 #include "hessian.h"
 #include "misc.h"
 #include <stdio.h>
 #include <stdlib.h>
 #include <assert.h>
 #include <string.h>
 #include <math.h>

typedef struct StrDblStack StrDblStack;
typedef struct StrDbl StrDbl;

struct StrDbl {
    char str[100];
    double val;
};

// A FIFO stack. New values are pushed onto the tail. Old ones are
// popped off of the head.
struct StrDblStack {
    struct StrDblStack *next;
    struct StrDbl strdbl;
};

void usage(void);
StrDblStack *StrDblStack_push(StrDblStack *prev, char *str, double val);
StrDblStack *StrDblStack_pop(StrDblStack *self, StrDbl *strdbl);
StrDblStack *StrDblStack_free(StrDblStack *self);
int          StrDblStack_length(StrDblStack *self);
void         StrDblStack_print(StrDblStack *self, FILE *fp);
int          StrDblStack_compare(StrDblStack *lhs, StrDblStack *rhs);
StrDblStack *parseLegofit(const char *fname);
void make_covar_matrix(int nfiles, int npar, double array[nfiles][npar],
                       double cov[npar][npar]);

// Push a value onto the tail of the stack. Return pointer to new
// head. Example:
//
// StrDblStack *stack=NULL;
// stack = StrDblStack_push(stack, "name1", 1.0);
// stack = StrDblStack_push(stack, "name2", 2.0);
StrDblStack *StrDblStack_push(StrDblStack *self, char *str, double val) {
    if(self != NULL) {
        self->next = StrDblStack_push(self->next, str, val);
        return self;
    }
    StrDblStack *new = malloc(sizeof(StrDblStack));
    CHECKMEM(new);
    new->strdbl.val = val;
    int status = snprintf(new->strdbl.str, sizeof(new->strdbl.str),
                          "%s", str);
    if(status > sizeof(new->strdbl.str)) {
        fprintf(stderr, "%s:%d: buffer overflow\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    new->next = NULL;
    return new;
}

// Pop a value off the head of the stack. Return pointer to new
// head. Example:
//
// StrDblStack *stack=NULL;
// stack = StrDblStack_push(stack, "name1", 1.0);
// stack = StrDblStack_push(stack, "name2", 2.0);
//
// StrDbl x;
// stack = StrDblStack_pop(stack, &x);  // x={"name1", 1.0}
// stack = StrDblStack_pop(stack, &x);  // x={"name2", 2.0}
StrDblStack *StrDblStack_pop(StrDblStack *self, StrDbl *strdbl) {
    if(self==NULL)
        return NULL;
    strdbl->val = self->strdbl.val;
    strcpy(strdbl->str, self->strdbl.str);
    StrDblStack *next = self->next;
    free(self);
    return next;
}

int StrDblStack_length(StrDblStack *self) {
    if(self==NULL)
        return 0;
    return 1 + StrDblStack_length(self->next);
}

StrDblStack *StrDblStack_free(StrDblStack *self) {
    if(self) {
        self->next = StrDblStack_free(self->next);
        free(self);
    }
    return NULL;
}

void StrDblStack_print(StrDblStack *self, FILE *fp) {
    while(self) {
        fprintf(fp,"%s = %lg\n", self->strdbl.str, self->strdbl.val);
        self = self->next;
    }
}

/**
 * Compare the str fields in two StrDblStack objects. Return -1, 0, or 1
 * if the lhs is less than, equal to, or greater than rhs.
 */
int StrDblStack_compare(StrDblStack *lhs, StrDblStack *rhs) {
    while(lhs && rhs) {
        int cmp = strcmp(lhs->strdbl.str, rhs->strdbl.str);
        if(cmp)
            return cmp;
        lhs = lhs->next;
        rhs = rhs->next;
    }
    if(lhs)
        return 1;
    if(rhs)
        return -1;
    return 0;
}

// Parse a legofit output file. Return an object of type StrDblStack,
// which contains the number of parameters, their names, and their values.
StrDblStack *parseLegofit(const char *fname) {
    FILE *fp = fopen(fname, "r");
    if(fp==NULL) {
        fprintf(stderr,"%s:%d: can's read file \"%s\"\n",
                __FILE__,__LINE__,fname);
        exit(EXIT_FAILURE);
    }
    char buff[500];
    int got_fitted=0;
    StrDblStack *stack=NULL;
    while(1) {
        if(NULL == fgets(buff, sizeof buff, fp)) {
            break;
        }
        if(NULL == strchr(buff, '\n') && !feof(stdin)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        if(!got_fitted) {
            if(0 == strncmp("Fitted", buff, 6))
                got_fitted=1;
            continue;
        }
        char *valstr = buff;
        char *name = strsep(&valstr, "=");
        if(name==NULL || valstr==NULL)
            continue;
        name = stripWhiteSpace(name);
        valstr = stripWhiteSpace(valstr);
        stack=StrDblStack_push(stack, name, strtod(valstr, NULL) );
    }
    return stack;
}

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
   On input, nfiles and npar are the number of rows and columns in
   the input data matrix "array". The npar X npar matrix "cov" should
   be allocated in the calling function.

   On return, cov[i][j] is the covariance between the i'th and j'th
   columns of "array".
 */
void make_covar_matrix(int nfiles, int npar, double array[nfiles][npar],
                       double cov[npar][npar]){
    double param_averages[npar];
    double covar_sum;
    int i, j, k;

    for(i = 0; i < npar; i++){
        param_averages[i] = 0;
        for (j = 0; j < nfiles; j++){
            param_averages[i] += array[j][i];
        }
        param_averages[i] /= nfiles;
    }

    for(i = 0; i < npar; i++){
        for (j = 0; j < npar; j++){
            covar_sum = 0;
            for (k = 0; k < nfiles; k++){
                covar_sum += (array[k][j] - param_averages[j])
                    * (array[k][i] - param_averages[i]);
            }
            cov[i][j] = (covar_sum/nfiles);
        }
    }
}

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
        stack[i] = parseLegofit(bootfname[i]);
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
    double c_matrix[npar][npar];
    make_covar_matrix(nfiles, npar, datmat, c_matrix);

    // Print it
    for (j = 0; j < npar; j++)
        printf(" %8s", parname[j]);
    putchar('\n');
    for (i = 0; i < npar; i++){
        for (j = 0; j < npar; j++){
            printf(" %8.2lg", c_matrix[i][j]);
        }
        printf("\n");
    }

    Hessian hesobj = hessian(ptsfname);
    return 0;
 }

