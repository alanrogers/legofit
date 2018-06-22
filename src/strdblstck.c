/**
 * @file strdblstck.c
 * @author Daniel R. Tabin and Alan R. Rogers
 * @brief Functions for Composite Likelihood Information Criterion.
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "hessian.h"
#include "misc.h"
#include "strdblstck.h"
#include <stdbool.h>

// Push a value onto the tail of the stack. Return pointer to new
// head. Example:
//
// StrDblStack *stack=NULL;
// stack = StrDblStack_push(stack, "name1", 1.0);
// stack = StrDblStack_push(stack, "name2", 2.0);
StrDblStack *StrDblStack_push(StrDblStack *self, const char *str, double val) {
    if(self != NULL) {
        self->next = StrDblStack_push(self->next, str, val);
        return self;
    }
    StrDblStack *new = malloc(sizeof(StrDblStack));
    CHECKMEM(new);
    new->strdbl.val = val;
    int status = snprintf(new->strdbl.str, sizeof(new->strdbl.str), "%s", str);
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
// StrDbl x;            NOTE: Make sure this is not NULL
// stack = StrDblStack_pop(stack, &x);  // x={"name1", 1.0}
// stack = StrDblStack_pop(stack, &x);  // x={"name2", 2.0}
StrDblStack *StrDblStack_pop(StrDblStack *self, StrDbl *strdbl) {
    if(self==NULL)
        return NULL;
    assert(strdbl != NULL);
    strdbl->val = self->strdbl.val;
    strcpy(strdbl->str, self->strdbl.str);
    StrDblStack *next = self->next;
    free(self);
    return next;
}

// //get a strdbl
// // NOTE: This is inefficient, we should work on making this faster
//
// StrDbl *StrDblStack_get(StrDblStack *self, StrDbl *strdbl, int index) {
//     if(self==NULL)
//         return NULL;
//
//     StrDblStack *temp;
//
//     for (int i = 0; i < index; i++)
//       temp = temp->next;
//
//     strdbl = &(temp->strdbl);
//     return strdbl;
// }

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

// Parse a legofit output file for CLIC. Return an object of type
// StrDblStack, which contains the number of parameters, their names,
// and their values.
StrDblStack *parseLegofit_CLIC(const char *fname) {
    FILE *fp = fopen(fname, "r");
    if(fp==NULL) {
        fprintf(stderr,"%s:%d: can't read file \"%s\"\n",
                __FILE__,__LINE__,fname);
        exit(EXIT_FAILURE);
    }
    char buff[2000];
    int got_fitted=0;
    StrDblStack *stack=NULL;
    while(1) {
        if(fgets(buff, sizeof buff, fp) == NULL) {
            break;
        }
        if(strchr(buff, '\n') == NULL && !feof(stdin)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        if(!got_fitted) {
            if(strncmp("Fitted", buff, 6) == 0)
                got_fitted=1;
            continue;
        }else if(got_fitted) {
            if(NULL != strstr(buff, "constrained")) {
                got_fitted = 0;
                continue;
            }
        }
        char *valstr = buff;
        char *name = strsep(&valstr, "=");
        if(name==NULL || valstr==NULL)
            continue;
        name = stripWhiteSpace(name);
        valstr = stripWhiteSpace(valstr);
        stack=StrDblStack_push(stack, name, strtod(valstr, NULL) );
    }
    assert(StrDblStack_length(stack) > 0);
    return stack;
}

// Parse a data file for BEPE. Return an object of type
// StrDblStack, which contains the number of parameters, their names,
// and their values.

StrDblStack *parseSitPat(const char *fname) {
    FILE *fp = fopen(fname, "r");
    if(fp==NULL) {
        fprintf(stderr,"%s:%d: can't read file \"%s\"\n",
                __FILE__,__LINE__,fname);
        exit(EXIT_FAILURE);
    }
    char buff[2000];
    bool got_sitepat = false;
    StrDblStack *stack=NULL;
    while(1) {
        if(fgets(buff, sizeof buff, fp) == NULL) {
            break;
        }
        if(strchr(buff, '\n') == NULL && !feof(stdin)) {
            fprintf(stderr, "%s:%d: Buffer overflow. size=%zu\n",
                    __FILE__, __LINE__, sizeof(buff));
            exit(EXIT_FAILURE);
        }
        if(!got_sitepat) {
            char* no_spaces_buff = stripInternalWhiteSpace(buff);
            if(strncmp("#SitePat", no_spaces_buff, 8) == 0)
                got_sitepat=true;
            continue;
        }

        char* temp = buff;
        char* name = strtok_r(temp, " ", &temp);
        char* valstr = strtok_r(temp, " ", &temp);

        if(name==NULL || valstr==NULL)
            continue;
        name = stripWhiteSpace(name);
        valstr = stripWhiteSpace(valstr);
        stack=StrDblStack_push(stack, name, strtod(valstr, NULL) );
    }
    return stack;
}

void StrDblStack_normalize(StrDblStack* self){
  int length = StrDblStack_length(self);
  double total = 0;
  StrDblStack* temp;

  for(temp = self; temp; temp = temp->next){
    total += temp->strdbl.val;
  }

  temp = self;
  for(int i = 0; i < length; i++){
    temp->strdbl.val = ((temp->strdbl.val)/total);
    temp = temp->next;
  }
}

// On input, nfiles and npar are the number of rows and columns in
// the input data matrix "array". The npar X npar matrix "cov" should
// be allocated in the calling function.
//
// On return, cov[i][j] is the covariance between the i'th and j'th
// columns of "array".

void make_covar_matrix(int nfiles, int npar, double array[nfiles][npar],
                      gsl_matrix *cov){
   assert(npar == cov->size1);
   assert(npar == cov->size2);
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
           gsl_matrix_set(cov, i, j, covar_sum/nfiles);
       }
   }
}
