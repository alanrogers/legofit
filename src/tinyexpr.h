/**
 * This version of tinyexpr has been modified from Lewis Van Winkle's original. 
 * The revised version uses hash tables to map names to functions and to
 * variables.
 * 
 * Alan Rogers
 *
 * TINYEXPR - Tiny recursive descent parser and evaluation engine in C
 *
 * Copyright (c) 2015, 2016 Lewis Van Winkle
 *
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 * claim that you wrote the original software. If you use this software
 * in a product, an acknowledgement in the product documentation would be
 * appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 * misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#ifndef __TINYEXPR_H__
#define __TINYEXPR_H__

#define TE_NAT_LOG

#include "strptrmap.h"
#include <stdio.h>

typedef struct te_expr {
    int type;
    union {
        double value;
        const double *bound;
        void *function;
    };
    void *parameters[1];
} te_expr;

enum {
    TE_VARIABLE = 0,

    TE_FUNCTION0 = 8, TE_FUNCTION1, TE_FUNCTION2, TE_FUNCTION3,
    TE_FUNCTION4, TE_FUNCTION5, TE_FUNCTION6, TE_FUNCTION7,

    TE_CLOSURE0 = 16, TE_CLOSURE1, TE_CLOSURE2, TE_CLOSURE3,
    TE_CLOSURE4, TE_CLOSURE5, TE_CLOSURE6, TE_CLOSURE7,

    TE_FLAG_PURE = 32
};

typedef struct te_variable {
    void *address;
    int type;
    void *context;
} te_variable;

te_variable *te_variable_new(void *address);
void te_free_variables(StrPtrMap *pars);

/** 
 * Parses the input expression, evaluates it, and frees it.
 * Returns NaN on error. */
double te_interp(const char *expression, int *error);

/** Parses the input expression and binds variables.
 * Returns NULL on error. 
 */
te_expr *te_compile(const char *expression, StrPtrMap * varmap,
                    int *error);

/** Evaluates the expression. */
double te_eval(const te_expr * n);

/** Prints debugging information on the syntax tree. */
void te_print(const te_expr * n, FILE * fp);

/**
 * Fill array "ptr" with pointers to the variables on which this
 * expression depends. Return the number of dependent variables, which
 * should be less than or equal to len. Abort if len is smaller than
 * the number of dependent variables.
 */ 
int te_dependencies(const te_expr *self, int len, const double *ptr[len]);

/** Frees the expression. */
/* This is safe to call on NULL pointers. */
void te_free(te_expr * n);

/** Free hash table used for function names. */
void te_free_func_map(void);

#endif /*__TINYEXPR_H__*/
