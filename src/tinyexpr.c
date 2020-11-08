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

/* COMPILE TIME OPTIONS */

/* Exponentiation associativity:
For a^b^c = (a^b)^c and -a^b = (-a)^b do nothing.
For a^b^c = a^(b^c) and -a^b = -(a^b) uncomment the next line.*/

/* #define TE_POW_FROM_RIGHT */

/* Logarithms
For log = base 10 log do nothing
For log = natural log uncomment the next line. */

#define TE_NAT_LOG

#include "tinyexpr.h"
#include "misc.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>
#include <ctype.h>
#include <pthread.h>

#ifndef NAN
#define NAN (0.0/0.0)
#endif

#ifndef INFINITY
#define INFINITY (1.0/0.0)
#endif

typedef double (*te_fun2) (double, double);

enum {
    TOK_NULL = TE_CLOSURE7 + 1, TOK_ERROR, TOK_END, TOK_SEP,
    TOK_OPEN, TOK_CLOSE, TOK_NUMBER, TOK_VARIABLE, TOK_INFIX
};

enum { TE_CONSTANT = 1 };

typedef struct state {
    const char *start;
    const char *next;
    int         type;
    union {
        double      value;
        const double *bound;
        void *function;
    };
    void       *context;

    StrPtrMap *lookup;
} state;

#define TYPE_MASK(TYPE) ((TYPE)&0x0000001F)

#define IS_PURE(TYPE) (((TYPE) & TE_FLAG_PURE) != 0)
#define IS_FUNCTION(TYPE) (((TYPE) & TE_FUNCTION0) != 0)
#define IS_CLOSURE(TYPE) (((TYPE) & TE_CLOSURE0) != 0)
#define ARITY(TYPE) ( ((TYPE) & (TE_FUNCTION0 | TE_CLOSURE0)) ? ((TYPE) & 0x00000007) : 0 )
#define NEW_EXPR(type, ...) new_expr((type), (const te_expr*[]){__VA_ARGS__})

void te_free_parameters(te_expr * n);
void next_token(state * s);
static double pi(void);
static double e(void);
static te_variable *new_func(void *address, int type, void *context);
static void init_func_map(void);
int string_copy(int dlen, char dest[dlen], int slen, const char src[slen]);

te_variable *te_variable_new(void *address) {
    te_variable *self = malloc(sizeof(te_variable));
    if(self==NULL)
        return NULL;
    memset(self, 0, sizeof(te_variable));
    self->address = address;
    self->type = TE_VARIABLE;
    return self;
}

// Free all variables and the hash map that holds them.
void te_free_variables(StrPtrMap *pars) {

    unsigned len = StrPtrMap_size(pars);
    if(len > 0) {
        void *v[len];
        StrPtrMap_ptrArray(pars, len, v);
        for(unsigned i=0; i<len; ++i) {
            te_variable *var = v[i];
            free(var);
        }
    }
    StrPtrMap_free(pars);
}

static te_expr *new_expr(const int type, const te_expr * parameters[]) {
    const int   arity = ARITY(type);
    const int   psize = sizeof(void *) * arity;
    const int   size =
        (sizeof(te_expr) - sizeof(void *)) + psize +
        (IS_CLOSURE(type) ? sizeof(void *) : 0);
    te_expr    *ret = malloc(size);
    memset(ret, 0, size);
    if(arity && parameters) {
        memcpy(ret->parameters, parameters, psize);
    }
    ret->type = type;
    ret->bound = 0;
    return ret;
}

void te_free_parameters(te_expr * n) {
    if(!n)
        return;
    switch (TYPE_MASK(n->type)) {
    case TE_FUNCTION7:
    case TE_CLOSURE7:
        te_free(n->parameters[6]);
    case TE_FUNCTION6:
    case TE_CLOSURE6:
        te_free(n->parameters[5]);
    case TE_FUNCTION5:
    case TE_CLOSURE5:
        te_free(n->parameters[4]);
    case TE_FUNCTION4:
    case TE_CLOSURE4:
        te_free(n->parameters[3]);
    case TE_FUNCTION3:
    case TE_CLOSURE3:
        te_free(n->parameters[2]);
    case TE_FUNCTION2:
    case TE_CLOSURE2:
        te_free(n->parameters[1]);
    case TE_FUNCTION1:
    case TE_CLOSURE1:
        te_free(n->parameters[0]);
    }
}

void te_free(te_expr * n) {
    if(!n)
        return;
    te_free_parameters(n);
    free(n);
}

static double pi(void) {
    return 3.14159265358979323846;
}
static double e(void) {
    return 2.71828182845904523536;
}
static double fac(double a) {   /* simplest version of fac */
    if(a < 0.0)
        return NAN;
    if(a > UINT_MAX)
        return INFINITY;
    unsigned int ua = (unsigned int) (a);
    unsigned long int result = 1, i;
    for(i = 1; i <= ua; i++) {
        if(i > ULONG_MAX / result)
            return INFINITY;
        result *= i;
    }
    return (double) result;
}
static double ncr(double n, double r) {
    if(n < 0.0 || r < 0.0 || n < r)
        return NAN;
    if(n > UINT_MAX || r > UINT_MAX)
        return INFINITY;
    unsigned long int un = (unsigned int) n, ur = (unsigned int) r, i;
    unsigned long int result = 1;
    if(ur > un / 2)
        ur = un - ur;
    for(i = 1; i <= ur; i++) {
        if(result > ULONG_MAX / (un - ur + i))
            return INFINITY;
        result *= un - ur + i;
        result /= i;
    }
    return result;
}
static double npr(double n, double r) {
    return ncr(n, r) * fac(r);
}

static pthread_mutex_t func_map_lock = PTHREAD_MUTEX_INITIALIZER;
static StrPtrMap *func_map=NULL;

/// Free func_map
void te_free_func_map(void) {
    int status = pthread_mutex_lock(&func_map_lock);
    if(status)
        ERR(status, "lock");

    unsigned len = StrPtrMap_size(func_map);
    if(len > 0) {
        void *v[len];
        StrPtrMap_ptrArray(func_map, len, v);
        for(unsigned i=0; i<len; ++i) {
            te_variable *var = v[i];
            free(var);
        }
    }
    StrPtrMap_free(func_map);

    status = pthread_mutex_unlock(&func_map_lock);
    if(status)
        ERR(status, "unlock");

    func_map = NULL;
}

/// Initialize func_map; abort on failure.
static void init_func_map(void) {
    int status = pthread_mutex_lock(&func_map_lock);
    if(status)
        ERR(status, "lock");

    func_map =  StrPtrMap_new();
    if(func_map == NULL) {
        fprintf(stderr,"%s:%d: can't allocate func_map\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    status = 0;

    status |= StrPtrMap_insert(func_map, "abs",
                               new_func(fabs, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "acos",
                               new_func(acos, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "asin",
                               new_func(asin, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "atan",
                               new_func(atan, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "atan2",
                               new_func(atan2, TE_FUNCTION2 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "ceil",
                               new_func(ceil, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "cos",
                               new_func(cos, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "cosh",
                               new_func(cosh, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "e",
                               new_func(e, TE_FUNCTION0 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "exp",
                               new_func(exp, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "fac",
                               new_func(fac, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "floor",
                               new_func(floor, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "ln",
                               new_func(log, TE_FUNCTION1 | TE_FLAG_PURE, 0));
#ifdef TE_NAT_LOG
    status |= StrPtrMap_insert(func_map, "log",
                               new_func(log, TE_FUNCTION1 | TE_FLAG_PURE, 0));
#else
    status |= StrPtrMap_insert(func_map, "log",
                               new_func(log10, TE_FUNCTION1 | TE_FLAG_PURE, 0));
#endif
    status |= StrPtrMap_insert(func_map, "log10",
                               new_func(log10, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "ncr",
                               new_func(ncr, TE_FUNCTION2 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "npr",
                               new_func(npr, TE_FUNCTION2 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "pi",
                               new_func(pi, TE_FUNCTION0 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "pow",
                               new_func(pow, TE_FUNCTION2 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "sin",
                               new_func(sin, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "sinh",
                               new_func(sinh, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "sqrt",
                               new_func(sqrt, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "tan",
                               new_func(tan, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    status |= StrPtrMap_insert(func_map, "tanh",
                               new_func(tanh, TE_FUNCTION1 | TE_FLAG_PURE, 0));
    if(status) {
        fprintf(stderr,"%s:%s:%d: can't initialize func_map\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }

    status = pthread_mutex_unlock(&func_map_lock);
    if(status)
        ERR(status, "unlock");
}

static te_variable *new_func(void *address, int type,
                                 void *context) {
    te_variable *new = malloc(sizeof(te_variable));
    if(new == NULL)
        return NULL;
    new->address = address;
    new->type = type;
    new->context = context;
    return new;
}

static const te_variable *find_builtin(const char *name) {
    if(func_map == NULL)
        init_func_map();

    int status = pthread_mutex_lock(&func_map_lock);
    if(status)
        ERR(status, "lock");

    te_variable *rval = StrPtrMap_get(func_map, name);

    status = pthread_mutex_unlock(&func_map_lock);
    if(status)
        ERR(status, "unlock");

    return rval;
}

static double add(double a, double b) {
    return a + b;
}
static double sub(double a, double b) {
    return a - b;
}
static double mul(double a, double b) {
    return a * b;
}
static double divide(double a, double b) {
    return a / b;
}
static double negate(double a) {
    return -a;
}
static double comma(double a, double b) {
    (void) a;
    return b;
}

/// Copy slen bytes from src into dest, and add a terminating '\0' character.
/// Return 0 on success, 1 on bufferflow.
int string_copy(int dlen, char dest[dlen], int slen, const char src[slen]) {
    if(dlen < slen + 1)
        return 1;
    int i;
    for(i=0; i < slen; ++i)
        dest[i] = src[i];
    dest[i] = '\0';
    return 0;
}

void next_token(state * s) {
    s->type = TOK_NULL;

    char name[100];

    do {

        if(!*s->next) {
            s->type = TOK_END;
            return;
        }

        /* Try reading a number. */
        char *end;
        errno = 0;
        s->value = strtod(s->next, &end);
        if(errno) {
            /* overflow or underflow */
            s->type = TOK_ERROR;
            s->next = end;
        }else if(end != s->next) {
            s->next = end;
            s->type = TOK_NUMBER;
        } else {
            /* Look for a variable or builtin function call. */
            if(isalpha(*s->next)) {
                const char *start;
                start = s->next;
                while(isalnum(*s->next) || *s->next == '_')
                    s->next++;
                int status = string_copy(sizeof(name), name,
                                         s->next-start, start);
                if(status) {
                    fprintf(stderr,"%s:%d: buffer overflow\n",
                            __FILE__,__LINE__);
                    exit(EXIT_FAILURE);
                }
                const te_variable *var = StrPtrMap_get(s->lookup, name);
                if(!var)
                    var = find_builtin(name);

                if(!var) {
                    s->type = TOK_ERROR;
                } else {
                    switch (TYPE_MASK(var->type)) {
                    case TE_VARIABLE:
                        s->type = TOK_VARIABLE;
                        s->bound = var->address;
                        break;

                    case TE_CLOSURE0:
                    case TE_CLOSURE1:
                    case TE_CLOSURE2:
                    case TE_CLOSURE3:
                    case TE_CLOSURE4:
                    case TE_CLOSURE5:
                    case TE_CLOSURE6:
                    case TE_CLOSURE7:
                        s->context = var->context;

                    case TE_FUNCTION0:
                    case TE_FUNCTION1:
                    case TE_FUNCTION2:
                    case TE_FUNCTION3:
                    case TE_FUNCTION4:
                    case TE_FUNCTION5:
                    case TE_FUNCTION6:
                    case TE_FUNCTION7:
                        s->type = var->type;
                        s->function = var->address;
                        break;
                    }
                }

            } else {
                /* Look for an operator or special character. */
                switch (s->next++[0]) {
                case '+':
                    s->type = TOK_INFIX;
                    s->function = add;
                    break;
                case '-':
                    s->type = TOK_INFIX;
                    s->function = sub;
                    break;
                case '*':
                    s->type = TOK_INFIX;
                    s->function = mul;
                    break;
                case '/':
                    s->type = TOK_INFIX;
                    s->function = divide;
                    break;
                case '^':
                    s->type = TOK_INFIX;
                    s->function = pow;
                    break;
                case '%':
                    s->type = TOK_INFIX;
                    s->function = fmod;
                    break;
                case '(':
                    s->type = TOK_OPEN;
                    break;
                case ')':
                    s->type = TOK_CLOSE;
                    break;
                case ',':
                    s->type = TOK_SEP;
                    break;
                case ' ':
                case '\t':
                case '\n':
                case '\r':
                    break;
                default:
                    s->type = TOK_ERROR;
                    break;
                }
            }
        }
    } while(s->type == TOK_NULL);
}

static te_expr *list(state * s);
static te_expr *expr(state * s);
static te_expr *power(state * s);

static te_expr *base(state * s) {
    /* <base>  =   <constant> | <variable> | <function-0> {"(" ")"} |
       <function-1> <power> | <function-X> "(" <expr> {"," <expr>} ")" |
       "(" <list> ")" */
    te_expr    *ret;
    int         arity;

    switch (TYPE_MASK(s->type)) {
    case TOK_NUMBER:
        ret = new_expr(TE_CONSTANT, 0);
        ret->value = s->value;
        next_token(s);
        break;

    case TOK_VARIABLE:
        ret = new_expr(TE_VARIABLE, 0);
        ret->bound = s->bound;
        next_token(s);
        break;

    case TE_FUNCTION0:
    case TE_CLOSURE0:
        ret = new_expr(s->type, 0);
        ret->function = s->function;
        if(IS_CLOSURE(s->type))
            ret->parameters[0] = s->context;
        next_token(s);
        if(s->type == TOK_OPEN) {
            next_token(s);
            if(s->type != TOK_CLOSE) {
                s->type = TOK_ERROR;
            } else {
                next_token(s);
            }
        }
        break;

    case TE_FUNCTION1:
    case TE_CLOSURE1:
        ret = new_expr(s->type, 0);
        ret->function = s->function;
        if(IS_CLOSURE(s->type))
            ret->parameters[1] = s->context;
        next_token(s);
        ret->parameters[0] = power(s);
        break;

    case TE_FUNCTION2:
    case TE_FUNCTION3:
    case TE_FUNCTION4:
    case TE_FUNCTION5:
    case TE_FUNCTION6:
    case TE_FUNCTION7:
    case TE_CLOSURE2:
    case TE_CLOSURE3:
    case TE_CLOSURE4:
    case TE_CLOSURE5:
    case TE_CLOSURE6:
    case TE_CLOSURE7:
        arity = ARITY(s->type);

        ret = new_expr(s->type, 0);
        ret->function = s->function;
        if(IS_CLOSURE(s->type))
            ret->parameters[arity] = s->context;
        next_token(s);

        if(s->type != TOK_OPEN) {
            s->type = TOK_ERROR;
        } else {
            int         i;
            for(i = 0; i < arity; i++) {
                next_token(s);
                ret->parameters[i] = expr(s);
                if(s->type != TOK_SEP) {
                    break;
                }
            }
            if(s->type != TOK_CLOSE || i != arity - 1) {
                s->type = TOK_ERROR;
            } else {
                next_token(s);
            }
        }

        break;

    case TOK_OPEN:
        next_token(s);
        ret = list(s);
        if(s->type != TOK_CLOSE) {
            s->type = TOK_ERROR;
        } else {
            next_token(s);
        }
        break;

    default:
        ret = new_expr(0, 0);
        s->type = TOK_ERROR;
        ret->value = NAN;
        break;
    }

    return ret;
}

static te_expr *power(state * s) {
    /* <power>     =    {("-" | "+")} <base> */
    int         sign = 1;
    while(s->type == TOK_INFIX && (s->function == add || s->function == sub)) {
        if(s->function == sub)
            sign = -sign;
        next_token(s);
    }

    te_expr    *ret;

    if(sign == 1) {
        ret = base(s);
    } else {
        ret = NEW_EXPR(TE_FUNCTION1 | TE_FLAG_PURE, base(s));
        ret->function = negate;
    }

    return ret;
}

#ifdef TE_POW_FROM_RIGHT
static te_expr *factor(state * s) {
    /* <factor>    =    <power> {"^" <power>} */
    te_expr    *ret = power(s);

    int         neg = 0;
    te_expr    *insertion = 0;

    if(ret->type == (TE_FUNCTION1 | TE_FLAG_PURE) && ret->function == negate) {
        te_expr    *se = ret->parameters[0];
        free(ret);
        ret = se;
        neg = 1;
    }

    while(s->type == TOK_INFIX && (s->function == pow)) {
        te_fun2     t = s->function;
        next_token(s);

        if(insertion) {
            /* Make exponentiation go right-to-left. */
            te_expr    *insert =
                NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE,
                         insertion->parameters[1], power(s));
            insert->function = t;
            insertion->parameters[1] = insert;
            insertion = insert;
        } else {
            ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, power(s));
            ret->function = t;
            insertion = ret;
        }
    }

    if(neg) {
        ret = NEW_EXPR(TE_FUNCTION1 | TE_FLAG_PURE, ret);
        ret->function = negate;
    }

    return ret;
}
#else
static te_expr *factor(state * s) {
    /* <factor>    =    <power> {"^" <power>} */
    te_expr    *ret = power(s);

    while(s->type == TOK_INFIX && (s->function == pow)) {
        te_fun2     t = s->function;
        next_token(s);
        ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, power(s));
        ret->function = t;
    }

    return ret;
}
#endif

static te_expr *term(state * s) {
    /* <term>      =    <factor> {("*" | "/" | "%") <factor>} */
    te_expr    *ret = factor(s);

    while(s->type == TOK_INFIX
          && (s->function == mul || s->function == divide
              || s->function == fmod)) {
        te_fun2     t = s->function;
        next_token(s);
        ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, factor(s));
        ret->function = t;
    }

    return ret;
}

static te_expr *expr(state * s) {
    /* <expr>      =    <term> {("+" | "-") <term>} */
    te_expr    *ret = term(s);

    while(s->type == TOK_INFIX && (s->function == add || s->function == sub)) {
        te_fun2     t = s->function;
        next_token(s);
        ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, term(s));
        ret->function = t;
    }

    return ret;
}

static te_expr *list(state * s) {
    /* <list>      =    <expr> {"," <expr>} */
    te_expr    *ret = expr(s);

    while(s->type == TOK_SEP) {
        next_token(s);
        ret = NEW_EXPR(TE_FUNCTION2 | TE_FLAG_PURE, ret, expr(s));
        ret->function = comma;
    }

    return ret;
}

#define TE_FUN(...) ((double (*) (__VA_ARGS__)) n->function)
#define M(e) te_eval(n->parameters[e])

double te_eval(const te_expr * n) {
    if(!n)
        return NAN;

    switch (TYPE_MASK(n->type)) {
    case TE_CONSTANT:
        return n->value;
    case TE_VARIABLE:
        return *n->bound;

    case TE_FUNCTION0:
    case TE_FUNCTION1:
    case TE_FUNCTION2:
    case TE_FUNCTION3:
    case TE_FUNCTION4:
    case TE_FUNCTION5:
    case TE_FUNCTION6:
    case TE_FUNCTION7:
        switch (ARITY(n->type)) {
        case 0:
            return TE_FUN(void) ();
        case 1:
            return TE_FUN(double) (M(0));
        case 2:
            return TE_FUN(double, double) (M(0), M(1));
        case 3:
            return TE_FUN(double, double, double) (M(0), M(1), M(2));
        case 4:
            return TE_FUN(double, double, double, double) (M(0), M(1), M(2),
                                                           M(3));
        case 5:
            return TE_FUN(double, double, double, double, double) (M(0), M(1),
                                                                   M(2), M(3),
                                                                   M(4));
        case 6:
            return TE_FUN(double, double, double, double, double,
                          double) (M(0), M(1), M(2), M(3), M(4), M(5));
        case 7:
            return TE_FUN(double, double, double, double, double, double,
                          double) (M(0), M(1), M(2), M(3), M(4), M(5), M(6));
        default:
            return NAN;
        }

    case TE_CLOSURE0:
    case TE_CLOSURE1:
    case TE_CLOSURE2:
    case TE_CLOSURE3:
    case TE_CLOSURE4:
    case TE_CLOSURE5:
    case TE_CLOSURE6:
    case TE_CLOSURE7:
        switch (ARITY(n->type)) {
        case 0:
            return TE_FUN(void *) (n->parameters[0]);
        case 1:
            return TE_FUN(void *, double) (n->parameters[1], M(0));
        case 2:
            return TE_FUN(void *, double, double) (n->parameters[2], M(0),
                                                   M(1));
        case 3:
            return TE_FUN(void *, double, double, double) (n->parameters[3],
                                                           M(0), M(1), M(2));
        case 4:
            return TE_FUN(void *, double, double, double,
                          double) (n->parameters[4], M(0), M(1), M(2), M(3));
        case 5:
            return TE_FUN(void *, double, double, double, double,
                          double) (n->parameters[5], M(0), M(1), M(2), M(3),
                                   M(4));
        case 6:
            return TE_FUN(void *, double, double, double, double, double,
                          double) (n->parameters[6], M(0), M(1), M(2), M(3),
                                   M(4), M(5));
        case 7:
            return TE_FUN(void *, double, double, double, double, double,
                          double, double) (n->parameters[7], M(0), M(1), M(2),
                                           M(3), M(4), M(5), M(6));
        default:
            return NAN;
        }

    default:
        return NAN;
    }

}

#undef TE_FUN
#undef M

static void optimize(te_expr * n) {
    /* Evaluates as much as possible. */
    if(n->type == TE_CONSTANT)
        return;
    if(n->type == TE_VARIABLE)
        return;

    /* Only optimize out functions flagged as pure. */
    if(IS_PURE(n->type)) {
        const int   arity = ARITY(n->type);
        int         known = 1;
        int         i;
        for(i = 0; i < arity; ++i) {
            optimize(n->parameters[i]);
            if(((te_expr *) (n->parameters[i]))->type != TE_CONSTANT) {
                known = 0;
            }
        }
        if(known) {
            const double value = te_eval(n);
            te_free_parameters(n);
            n->type = TE_CONSTANT;
            n->value = value;
        }
    }
}

te_expr *te_compile(const char *expression, StrPtrMap * varmap,
                    int *error) {
    state       s;
    s.start = s.next = expression;
    s.lookup = varmap;

    next_token(&s);
    te_expr    *root = list(&s);

    if(s.type != TOK_END) {
        te_free(root);
        if(error) {
            *error = (s.next - s.start);
            if(*error == 0)
                *error = 1;
        }
        return 0;
    } else {
        optimize(root);
        if(error)
            *error = 0;
        return root;
    }
}

double te_interp(const char *expression, int *error) {
    te_expr    *n = te_compile(expression, 0, error);
    double      ret;
    if(n) {
        ret = te_eval(n);
        te_free(n);
    } else {
        ret = NAN;
    }
    return ret;
}

static void pn(const te_expr * n, int depth, FILE *fp) {
    int         i, arity;
    fprintf(fp, "%*s", depth, "");

    switch (TYPE_MASK(n->type)) {
    case TE_CONSTANT:
        fprintf(fp,"%f\n", n->value);
        break;
    case TE_VARIABLE:
        fprintf(fp,"bound %p\n", n->bound);
        break;

    case TE_FUNCTION0:
    case TE_FUNCTION1:
    case TE_FUNCTION2:
    case TE_FUNCTION3:
    case TE_FUNCTION4:
    case TE_FUNCTION5:
    case TE_FUNCTION6:
    case TE_FUNCTION7:
    case TE_CLOSURE0:
    case TE_CLOSURE1:
    case TE_CLOSURE2:
    case TE_CLOSURE3:
    case TE_CLOSURE4:
    case TE_CLOSURE5:
    case TE_CLOSURE6:
    case TE_CLOSURE7:
        arity = ARITY(n->type);
        fprintf(fp,"f%d", arity);
        for(i = 0; i < arity; i++) {
            fprintf(fp," %p", n->parameters[i]);
        }
        fprintf(fp,"\n");
        for(i = 0; i < arity; i++) {
            pn(n->parameters[i], depth + 1, fp);
        }
        break;
    }
}

void te_print(const te_expr * n, FILE *fp) {
    pn(n, 0, fp);
}

/*
 * Fill array "ptr" with pointers to the variables on which this
 * expression depends. Return the number of dependent variables, which
 * should be less than or equal to len. Abort if len is smaller than
 * the number of dependent variables.
 */
int te_dependencies(const te_expr *self, int len, const double *ptr[len]) {
    int i, arity, n=0;
    switch(TYPE_MASK(self->type)) {
    case TE_VARIABLE:
        if(len == 0) {
            fprintf(stderr,"%s:%s:%d: buffer overflow\n",
                    __FILE__,__func__,__LINE__);
            exit(EXIT_FAILURE);
        }
        ptr[0] = self->bound;
        n = 1;
        break;
    case TE_FUNCTION0:
    case TE_FUNCTION1:
    case TE_FUNCTION2:
    case TE_FUNCTION3:
    case TE_FUNCTION4:
    case TE_FUNCTION5:
    case TE_FUNCTION6:
    case TE_FUNCTION7:
    case TE_CLOSURE0:
    case TE_CLOSURE1:
    case TE_CLOSURE2:
    case TE_CLOSURE3:
    case TE_CLOSURE4:
    case TE_CLOSURE5:
    case TE_CLOSURE6:
    case TE_CLOSURE7:
        arity = ARITY(self->type);
        for(i=0; i < arity; ++i)
            n += te_dependencies(self->parameters[i], len-n, ptr+n);
        break;
    default:
        break;
    }
    return n;
}
