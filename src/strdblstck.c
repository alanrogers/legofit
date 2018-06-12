#include "hessian.h"
#include "misc.h"
#include "strdblstck.h"

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
