/**
@file mkfitted.c
@page mkfitted
@author Alan R. Rogers
@brief Make a fitted .lgo file from .legofit output

# mkfitted: make a fitted .lgo file from .legofit output

usage: mkfitted <original.lgo> <fitted.legofit>

@copyright Copyright (c) 2023, Alan R. Rogers
<rogers@anthro.utah.edu>. This file is released under the Internet
Systems Consortium License, which can be found in file "LICENSE".
*/

#include "strdblqueue.h"
#include "strdblmap.h"
#include "misc.h"
#include <strings.h>
#include <time.h>
#include <ctype.h>

#define BUFSIZE 4096
#define VARNAMESIZE 512

struct LineReader {
    char buff[BUFSIZE];
    char varname[VARNAMESIZE];
    char *eqpos; // position of "=" (if any) w/i buff
};

enum linetype {assign, noassign, eof};

LineReader *LineReader_new(void);
enum linetype LineReader_next(LineReader *self, FILE *fp);
void LineReader_free(LineReader *self);
int LineReader_pr_orig(LineReader *self);
int LineReader_pr_newval(LineReader *self, double val);

// Function prototypes
void usage(void);

// Print usage message and exit.
void usage(void) {
    fputs("usage: mkfitted <original.lgo> <fitted.legofit>\n",
          stderr);
    exit(EXIT_FAILURE);
}

void LineReader_free(LineReader *self) {
    free(self);
}

LineReader *LineReader_new(void) {
    LineReader *self = malloc(sizeof(LineReader));
    CHECKMEM(self);
    self->buff[0] = '\0';
    self->eqpos = NULL;
    return self;
};

/// Read next line. Return 0 if input is not an assignment statement,
/// 1 if it is an assignment, or EOF on end of file.
enum linetype LineReader_next(LineReader *self, FILE *fp) {
    self->buff[0] = '\0';
    self->eqpos = NULL;

    char *s = self->buff;
    int c;
    int comment = 0;
    int operator = 0;
    char *eqpos = NULL;
    enum linetype ltype = noassign;

    // Outer do-while loop continues until we have accumulated a
    // complete logical line. A logical line may include several
    // physical lines, provided that the incomplete physical lines end
    // with an operator (+-*/).  In deciding whether a physical line
    // ends with an operator, the algorithm ignores whitespace and
    // #-initiated comments at the end of the physical line.
    do{
        // Each pass through the while loop handles a single character.
        while(1){
            c = getc(fp);
            switch(c) {
            case EOF:
                if(operator) {
                    fprintf(stderr, "%s:%d: .lgo file must not end with an"
                            " operator: +,-,*,/\n",
                            __FILE__,__LINE__);
                    exit(EXIT_FAILURE);
                }
                goto end_physical_line;
            case '#':
                comment = 1;
                break;
            case '=':
                if( !comment )
                    eqpos = s;
                break;
            case '\n':
                comment = 0;
                *s++ = c;
                goto end_physical_line;
            default:
                if(!comment) {
                    if( strchr("+-*/", c) )
                        operator = 1;
                    else if( !isspace(c) )
                        operator = 0;
                }
                break;
            } // end switch
            if(s - self->buff >= BUFSIZE) {
                fprintf(stderr,"%s:%d: buffer overflow\n",
                        __FILE__,__LINE__);
                exit(EXIT_FAILURE);
            }
            *s++ = c;
        }
    end_physical_line:
        continue;
    }while(operator);
    if(s - self->buff >= BUFSIZE) {
        fprintf(stderr,"%s:%d: buffer overflow\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    *s = '\0';

    if(s == self->buff) { // no input
        return eof;
    }

    // Classify statement
    if( 0 == strncmp(self->buff, "time", 4)
        || 0 == strncmp(self->buff, "twoN", 4)
        || 0 == strncmp(self->buff, "mixFrac", 7)
        || 0 == strncmp(self->buff, "param", 5) ) {
        ltype = assign;
    }else{
        ltype = noassign;
    }

    if(ltype == noassign)
        return ltype;

    // We've read an assignment statement.

    // Copy variable name into self->varname.
    char *end = eqpos;
    while(end > self->buff && isspace(*(end-1)) )
        --end;
    char *name = end;
    while(name > self->buff && !isspace(*(name - 1)) )
        --name;
    int len = end - name;
    if(len >= VARNAMESIZE) {
        fprintf(stderr,"%s:%d: buffer overflow\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }
    int status = snprintf(self->varname, len+1, "%s", name);
    assert(status >= 0);
    if( !legalName(self->varname) ) {
        fprintf(stderr, "%s:%d: %s is not a legal variable name\n",
                __FILE__,__LINE__, self->varname);
        exit(EXIT_FAILURE);
    }

    self->eqpos = eqpos;

    return assign;
}

/// Print input line unchanged. Return non-negative integer on success,
/// EOF on failure.
int LineReader_pr_orig(LineReader *self) {
    int status =  fputs(self->buff, stdout);
    if(status == EOF)
        return EOF;
    return 0;
}

/// Print assignment statement with new value. Return non-negative integer
/// on success, EOF on failure.
int LineReader_pr_newval(LineReader *self, double val) {
    int status, eol=0, lastchr=0;
    char *s;

    if(self->eqpos == NULL) {
        fprintf(stderr,"%s:%d: can't rewrite line because it's not an"
                " assignment statement\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    // Print up to and including '=' character.
    for(s=self->buff; s < self->buff + BUFSIZE; ++s) {
        lastchr = *s;
        status = putchar(*s);
        if(status == EOF)
            return EOF;
        if(*s == '\0') {
            fprintf(stderr,"%s:%d: line ended too soon\n",
                    __FILE__,__LINE__);
            exit(EXIT_FAILURE);
        }
        if( s == self->eqpos ) {
            ++s;
            break;
        }
    }

    // Print whitespace following '=' character.
    while(s < self->buff + BUFSIZE && isspace(*s)) {
        lastchr = *s;
        status = putchar(*s++);
        if(status == EOF)
            return EOF;
    }

    // Print new value
    printf("%lg", val);

    // Print comments if any. Outer loop continues until we
    // hit a '\0' character or exhaust the buffer.
    while(s < self->buff + BUFSIZE) {
        // skip until next '#' character
        while(s < self->buff + BUFSIZE && *s != '#') {
            if(*s == '\0')
                goto pr_newval_done;
            if(*s == '\n') {
                eol = 1;
                if(lastchr != '\n') {
                    lastchr = '\n';
                    putchar('\n');
                }
            }
            ++s;
        }

        // space between value and comment
        if(s < self->buff + BUFSIZE && *s == '#'){
            lastchr = ' ';
            putchar(' ');
        }
        
        // print until end of line.
        while(s < self->buff + BUFSIZE && !eol) {
            if(*s == '\0')
                goto pr_newval_done;
            if(*s == '\n') {
                eol = 1;
            }
            if(*s != '\n' || lastchr != '\n') {
                lastchr = *s;
                putchar(*s++);
            }
        }
    }
 pr_newval_done:
    return 0;
}

int main(int argc, char **argv) {

    time_t currtime = time(NULL);
    int status, i;
    int num_rewritten = 0;

    hdr("mkfitted: make a fitted .lgo file from legofit output");
#if defined(__DATE__) && defined(__TIME__)
    printf("# Program was compiled: %s %s\n", __DATE__, __TIME__);
#endif
    printf("# Program was run: %s", ctime(&currtime));
    printf("# Cmd:");
    for(i=0; i<argc; ++i)
        printf(" %s", argv[i]);
    putchar('\n');

    if(argc != 3)
        usage();
    for(i=1; i<argc; ++i) {
        if(argv[i][0] == '-')
            usage();
    }

    const char *lgofname = argv[1]; // .lgo file to be rewritten
    const char *fitfname = argv[2]; // .legofit file

    // Read .legofit file to create a queue of key-value pairs,
    // which contain the fitted parameter names and values.
    StrDblQueue *fitqueue = StrDblQueue_parseLegofit(fitfname, 0);
    if(fitqueue == NULL) {
        fprintf(stderr, "%s:%d: can't parse file \"%s\" as legofit output.\n",
                __FILE__, __LINE__, fitfname);
        exit(EXIT_FAILURE);
    }
    int queue_len = StrDblQueue_length(fitqueue);
    //StrDblQueue_print(fitqueue, stdout);

    StrDblMap *parmap = StrDblMap_new(128);
    if(parmap == NULL) {
        fprintf(stderr, "%s:%d: can't allocate StrDblMap.\n",
                __FILE__,__LINE__);
        exit(EXIT_FAILURE);
    }

    // Use fitqueue to initialize parmap
    while(fitqueue != NULL) {
        StrDbl strdbl;
        fitqueue = StrDblQueue_pop(fitqueue, &strdbl);
        status = StrDblMap_insert(parmap, strdbl.str, strdbl.val);
        switch(status) {
        case 0:
            continue;
        case 1:
            fprintf(stderr,"%s:%d: can't insert same parameter %s twice.\n",
                    __FILE__,__LINE__,strdbl.str);
            exit(EXIT_FAILURE);
        case ENOMEM:
            fprintf(stderr,"%s:%d: can't resize hash table.\n",
                    __FILE__,__LINE__);
            exit(EXIT_FAILURE);
        default:
            fprintf(stderr,"%s:%d: this should not happen.\n",
                    __FILE__,__LINE__);
            exit(EXIT_FAILURE);
        }
    }

    // Use parmap to edit .lgo file
    FILE *lgofile = fopen(lgofname, "r");
    if(lgofile == NULL) {
        fprintf(stderr, "%s:%d: can't read file %s\n",
                __FILE__,__LINE__, lgofname);
        exit(EXIT_FAILURE);
    }
    double val;
    LineReader *lr = LineReader_new();
    enum linetype ltype;
    do{
        ltype = LineReader_next(lr, lgofile);
        switch(ltype) {
        case eof:
            break;
        case noassign:
            LineReader_pr_orig(lr);
            break;
        case assign:
            val = StrDblMap_get(parmap, lr->varname, &status);
            if(status == 0) {
                // variable gets new value
                num_rewritten += 1;
                LineReader_pr_newval(lr, val);
            }else{
                // variable not in parmap: leave unchanged
                LineReader_pr_orig(lr);
            }
            break;
        default:
            fprintf(stderr,"%s:%d: this shouldn't happen. ltype=%d\n",
                    __FILE__,__LINE__, ltype);
            exit(EXIT_FAILURE);
        }
    }while(ltype != eof);

    fprintf(stderr, "mkfitted: Read %d fitted parameter values.\n",
            queue_len);
    fprintf(stderr,"mkfitted: Rewrote %d assignment statements.\n",
            num_rewritten);

    fclose(lgofile);
    return 0;
}
