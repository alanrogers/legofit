/**
@file rafmain.c
@page rafmain
@brief Wrapper for raf function.

This file defines a `main` function that does nothing but call `raf`
and return a status code. The excutable produced by this file takes
no arguments, reads from standard input, and writes to standard output.
See documentation for `raf`.
**/

#include "raf.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {

    if(argc != 1) {
        fprintf(stderr, "Usage: raf\n");
        fprintf(stderr, "       Reads standard input;"
                " writes to standard output.\n");
        exit(EXIT_FAILURE);
    }

    int status = raf(stdin, stdout);

    return status==0 ? 0 : EXIT_FAILURE;
}
