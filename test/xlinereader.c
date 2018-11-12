/**
 * @file xlinereader.c
 * @author Alan R. Rogers
 * @brief Test linereader.c
 * @copyright Copyright (c) 2018, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */

#include "linereader.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef NDEBUG
#error "Unit tests must be compiled without -DNDEBUG flag"
#endif

const char *tstInput =
    "#chr  pos  ref  alt  raf\n"
    "1  1  a  t  0\n"
    "onefield";

int main(int argc, char **argv) {
    size_t      buffsize = 2;

    if(argc != 1) {
        fprintf(stderr,"usage: xlinereader\n");
        exit(EXIT_FAILURE);
    }

    const char *fname = "linereader.tmp";
    FILE *fp = fopen(fname,"w");
    fputs(tstInput, fp);
    fclose(fp);
    fp = fopen(fname,"r");

    LineReader *lr = LineReader_new(buffsize);
    assert(lr);


    char *buff;

    buff = LineReader_next(lr, fp);
    assert(buff != NULL);
    assert(0 == strcmp(buff, "#chr  pos  ref  alt  raf\n"));

    buff = LineReader_next(lr, fp);
    assert(buff != NULL);
    assert(0 == strcmp(buff, "1  1  a  t  0\n"));

    buff = LineReader_next(lr, fp);
    assert(buff != NULL);
    assert(0 == strcmp(buff, "onefield"));

    assert(feof(fp));
    assert(NULL == LineReader_next(lr, fp));

    LineReader_free(lr);
    fclose(fp);
    remove(fname);
    unitTstResult("LineReader", "OK");

    return 0;
}
