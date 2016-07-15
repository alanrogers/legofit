#include "vcfreader.h"
#include "tokenizer.h"
#include "misc.h"
#include <string.h>

#define VCF_MAXFIELDS 200
#define VCF_MAXLINE 500

struct VCFReader {
    FILE *fp;
    char buff[VCF_MAXLINE];
    char *reference;
    Tokenizer *tkz;
    unsigned long snpid;
};

VCFReader *VCFReader_new(const char *fname) {
    if(fp == NULL) {
        fprintf(stderr,"%s:%s:%d: input stream is NULL\n",
                __FILE__,__func__,__LINE__);
        exit(EXIT_FAILURE);
    }
    VCFReader *self = malloc(sizeof(*self));
    CHECKMEM(self);
    self->fname = strdup(fname);
    CHECKMEM(self->fname);
    self->fp = fopen(self->fname);
    if(self->fp == NULL) {
        fprintf(stderr,"%s:%s:%d: can't open \"%s\" for input.\n",
                __FILE__,__func__,__LINE__, self->fname);
        exit(EXIT_FAILURE);
    }
    self->tkz = Tokenizer_new(VCF_MAXFIELDS);
    VCFReader_parseHdr(self);
    return self;
}

void VCFReader_free(VCFReader *self) {
    fclose(self->fp);
    free(self->fname);
    free(self->reference);
    Tokenizer_free(self->tkz);
}

void VCFReader_parseHdr(VCFReader *self) {
    while(1) {
        if(fgets(self->buff, sizeof(self->buff), self->fp) == NULL)
            break;
        if(strcmp(buff, "##reference") == 0){
            Tokenizer_split(self->tkz, self->buff, "=");
            Tokenizer_strip(self->tkz, " \t\n");
            assert(Tokenizer_ntokens(self->tkz) == 2);
            self->reference = strdup(Tokenizer_token(self->tkz,1));
            CHECKMEM(self->reference);
        }else if(buff[0] == '#')
            continue;
        else if(NULL != strchr(',', self->buff))
            continue;
        else
            break;
    }
    self->snpid = 0;
}

