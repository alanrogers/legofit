#ifndef POPNODETAB_H
#  define POPNODETAB_H
#  include "typedefs.h"

PopNodeTab *PopNodeTab_new(void);
void        PopNodeTab_free(PopNodeTab * self);
PopNode    *PopNodeTab_get(PopNodeTab * self, const char *key);
int         PopNodeTab_insert(PopNodeTab * self, const char *key,
                              PopNode * node);
unsigned long PopNodeTab_size(PopNodeTab * self);
void        PopNodeTab_print(PopNodeTab * self);
PopNode    *PopNodeTab_check_and_root(PopNodeTab *self, const char *file,
                                      int line);

#endif
