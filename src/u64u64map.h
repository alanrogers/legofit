#ifndef U64U64MAP_H
#  define U64U64MAP_H
#  include "typedefs.h"
#  include <stdio.h>
#  include <stdint.h>

U64U64Map *U64U64Map_new(void);
void       U64U64Map_free(U64U64Map * self);
int        U64U64Map_get(U64U64Map * self, uint64_t key, uint64_t *value);
int        U64U64Map_insert(U64U64Map * self, uint64_t key, uint64_t value);
unsigned long U64U64Map_size(U64U64Map * self);
void       U64U64Map_print(U64U64Map * self, FILE *fp);

#endif
