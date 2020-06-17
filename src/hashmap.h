/**
Generic code for the include file of a hash map. To create a concrete
version, include this into a .h file that defines MAPTYPE, KEYTYPE,
and VALTYPE.
**/

#define INNER(CLASS, NAME) CLASS ## _ ## NAME
#define FUNC(CLASS, NAME) INNER(CLASS, NAME)

MAPTYPE * FUNC(MAPTYPE, new)(void);
void      FUNC(MAPTYPE, free)(MAPTYPE * self);
int       FUNC(MAPTYPE, get)(MAPTYPE * self, KEYTYPE key, VALTYPE *value);
int       FUNC(MAPTYPE, insert)(MAPTYPE * self, KEYTYPE key, VALTYPE value);
unsigned long FUNC(MAPTYPE, size)(MAPTYPE * self);
int       FUNC(MAPTYPE, keys)(MAPTYPE *self, unsigned size, KEYTYPE keys[size]);
