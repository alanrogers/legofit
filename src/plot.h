/**
Generic code for header for plot function. To create a concrete
version, include this into a .h file that defines NODETYPE.
**/

#include "sampndx.h"

#define INNER(CLASS, NAME) CLASS ## _ ## NAME
#define FUNC(CLASS, NAME) INNER(CLASS, NAME)

// PLOT expands to something like Foo_plot, where "Foo" is the value
// of the NODETYPE macro.
#define PLOT FUNC(NODETYPE, plot)

void PLOT(NODETYPE *root, SampNdx *sndx, FILE *fp);

#undef INNER
#undef FUNC
#undef PLOT
