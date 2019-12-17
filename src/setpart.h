#ifndef ARR_SETPART_H
#define ARR_SETPART_H

int traverseSetPartitions(unsigned nelements, unsigned nsubdivisions,
                          int (*visit)(unsigned n, unsigned a[n],
                                       void *data), void *data);

#endif
