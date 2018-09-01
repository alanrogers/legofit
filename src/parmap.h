/**
   This file is meant to be included by other files, which define macros
   METHOD, CLASS, and KEYTYPE. It can be included by more than one such
   file in order to define more than one class.
 */

#include "typedefs.h"
#include <stdio.h>

Param      *METHOD(search)(CLASS *root, KEYTYPE key);
CLASS      *METHOD(insert)(CLASS *root, Param *par);
void        METHOD(free)(CLASS *h);
void        METHOD(print)(CLASS *h, FILE *fp, int indent);

