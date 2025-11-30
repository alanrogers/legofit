/**
 * @file ptru32map.c
 * @author Alan R. Rogers
 * @brief Map pointer to unsigned int
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "ptru32map.h"
#include "binary.h"
#include "misc.h"
#include "error.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "hashmap.c"
