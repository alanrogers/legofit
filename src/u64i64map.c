/**
 * @file u64i64map.c
 * @author Alan R. Rogers
 * @brief Map uint64_t to int64_t
 * @copyright Copyright (c) 2020, Alan R. Rogers
 * <rogers@anthro.utah.edu>. This file is released under the Internet
 * Systems Consortium License, which can be found in file "LICENSE".
 */
#include "u64i64map.h"
#include "binary.h"
#include "misc.h"
#include "error.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "hashmap.c"
