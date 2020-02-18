#ifndef DC_EXPAND_H
#define DC_EXPAND_H

#include "utils.h"
#include "solution.h"

solution **new_expand_solutions(const rundata *run,
        solution **sols, int n_sols, int *out_n_sols);

#endif