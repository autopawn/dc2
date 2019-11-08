#ifndef DC_CONSTRUCTION_H
#define DC_CONSTRUCTION_H

#include "utils.h"
#include "problem.h"
#include "solution.h"
#include "reduction.h"
#include "expand.h"
#include "bnb.h"

solution **new_find_best_solutions(problem *prob, redstrategy *rstrats, int n_rstrats,
        int *out_n_sols, int *out_n_iterations);

#endif