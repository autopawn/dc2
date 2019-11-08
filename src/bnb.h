#ifndef DC_BNB_H
#define DC_BNB_H

#include "utils.h"
#include "problem.h"
#include "solution.h"

void branch_and_bound(problem *prob, solution **sols, int *n_sols);

#endif