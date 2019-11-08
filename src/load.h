#ifndef DC_LOAD_H
#define DC_LOAD_H

#include "utils.h"
#include "problem.h"

// Loads a problem from a given file and performs precomputations.
problem *new_problem_load(const char *file);

#endif