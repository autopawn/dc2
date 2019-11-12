#ifndef DC_OUTPUT_H
#define DC_OUTPUT_H

#include "utils.h"
#include "problem.h"
#include "solution.h"
#include <sys/time.h>

// Convenience function to get the delta in seconds between two timevals
float get_delta_seconds(struct timeval tv1, struct timeval tv2);

// Saves the results of an execution, with the problem information and
// the resulting solutions
void save_solutions(const char *file, 
        const problem *prob, solution **sols, int n_sols,
        const char *input_file, float seconds, float elapsed, int n_iters,
        const redstrategy *strategies, int n_strategies);


#endif