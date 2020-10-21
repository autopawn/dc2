#ifndef DC2_RUNINFO_H
#define DC2_RUNINFO_H

#include "problem.h"

typedef struct {
    // | Total number of iterations performed. Just the first restart when using restarts.
    int firstr_n_iterations;
    // | Number of solutions after expansion. Just the first restart when using restarts.
    int *firstr_per_size_n_sols;
    // | Number of solutions after reduction. Just the first restart when using restarts.
    int *firstr_per_size_n_sols_after_red;
    // | Number of different local optima found. Just the first restart when using restarts.
    int *firstr_per_size_n_local_optima;
    // | Total number of iterations
    int total_n_iterations;
    // | Total number of local searches
    long long int n_local_searches;
    // | Total number of local search movements
    long long int n_local_search_movements;
    // | CPU time performing local search:
    double local_search_seconds;
    // | Time taken on each restart
    double *restart_times;
    // | Values on each restart
    double *restart_values;
    // | CPU time performing path relinking:
    double path_relinking_seconds;
} runinfo;

runinfo *runinfo_init(const problem *prob, int n_restarts);
void runinfo_free(runinfo *rinf);

#endif
