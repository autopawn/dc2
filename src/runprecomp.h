#ifndef DC2_RUNPRECOMP_H
#define DC2_RUNPRECOMP_H

#include "problem.h"
#include "redstrategy.h"

typedef struct {
    // | Number of facilitites and client to keep the struct independent.
    int n_facs, n_clis;
    // | Precomputed optimal gain from all clients
    double precomp_client_optimal_gain;
    // | Precomputed distance matrices between facilities (for each mode)
    double **facs_distance[N_FACDIS_MODES];
    // | Precomputed facility indexes by proximity for each client for Resende and Werneck's local search
    int **nearly_indexes;
} runprecomp;

runprecomp *runprecomp_init(const problem *prob, redstrategy *rstrats, int n_rstrats, int precomp_nearly_indexes, int n_threads, int verbose);

void runprecomp_free(runprecomp *pcomp);

#endif
