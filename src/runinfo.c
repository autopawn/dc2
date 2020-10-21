#include "runinfo.h"

runinfo *runinfo_init(const problem *prob, int n_restarts){
    runinfo *rinf = safe_malloc(sizeof(runinfo));

    rinf->firstr_n_iterations      = 0;
    rinf->total_n_iterations       = 0;
    rinf->n_local_searches         = 0;
    rinf->n_local_search_movements = 0;
    rinf->local_search_seconds     = 0;
    rinf->path_relinking_seconds   = 0;

    // First restart data
    rinf->firstr_per_size_n_sols = safe_malloc(sizeof(int)*(prob->n_facs+2));
    memset(rinf->firstr_per_size_n_sols,0,     sizeof(int)*(prob->n_facs+2));
    rinf->firstr_per_size_n_sols_after_red = safe_malloc(sizeof(int)*(prob->n_facs+2));
    memset(rinf->firstr_per_size_n_sols_after_red,0,     sizeof(int)*(prob->n_facs+2));
    rinf->firstr_per_size_n_local_optima = safe_malloc(sizeof(int)*(prob->n_facs+2));
    memset(rinf->firstr_per_size_n_local_optima,0,     sizeof(int)*(prob->n_facs+2));

    // Restart data
    rinf->restart_times  = safe_malloc(sizeof(double)*n_restarts);
    rinf->restart_values = safe_malloc(sizeof(double)*n_restarts);
    for(int r=0;r<n_restarts;r++){
        rinf->restart_times[r] = NAN;
        rinf->restart_values[r] = -INFINITY;
    }

    return rinf;
}

void runinfo_free(runinfo *rinf){
    // Free first restart data
    free(rinf->firstr_per_size_n_sols);
    free(rinf->firstr_per_size_n_sols_after_red);
    free(rinf->firstr_per_size_n_local_optima);
    // Free restart data
    free(rinf->restart_times);
    free(rinf->restart_values);
    //
    free(rinf);
}
