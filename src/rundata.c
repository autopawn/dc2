#include "rundata.h"

const char *filter_names[] = {
    "NO_FILTER",
    "BETTER_THAN_EMPTY",
    "BETTER_THAN_ONE_PARENT",
    "BETTER_THAN_ALL_PARENTS",
    "BETTER_THAN_SUBSETS",
};

const char *local_search_names[] = {
    "NO_LOCAL_SEARCH",
    "SWAP_BEST_IMPROVEMENT",
    "SWAP_FIRST_IMPROVEMENT",
    "SWAP_RESENDE_WERNECK",
};

const char *path_relinking_names[] = {
    "NO_PATH_RELINKING",
    "PATH_RELINKING_1_ITER",
    "PATH_RELINKING_UNTIL_NO_BETTER",
};

void rundata_free(rundata *data){
    // Free precomputations
    runprecomp_free(data->precomp);
    // Free run info
    runinfo_free(data->run_inf);
    // Free problem
    problem_free(data->prob);
    // Free rundata
    free(data);
}



// Creates a rundata for the given problem and performs precomputations
rundata *rundata_init(problem *prob, redstrategy *rstrats, int n_rstrats, int n_restarts, int precomp_nearly_indexes, int n_threads, int verbose){
    rundata *run = safe_malloc(sizeof(rundata));

    run->prob = problem_copy(prob);
    prob = run->prob;

    run->filter = FILTER_DEFAULT;
    run->branch_and_bound = BRANCH_AND_BOUND_DEFAULT;
    run->random_seed = 42;
    run->n_restarts = n_restarts;

    // Default values
    run->bnb_lower_bound = -INFINITY;
    run->verbose     = verbose;
    run->branching_factor     = DEFAULT_BRANCHING_FACTOR;
    run->branching_correction = DEFAULT_BRANCHING_CORRECTION;
    run->path_relinking   = DEFAULT_PATH_RELINKING;

    run->target_sols  = DEFAULT_TARGET_SOLS;
    run->n_threads    = n_threads;
    run->local_search = DEFAULT_LOCAL_SEARCH;
    run->local_search_pr = DEFAULT_LOCAL_SEARCH;
    run->local_search_before_select = DEFAULT_LOCAL_SEARCH_BEFORE_SELECT;
    run->select_only_terminal = DEFAULT_SELECT_ONLY_TERMINAL;
    run->local_search_rem_movement = DEFAULT_LOCAL_SEARCH_SIZE_CHANGE_MOVEMENTS_ENABLED;
    run->local_search_add_movement = DEFAULT_LOCAL_SEARCH_SIZE_CHANGE_MOVEMENTS_ENABLED;

    // Initialize runinfo
    run->run_inf = runinfo_init(prob,n_restarts);

    // Initialize precomputations and perform them
    run->precomp = runprecomp_init(prob,rstrats,n_rstrats,precomp_nearly_indexes,run->n_threads,run->verbose);

    return run;
}

// Prints a briefing of the rundata parameters
void rundata_print(const rundata *run, FILE *fp){
    const problem *prob = run->prob;

    fprintf(fp,"== PROBLEM ==\n");
    fprintf(fp,"# N_FACILITIES: %d\n",prob->n_facs);
    fprintf(fp,"# N_CLIENTS: %d\n",prob->n_clis);
    fprintf(fp,"# SIZE_RESTRICTION_MINIMUM: %d\n",prob->size_restriction_minimum);
    fprintf(fp,"# SIZE_RESTRICTION_MAXIMUM: %d\n",prob->size_restriction_maximum);
    fprintf(fp,"\n");
    fprintf(fp,"== RUN DATA ==\n");
    fprintf(fp,"# FILTER: %s (%d)\n",filter_names[run->filter],(int)run->filter);
    fprintf(fp,"# TARGET_SOLS: %d\n",run->target_sols);
    fprintf(fp,"# PRECOMP_CLIENT_OPTIMAL_GAIN: %lf\n",run->precomp->precomp_client_optimal_gain);
    fprintf(fp,"# N_THREADS: %d\n",run->n_threads);
    fprintf(fp,"# SELECT_ONLY_TERMINAL: %d\n",run->select_only_terminal);
    fprintf(fp,"# LOCAL_SEARCH: %s\n",local_search_names[run->local_search]);
    fprintf(fp,"# LOCAL_SEARCH_BEFORE_SELECT: %d\n",run->local_search_before_select);
    fprintf(fp,"# LOCAL_SEARCH_REM_MOVEMENT: %d\n",run->local_search_rem_movement);
    fprintf(fp,"# LOCAL_SEARCH_ADD_MOVEMENT: %d\n",run->local_search_add_movement);
    fprintf(fp,"# BRANCHING_FACTOR: %d\n",run->branching_factor);
    fprintf(fp,"# BRANCHING_FACTOR_CORRECTION: %d\n",run->branching_correction);
    fprintf(fp,"# PATH_RELINKING: %s\n",path_relinking_names[run->path_relinking]);
    fprintf(fp,"# PATH_RELINKING_LOCAL_SEARCH: %s\n",local_search_names[run->local_search_pr]);
    fprintf(fp,"# RANDOM_SEED: %d\n",run->random_seed);
    fprintf(fp,"# RESTARTS: %d\n",run->n_restarts);
    fprintf(fp,"# VERBOSE: %d\n",run->verbose);
    fprintf(fp,"# BRANCH_AND_BOUND: %d (deprecated)\n",run->branch_and_bound);
    fprintf(fp,"# LOWER_BOUND: %lf \n",run->bnb_lower_bound);
}
