#include "problem.h"

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
};

problem *problem_init(int n_facs, int n_clis){
    problem *prob = safe_malloc(sizeof(problem));
    prob->n_facs = n_facs;
    prob->n_clis = n_clis;
    //
    prob->client_weight = safe_malloc(sizeof(double)*prob->n_clis);
    memset(prob->client_weight,0,     sizeof(double)*prob->n_clis);
    //
    prob->facility_cost = safe_malloc(sizeof(double)*prob->n_facs);
    memset(prob->facility_cost,0,     sizeof(double)*prob->n_facs);
    //
    prob->distance = safe_malloc(sizeof(double*)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++){
        prob->distance[i] = safe_malloc(sizeof(double)*prob->n_clis);
        memset(prob->distance[i],0,     sizeof(double)*prob->n_clis);
    }

    prob->size_restriction_minimum = -1;
    prob->size_restriction_maximum = -1;
    prob->transport_cost = 1;
    prob->client_gain = 0;
    prob->unassigned_cost = INFINITY;

    prob->verbose = 1;

    return prob;
}

problem *problem_copy(const problem *other){
    problem *prob = problem_init(other->n_facs,other->n_clis);
    //
    prob->transport_cost = other->transport_cost;
    prob->client_gain = other->client_gain;
    prob->unassigned_cost = other->unassigned_cost;
    prob->size_restriction_minimum = other->size_restriction_minimum;
    prob->size_restriction_maximum = other->size_restriction_maximum;
    prob->verbose = other->verbose;
    //
    memcpy(prob->client_weight,other->client_weight,sizeof(double)*prob->n_clis);
    //
    memcpy(prob->facility_cost,other->facility_cost,sizeof(double)*prob->n_facs);
    //
    for(int i=0;i<prob->n_facs;i++){
        for(int j=0;j<prob->n_clis;j++){
            prob->distance[i][j] = other->distance[i][j];
        }
    }
    return prob;
}

void problem_free(problem *prob){
    // Free facility-client distances
    for(int i=0;i<prob->n_facs;i++){
        free(prob->distance[i]);
    }
    free(prob->distance);
    // Free per facility and client arrays
    free(prob->facility_cost);
    free(prob->client_weight);
    // Free problem
    free(prob);
}

void rundata_free(rundata *data){
    // Free facility-facility distance matrices, if they are in use:
    for(int mode=0;mode<N_FACDIS_MODES;mode++){
        if(data->facs_distance[mode]){
            for(int i=0;i<data->prob->n_facs;i++){
                free(data->facs_distance[mode][i]);
            }
            free(data->facs_distance[mode]);
        }
    }

    // Free first restart data
    free(data->firstr_per_size_n_sols);
    free(data->firstr_per_size_n_sols_after_red);
    free(data->firstr_per_size_n_local_optima);
    // Free restart data
    free(data->restart_times);
    free(data->restart_values);
    // Free problem
    problem_free(data->prob);
    // Free rundata
    free(data);
}

typedef struct {
    int thread_id;
    rundata *run;
    int mode;
} precomp_facs_dist_thread_args;

void *precomp_facs_dist_thread_execution(void *arg){
    precomp_facs_dist_thread_args *args = (precomp_facs_dist_thread_args *) arg;
    const problem *prob = args->run->prob;
    // Compute facility-facility distances acording to mode
    for(int a=args->thread_id;a<prob->n_facs;a+=args->run->n_threads){
        for(int b=a;b<prob->n_facs;b++){
            double dist = 0;
            if(args->mode==FACDIS_SUM_OF_DELTAS){
                dist = 0;
                for(int j=0;j<prob->n_clis;j++){
                    double delta = problem_assig_value(prob,a,j) - problem_assig_value(prob,b,j);
                    if(delta<0) delta = -delta;
                    dist += delta;
                }
            }else if(args->mode==FACDIS_MIN_TRIANGLE){
                dist = INFINITY;
                for(int j=0;j<prob->n_clis;j++){
                    double dist_sum = prob->distance[a][j]+prob->distance[b][j];
                    if(dist_sum<dist) dist = dist_sum;
                }
            }else{
                assert(!"Asked to precompute valid facility distance.");
            }
            args->run->facs_distance[args->mode][a][b] = dist;
            args->run->facs_distance[args->mode][b][a] = dist;
        }
    }
    return NULL;
}

void rundata_precompute(rundata *run, redstrategy *rstrats, int n_rstrats){
    const problem *prob = run->prob;
    // Precompute value of empty solution
    double empty_val = 0;
    for(int j=0;j<prob->n_clis;j++){
        empty_val += problem_assig_value(prob,-1,j);
    }
    run->precomp_empty_value = empty_val;
    // Update lower bound
    if(empty_val > run->lower_bound) run->lower_bound = empty_val;

    // Precompute precomp_client_optimal_gain
    run->precomp_client_optimal_gain = 0;
    for(int k=0;k<prob->n_clis;k++){
        double best_val = problem_assig_value(prob,-1,k);
        for(int i=0;i<prob->n_facs;i++){
            double val = problem_assig_value(prob,i,k);
            if(best_val < val) best_val = val;
        }
        run->precomp_client_optimal_gain += best_val;
    }

    // Precompute facility distances
    for(int r=0;r<n_rstrats;r++){
        int mode = redstrategy_required_facdis_mode(rstrats[r]);
        if(mode!=FACDIS_NONE && run->facs_distance[mode]==NULL){
            // Allocate distance matrix between facilities
            run->facs_distance[mode] = safe_malloc(sizeof(double*)*prob->n_facs);
            for(int i=0;i<prob->n_facs;i++){
                run->facs_distance[mode][i] = safe_malloc(sizeof(double)*prob->n_facs);
            }
            // Allocate memory for threads and arguments
            pthread_t *threads = safe_malloc(sizeof(pthread_t)*run->n_threads);
            precomp_facs_dist_thread_args *targs = safe_malloc(sizeof(precomp_facs_dist_thread_args)*run->n_threads);
            // Call threads to compute facility-facility distances
            for(int i=0;i<run->n_threads;i++){
                targs[i].thread_id = i;
                targs[i].run = run;
                targs[i].mode = mode;
                int rc = pthread_create(&threads[i],NULL,precomp_facs_dist_thread_execution,&targs[i]);
                if(rc){
                    fprintf(stderr,"ERROR: Error %d on pthread_create\n",rc);
                    exit(1);
                }
            }
            // Join threads
            for(int i=0;i<run->n_threads;i++){ // Join threads
                pthread_join(threads[i],NULL);
            }
            // Free memory
            free(targs);
            free(threads);
        }
    }
}

// Creates a rundata for the given problem and performs precomputations
rundata *rundata_init(problem *prob, redstrategy *rstrats, int n_rstrats, int n_restarts){
    rundata *run = safe_malloc(sizeof(rundata));

    run->prob = problem_copy(prob);
    prob = run->prob;

    run->filter = FILTER_DEFAULT;
    run->branch_and_bound = 1;
    run->random_seed = 42;

    run->firstr_n_iterations      = 0;
    run->total_n_iterations       = 0;
    run->n_local_searches         = 0;
    run->n_local_search_movements = 0;
    run->local_search_seconds     = 0;

    // Default values
    run->lower_bound = -INFINITY;

    run->target_sols  = DEFAULT_TARGET_SOLS;
    run->n_threads    = DEFAULT_THREADS;
    run->local_search = DEFAULT_LOCAL_SEARCH;
    run->local_search_only_terminal = DEFAULT_LOCAL_SEARCH_ONLY_TERMINAL;
    run->local_search_remove_movement = DEFAULT_LOCAL_SEARCH_REMOVE_MOVEMENT;

    // Facility distances not yet computed
    for(int mode=0;mode<N_FACDIS_MODES;mode++){
        run->facs_distance[mode] = NULL;
    }

    // First restart data
    run->firstr_per_size_n_sols = safe_malloc(sizeof(int)*(prob->n_facs+2));
    memset(run->firstr_per_size_n_sols,0,     sizeof(int)*(prob->n_facs+2));
    run->firstr_per_size_n_sols_after_red = safe_malloc(sizeof(int)*(prob->n_facs+2));
    memset(run->firstr_per_size_n_sols_after_red,0,     sizeof(int)*(prob->n_facs+2));
    run->firstr_per_size_n_local_optima = safe_malloc(sizeof(int)*(prob->n_facs+2));
    memset(run->firstr_per_size_n_local_optima,0,     sizeof(int)*(prob->n_facs+2));

    // Restart data
    run->n_restarts = n_restarts;
    run->restart_times  = safe_malloc(sizeof(double)*run->n_restarts);
    run->restart_values = safe_malloc(sizeof(double)*run->n_restarts);
    for(int r=0;r<run->n_restarts;r++){
        run->restart_times[r] = NAN;
        run->restart_values[r] = -INFINITY;
    }

    // Perform precomputations
    rundata_precompute(run,rstrats,n_rstrats);

    return run;
}

// Prints a briefing of the rundata parameters
void rundata_print(const rundata *run, FILE *fp){
    const problem *prob = run->prob;

    fprintf(fp,"== PROBLEM ==\n");
    fprintf(fp,"# N_FACILITIES: %d\n",prob->n_facs);
    fprintf(fp,"# N_CLIENTS: %d\n",prob->n_clis);
    fprintf(fp,"# TRANSPORT_COST: %lf\n",prob->transport_cost);
    fprintf(fp,"# CLIENT_GAIN: %lf\n",prob->client_gain);
    fprintf(fp,"# UNASSIGNED_COST: %lf\n",prob->unassigned_cost);
    fprintf(fp,"# SIZE_RESTRICTION_MINIMUM: %d\n",prob->size_restriction_minimum);
    fprintf(fp,"# SIZE_RESTRICTION_MAXIMUM: %d\n",prob->size_restriction_maximum);
    fprintf(fp,"# VERBOSE: %d\n",prob->verbose);
    fprintf(fp,"\n");
    fprintf(fp,"== RUN DATA ==\n");
    fprintf(fp,"# FILTER: %s (%d)\n",filter_names[run->filter],(int)run->filter);
    fprintf(fp,"# BRANCH_AND_BOUND: %d\n",run->branch_and_bound);
    fprintf(fp,"# LOWER_BOUND: %lf\n",run->lower_bound);
    fprintf(fp,"# TARGET_SOLS: %d\n",run->target_sols);
    fprintf(fp,"# PRECOMP_CLIENT_OPTIMAL_GAIN: %lf\n",run->precomp_client_optimal_gain);
    fprintf(fp,"# PRECOMP_EMPTY_VALUE: %lf\n",run->precomp_empty_value);
    fprintf(fp,"# N_THREADS: %d\n",run->n_threads);
    fprintf(fp,"# LOCAL_SEARCH: %s\n",local_search_names[run->local_search]);
    fprintf(fp,"# LOCAL_SEARCH_ONLY_TERMINAL: %d\n",run->local_search_only_terminal);
    fprintf(fp,"# LOCAL_SEARCH_REMOVE_MOVEMENT: %d\n",run->local_search_remove_movement);
    fprintf(fp,"# RANDOM_SEED: %d\n",run->random_seed);
    fprintf(fp,"# RESTARTS: %d\n",run->n_restarts);
}
