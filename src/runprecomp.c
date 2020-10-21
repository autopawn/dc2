#include <assert.h>

#include "runprecomp.h"

// ============================================================================
// Facility-facility distance precomputation thread execution

typedef struct {
    runprecomp *pcomp;
    const problem *prob;
    int thread_id;
    int n_threads;
    int mode;
} precomp_facs_dist_thread_args;

void *precomp_facs_dist_thread_execution(void *arg){
    precomp_facs_dist_thread_args *args = (precomp_facs_dist_thread_args *) arg;
    const problem *prob = args->prob;
    // Compute facility-facility distances acording to mode
    for(int a=args->thread_id;a<prob->n_facs;a+=args->n_threads){
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
                    double dist_sum = prob->distance_cost[a][j]+prob->distance_cost[b][j];
                    if(dist_sum<dist) dist = dist_sum;
                }
            }else{
                assert(!"Asked to precompute valid facility distance.");
            }
            args->pcomp->facs_distance[args->mode][a][b] = dist;
            args->pcomp->facs_distance[args->mode][b][a] = dist;
        }
    }
    return NULL;
}

// ============================================================================
// Nearly indexes thread execution

typedef struct {
    runprecomp *pcomp;
    const problem *prob;
    int thread_id;
    int n_threads;
} precomp_nearly_indexes_args;

typedef struct {
    double value;
    int indx;
} distpair;

// Compare distpairs by distance
int distpair_cmp(const void *a,const void *b){
    const distpair *aa = a;
    const distpair *bb = b;
    double diff = bb->value - aa->value;
    assert(!isnan(diff));
    if(diff<0) return -1;
    if(diff==0) return 0;
    return 1;
}

void *precomp_nearly_indexes_thread_execution(void *arg){
    precomp_nearly_indexes_args *args = (precomp_nearly_indexes_args *) arg;
    const problem *prob = args->prob;
    runprecomp *pcomp = args->pcomp;
    // Compute facility indexes sorted by proximity, client-wise.
    for(int i=args->thread_id;i<prob->n_clis;i+=args->n_threads){
        // Initialize array of distpairs with distances and facility indexes
        distpair *pairs = safe_malloc(sizeof(distpair)*prob->n_facs);
        for(int f=0;f<prob->n_facs;f++){
            pairs[f].value = problem_assig_value(prob,f,i);
            pairs[f].indx = f;
        }
        // Sort pairs by distance
        qsort(pairs,prob->n_facs,sizeof(distpair),distpair_cmp);
        assert(prob->n_facs==0 || pairs[0].value >= pairs[prob->n_facs-1].value);
        // Initialize nearly indexes list
        for(int k=0;k<prob->n_facs;k++){
            pcomp->nearly_indexes[i][k] = pairs[k].indx;
        }
        assert(prob->n_facs==0 || problem_assig_value(prob,pcomp->nearly_indexes[i][0],i) >= problem_assig_value(prob,pcomp->nearly_indexes[i][prob->n_facs-1],i));
        // Free pairs data structure
        free(pairs);
    }
    return NULL;
}

// ============================================================================

runprecomp *runprecomp_init(const problem *prob, redstrategy *rstrats, int n_rstrats, int precomp_nearly_indexes, int n_threads, int verbose){

    runprecomp *pcomp = safe_malloc(sizeof(runprecomp));
    // Facility distances not yet computed
    for(int mode=0;mode<N_FACDIS_MODES;mode++){
        pcomp->facs_distance[mode] = NULL;
    }
    // Nearly indexes for each client not yet computed
    pcomp->nearly_indexes = NULL;

    pcomp->n_clis = prob->n_clis;
    pcomp->n_facs = prob->n_facs;

    // Precompute precomp_client_optimal_gain
    pcomp->precomp_client_optimal_gain = 0;
    for(int k=0;k<prob->n_clis;k++){
        double best_val = problem_assig_value(prob,-1,k);
        for(int i=0;i<prob->n_facs;i++){
            double val = problem_assig_value(prob,i,k);
            if(best_val < val) best_val = val;
        }
        pcomp->precomp_client_optimal_gain += best_val;
    }

    // Precompute facility-facility distances
    int notification = 0;
    for(int r=0;r<n_rstrats;r++){
        int mode = redstrategy_required_facdis_mode(rstrats[r]);
        // Compute distance matrix for the mode required by the reduction strategy (if hasn't been computed already)
        if(mode!=FACDIS_NONE && pcomp->facs_distance[mode]==NULL){
            if(!notification && verbose!=0){
                printf("\nPrecomputing facility-facility distances.\n");
                notification = 1;
            }
            // Allocate distance matrix between facilities
            pcomp->facs_distance[mode] = safe_malloc(sizeof(double*)*prob->n_facs);
            for(int i=0;i<prob->n_facs;i++){
                pcomp->facs_distance[mode][i] = safe_malloc(sizeof(double)*prob->n_facs);
            }
            // Allocate memory for threads and arguments
            pthread_t *threads = safe_malloc(sizeof(pthread_t)*n_threads);
            precomp_facs_dist_thread_args *targs = safe_malloc(sizeof(precomp_facs_dist_thread_args)*n_threads);
            // Call threads to compute facility-facility distances
            for(int i=0;i<n_threads;i++){
                targs[i].pcomp = pcomp;
                targs[i].prob  = prob;
                targs[i].thread_id = i;
                targs[i].n_threads = n_threads;
                targs[i].mode = mode;

                int rc = pthread_create(&threads[i],NULL,precomp_facs_dist_thread_execution,&targs[i]);
                if(rc){
                    fprintf(stderr,"ERROR: Error %d on pthread_create\n",rc);
                    exit(1);
                }
            }
            // Join threads
            for(int i=0;i<n_threads;i++){ // Join threads
                pthread_join(threads[i],NULL);
            }
            // Free memory
            free(targs);
            free(threads);
        }
    }

    // Precompute facility indexes by proximity to each client
    if(precomp_nearly_indexes){
        if(pcomp->nearly_indexes==NULL){
            if(verbose!=0){
                printf("\nPrecomputing nearly indexes.\n");
            }
            // Initialize memory for the nearly_indexes
            pcomp->nearly_indexes = safe_malloc(sizeof(int*)*prob->n_clis);
            for(int i=0;i<prob->n_clis;i++){
                pcomp->nearly_indexes[i] = safe_malloc(sizeof(int)*prob->n_facs);
            }
            // Allocate memory for threads and arguments
            pthread_t *threads = safe_malloc(sizeof(pthread_t)*n_threads);
            precomp_nearly_indexes_args *targs = safe_malloc(sizeof(precomp_nearly_indexes_args)*n_threads);
            // Call threads to compute facility-facility distances
            for(int i=0;i<n_threads;i++){
                targs[i].pcomp = pcomp;
                targs[i].prob  = prob;
                targs[i].thread_id = i;
                targs[i].n_threads = n_threads;
                int rc = pthread_create(&threads[i],NULL,precomp_nearly_indexes_thread_execution,&targs[i]);
                if(rc){
                    fprintf(stderr,"ERROR: Error %d on pthread_create\n",rc);
                    exit(1);
                }
            }
            // Join threads
            for(int i=0;i<n_threads;i++){ // Join threads
                pthread_join(threads[i],NULL);
            }
            // Free memory
            free(targs);
            free(threads);
        }
    }

    return pcomp;
}

void runprecomp_free(runprecomp *pcomp){
    // Free facility-facility distance matrices, if they are in use:
    for(int mode=0;mode<N_FACDIS_MODES;mode++){
        if(pcomp->facs_distance[mode]){
            for(int i=0;i<pcomp->n_facs;i++){
                free(pcomp->facs_distance[mode][i]);
            }
            free(pcomp->facs_distance[mode]);
        }
    }
    // Free nearly indexes if it is in use
    if(pcomp->nearly_indexes){
        for(int i=0;i<pcomp->n_clis;i++){
            free(pcomp->nearly_indexes[i]);
        }
        free(pcomp->nearly_indexes);
    }
    // Free the precomputation
    free(pcomp);
}
