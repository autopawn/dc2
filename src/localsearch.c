#include "localsearch.h"

void update_phi1_and_phi2(const problem *prob, const solution *sol, int f_ins, int f_rem,
        int *phi1, int *phi2, int *affected_mask){
    for(int i=0;i<prob->n_clis;i++){
        int before_phi1 = phi1[i];
        int before_phi2 = phi2[i];
        // Update phi2
        int near[3];
        if(f_ins==-1){
            near[0] = phi1[i];
            near[1] = phi2[i];
            near[2] = -2; // Unknown, may be f_ins or a phi3[i].
        }else{
            if(sol->assigns[i]==f_ins){
                near[0] = f_ins;
                near[1] = phi1[i];
                near[2] = phi2[i];
            }else{
                double phi2_val = problem_assig_value(prob,phi2[i],i);
                double best_ins_val = problem_assig_value(prob,f_ins,i);
                if(best_ins_val>phi2_val){
                    near[0] = phi1[i];
                    near[1] = f_ins;
                    near[2] = phi2[i];
                }else{
                    near[0] = phi1[i];
                    near[1] = phi2[i];
                    near[2] = -2; // Unknown, may be f_ins or a phi3[i].
                }
            }
        }
        // Find phi2
        int k = 0;
        for(int u=0;u<3;u++){
            if(near[u]==f_rem) continue;
            k += 1;
            if(k==2) phi2[i] = near[u];
        }
        if(phi2[i]==-2) phi2[i] = solution_client_2nd_nearest(prob,sol,i);

        // Update phi1
        phi1[i] = sol->assigns[i];
        assert(phi1[i]!=phi2[i] || phi1[i]==-1);
        assert(problem_assig_value(prob,phi1[i],i) >= problem_assig_value(prob,phi2[i],i));
        assert(affected_mask==NULL || affected_mask[i] || (before_phi1==phi1[i] && before_phi2==phi2[i]));
    }
}

void solutions_delete_repeated(solution **sols, int *n_sols){
    // If there are 0 solutions, do nothing.
    if(*n_sols==0) return;
    // Sort solutions for easier comparison
    qsort(sols,*n_sols,sizeof(solution *),solutionp_facs_cmp);
    // Delete repeated solutions
    int n_final = 1;
    for(int i=1;i<*n_sols;i++){
        if(solutionp_facs_cmp(&sols[n_final-1],&sols[i])!=0){
            sols[n_final] = sols[i];
            n_final++;
        }else{
            solution_free(sols[i]);
        }
    }
    *n_sols = n_final;
}

typedef struct {
    int thread_id;
    const rundata *run;
    solution **sols;
    int n_sols;
    int n_moves;
    shuffler *shuff;
} hillclimb_thread_args;

void *hillclimb_thread_execution(void *arg){
    hillclimb_thread_args *args = (hillclimb_thread_args *) arg;
    if(args->run->local_search==SWAP_RESENDE_WERNECK){
        if(args->thread_id<args->n_sols){
            fastmat *mat = fastmat_init(args->run->prob->n_facs,args->run->prob->n_facs);
            for(int r=args->thread_id;r<args->n_sols;r+=args->run->n_threads){
                // Perform local search on the given solution
                args->n_moves += solution_resendewerneck_hill_climbing(args->run,args->sols[r],mat);
            }
            fastmat_free(mat);
        }
    }else{
        for(int r=args->thread_id;r<args->n_sols;r+=args->run->n_threads){
            // Perform local search on the given solution
            args->n_moves += solution_whitaker_hill_climbing(args->run,args->sols[r],args->shuff);
        }
    }
    return NULL;
}

// Perform local searches (in parallel).
void solutions_hill_climbing(rundata *run, solution **sols, int n_sols){
    // Start measuring time
    clock_t start = clock();
    // Allocate memory for threads and arguments
    pthread_t *threads = safe_malloc(sizeof(pthread_t)*run->n_threads);
    hillclimb_thread_args *targs = safe_malloc(sizeof(hillclimb_thread_args)*run->n_threads);
    // Call all threads to perform local search
    for(int i=0;i<run->n_threads;i++){
        // Set arguments for the thread
        targs[i].thread_id = i;
        targs[i].run = run;
        targs[i].sols = sols;
        targs[i].n_sols = n_sols;
        targs[i].n_moves = 0;
        // Set random number generator for the thread
        if(run->local_search==SWAP_FIRST_IMPROVEMENT){
            targs[i].shuff = shuffler_init(run->prob->n_facs);
        }else{
            targs[i].shuff = NULL;
        }
        // Create thread
        int rc = pthread_create(&threads[i],NULL,hillclimb_thread_execution,&targs[i]);
        if(rc){
            fprintf(stderr,"ERROR: Error %d on pthread_create\n",rc);
            exit(1);
        }
    }
    // Join threads
    int n_moves = 0;
    for(int i=0;i<run->n_threads;i++){ // Join threads
        pthread_join(threads[i],NULL);
        n_moves += targs[i].n_moves;
    }
    // Free memory
    for(int i=0;i<run->n_threads;i++){
        if(targs[i].shuff!=NULL) shuffler_free(targs[i].shuff);
    }
    free(targs);
    free(threads);
    // End measuring time
    clock_t end = clock();
    double seconds = (double)(end - start) / (double)CLOCKS_PER_SEC;
    run->n_local_searches += n_sols;
    run->n_local_search_movements += n_moves;
    run->local_search_seconds += seconds;
}
