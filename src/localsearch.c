#include "localsearch.h"

int solution_whitaker_hill_climbing(const problem *prob, solution *sol, shuffler *shuff){
    if(sol->n_facs==0) return 0;
    // Is this first improvement?
    int first_improvement = prob->local_search==SWAP_FIRST_IMPROVEMENT;
    assert(shuff!=NULL || !first_improvement);
    // First and Second nearest facility to each client
    int *phi1 = safe_malloc(sizeof(int)*prob->n_clis);
    int *phi2 = safe_malloc(sizeof(int)*prob->n_clis);
    // Array if a facility appears in the solution
    int *used = safe_malloc(sizeof(int)*prob->n_facs);
    // Auxiliar array for solution_findout
    double *v = safe_malloc(sizeof(double)*prob->n_facs);
    // Initialize arrays
    for(int i=0;i<prob->n_facs;i++){
        used[i] = 0;
        v[i] = -INFINITY;
    }
    for(int i=0;i<sol->n_facs;i++){
        used[sol->facs[i]] = 1;
    }
    // For each client, find the second nearest facility in the solution
    for(int i=0;i<prob->n_clis;i++){
        phi1[i] = sol->assigns[i];
        phi2[i] = solution_client_2nd_nearest(prob,sol,i);
    }
    // Each movement:
    int n_moves = 0;
    while(1){
        // Insertion candidate:
        double best_delta = 0;
        int best_rem = -1;
        int best_ins = -1;

        if(first_improvement) shuffler_reshuffle(shuff);

        for(int k=0;k<prob->n_facs;k++){
            int f_ins = first_improvement? shuffler_next(shuff) : k;
            if(used[f_ins]) continue;
            // Find the best option for removal:
            int f_rem;
            double delta_profit;
            solution_findout(prob,sol,f_ins,v,phi2,&f_rem,&delta_profit);
            if(delta_profit>best_delta){
                best_delta = delta_profit;
                best_rem = f_rem;
                best_ins = f_ins;
                // If it's a first_improvement policy, break the loop now.
                if(first_improvement) break;
            }
        }
        // Stop when no movement results in a better solution:
        if(best_ins==-1) break;
        // Perform swap:
        assert(best_rem!=-1);
        double old_value = sol->value;
        solution_remove(prob,sol,best_rem,phi2);
        used[best_rem] = 0;
        solution_add(prob,sol,best_ins);
        used[best_ins] = 1;
        assert(sol->value>old_value);

        // Update phi1 and phi2
        update_phi1_and_phi2(prob,sol,best_ins,best_rem,phi1,phi2);

        // Count one move:
        n_moves += 1;
    }

    // Free memory
    free(v);
    free(used);
    free(phi2);
    free(phi1);
    return n_moves;
}

void update_phi1_and_phi2(const problem *prob, const solution *sol, int f_ins, int f_rem,
        int *phi1, int *phi2){
    for(int i=0;i<prob->n_clis;i++){
        // Update phi2
        int near[3];
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
        assert(phi1[i]!=phi2[i]);
    }
}

void solutions_delete_repeated(const problem *prob, solution **sols, int *n_sols){
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
    const problem *prob;
    solution **sols;
    int n_sols;
    int n_moves;
    shuffler *shuff;
} hillclimb_thread_args;

void *hillclimb_thread_execution(void *arg){
    hillclimb_thread_args *args = (hillclimb_thread_args *) arg;
    for(int r=args->thread_id;r<args->n_sols;r+=args->prob->n_threads){
        // Perform local search on the given solution
        args->n_moves += solution_whitaker_hill_climbing(args->prob,args->sols[r],args->shuff);
    }
    return NULL;
}

// Perform local searches (in parallel).
void solutions_hill_climbing(problem *prob, solution **sols, int n_sols){
    // Start measuring time
    clock_t start = clock();
    // Allocate memory for threads and arguments
    pthread_t *threads = safe_malloc(sizeof(pthread_t)*prob->n_threads);
    hillclimb_thread_args *targs = safe_malloc(sizeof(hillclimb_thread_args)*prob->n_threads);
    // Call all threads to perform local search
    for(int i=0;i<prob->n_threads;i++){
        // Set arguments for the thread
        targs[i].thread_id = i;
        targs[i].prob = prob;
        targs[i].sols = sols;
        targs[i].n_sols = n_sols;
        targs[i].n_moves = 0;
        // Set random number generator for the thread
        if(prob->local_search==SWAP_FIRST_IMPROVEMENT){
            targs[i].shuff = shuffler_init(prob->n_facs);
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
    for(int i=0;i<prob->n_threads;i++){ // Join threads
        pthread_join(threads[i],NULL);
        n_moves += targs[i].n_moves;
    }
    // Free memory
    for(int i=0;i<prob->n_threads;i++){
        if(targs[i].shuff!=NULL) shuffler_free(targs[i].shuff);
    }
    free(targs);
    free(threads);
    // End measuring time
    clock_t end = clock();
    double seconds = (double)(end - start) / (double)CLOCKS_PER_SEC;
    prob->n_local_searches += n_sols;
    prob->n_local_search_movements += n_moves;
    prob->local_search_seconds += seconds;
}