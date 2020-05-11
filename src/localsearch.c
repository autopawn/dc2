#include "localsearch.h"

#define NO_MOVEMENT (-2)

int solution_whitaker_hill_climbing(const rundata *run, solution *sol, shuffler *shuff){
    const problem *prob = run->prob;

    // Is this first improvement?
    int first_improvement = run->local_search==SWAP_FIRST_IMPROVEMENT;
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
        int best_rem = NO_MOVEMENT;
        int best_ins = NO_MOVEMENT;

        // Reshuffle shuffler
        if(first_improvement) shuffler_reshuffle(shuff);

        // Check movements that are allowed and not allowed
        int allow_size_decrease = run->local_search_rem_movement && (prob->size_restriction_minimum==-1 || sol->n_facs>prob->size_restriction_minimum);
        int allow_size_increase = run->local_search_add_movement && (prob->size_restriction_maximum==-1 || sol->n_facs<prob->size_restriction_maximum);

        // Consider adding a facility while removing the worst facility
        int k_ini = allow_size_decrease? -1 : 0; // also consider not adding a facility if allow_size_decrease
        for(int k=k_ini;k<prob->n_facs;k++){
            int f_ins;
            if(first_improvement && k!=-1){
                f_ins = shuffler_next(shuff);
            }else{
                f_ins = k;
            }

            // Ignore insertion if the facility is already present in the solution
            if(f_ins>=0 && used[f_ins]) continue;
            // Find best facility to remove after inserting f_ins, and profits
            int f_rem;
            double delta_profit, delta_profit_worem;
            solution_findout(prob,sol,f_ins,v,phi2,&f_rem,&delta_profit,&delta_profit_worem);
            // Update best removal and insertion
            int improvement = 0;
            if(delta_profit>best_delta){
                best_delta = delta_profit;
                best_rem = f_rem;
                best_ins = f_ins;
                improvement = 1;
            }
            if(delta_profit_worem>best_delta && allow_size_increase){
                best_delta = delta_profit_worem;
                best_rem = -1;
                best_ins = f_ins;
                improvement = 1;
            }
            // break the loop on first_improvement
            if(first_improvement && improvement) break;
        }
        // No inserting and no removing is equivalent to not performing a movement
        if(best_ins==-1 && best_rem==-1){
            best_ins = NO_MOVEMENT;
            best_rem = NO_MOVEMENT;
        }

        // Stop when no movement results in a better solution:
        if(best_ins==NO_MOVEMENT) break;

        // Perform swap:
        assert(best_rem!=NO_MOVEMENT);
        double old_value = sol->value;
        if(best_rem!=-1){
            solution_remove(prob,sol,best_rem,phi2);
            used[best_rem] = 0;
        }
        if(best_ins!=-1){
            solution_add(prob,sol,best_ins);
            used[best_ins] = 1;
        }
        assert(sol->value>old_value);
        #ifdef DEBUG
            double error_pred = (sol->value-old_value)-best_delta;
            if(error_pred<0) error_pred *= -1;
            assert(error_pred<1e-5 || isnan(error_pred));
        #endif

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
    for(int r=args->thread_id;r<args->n_sols;r+=args->run->n_threads){
        // Perform local search on the given solution
        args->n_moves += solution_whitaker_hill_climbing(args->run,args->sols[r],args->shuff);
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
