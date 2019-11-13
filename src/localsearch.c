#include "localsearch.h"

void solution_whitaker_hill_climbing(const problem *prob, solution *sol){
    if(sol->n_facs==0) return;
    // Second nearest facility to each client
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
    // Each movement:
    while(1){
        // Clients nearest to the solution
        for(int i=0;i<prob->n_clis;i++){
            phi2[i] = solution_client_2nd_nearest(prob,sol,i);
        }
        // Insertion candidate:
        double best_delta = 0;
        int best_rem = -1;
        int best_ins = -1;
        for(int f=0;f<prob->n_facs;f++){
            if(used[f]) continue;
            // Find the best option for removal:
            int f_rem;
            double delta_profit;
            solution_findout(prob,sol,f,v,phi2,&f_rem,&delta_profit);
            if(delta_profit>best_delta){
                best_delta = delta_profit;
                best_rem = f_rem;
                best_ins = f;
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
    }

    // Free memory
    free(v);
    free(used);
    free(phi2);
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
} hillclimb_thread_args;

void *hillclimb_thread_execution(void *arg){
    hillclimb_thread_args *args = (hillclimb_thread_args *) arg;
    for(int r=args->thread_id;r<args->n_sols;r+=args->prob->n_threads){
        // Perform local search on the given solution
        solution_whitaker_hill_climbing(args->prob,args->sols[r]);
    }
    return NULL;
}

// Perform local searches (in parallel).
void solutions_hill_climbing(const problem *prob, solution **sols, int n_sols){
    // Allocate memory for threads and arguments
    pthread_t *threads = safe_malloc(sizeof(pthread_t)*prob->n_threads);
    hillclimb_thread_args *targs = safe_malloc(sizeof(hillclimb_thread_args)*prob->n_threads);
    // Call all threads to perform local search
    for(int i=0;i<prob->n_threads;i++){
        targs[i].thread_id = i;
        targs[i].prob = prob;
        targs[i].sols = sols;
        targs[i].n_sols = n_sols;
        int rc = pthread_create(&threads[i],NULL,hillclimb_thread_execution,&targs[i]);
        if(rc){
            fprintf(stderr,"ERROR: Error %d on pthread_create\n",rc);
            exit(1);
        }
    }
    // Join threads
    for(int i=0;i<prob->n_threads;i++){ // Join threads
        pthread_join(threads[i],NULL);
    }
    // Free memory
    free(targs);
    free(threads);
}