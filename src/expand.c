#include "expand.h"

//#############################################################
// FUTURESOLS
//#############################################################

// Possible future solution that results from another one.
typedef struct {
    const solution *origin;
    int newf;
    uint hash;
    int n_facs;
    int facs[0]; // Flexible array member.
} futuresol;

// Compares futuresols so that similar are consecutive.
int futuresol_cmp(const void *a, const void *b){
    const futuresol *aa = (const futuresol *) a;
    const futuresol *bb = (const futuresol *) b;
    if(aa->hash>bb->hash) return +1;
    if(aa->hash<bb->hash) return -1;
    int nf_delta = aa->n_facs - bb->n_facs;
    if(nf_delta!=0) return nf_delta;
    for(int i=0;i<aa->n_facs;i++){
        int idx_delta = aa->facs[i]-bb->facs[i];
        if(idx_delta!=0) return idx_delta;
    }
    return 0;
}

// Inits a futuresol from a current sol and a new facility
int futuresol_init_from(futuresol *fsol, const solution *sol, int newf){
    fsol->origin = sol;
    fsol->newf = newf;
    fsol->n_facs = sol->n_facs;
    // Copy facilitites, check if f already exists, init hash.
    fsol -> hash = hash_int(newf);
    for(int k=0;k<sol->n_facs;k++){
        if(sol->facs[k]==newf) return 0; // If the solution candidate is not new
        fsol->facs[k] = sol->facs[k];
        // Include solution on the hash
        fsol->hash = fsol->hash ^ hash_int(fsol->facs[k]);
    }
    // Add the new facility to the inner array
    add_to_sorted(fsol->facs,&fsol->n_facs,newf);
    return 1;
}

//#############################################################
// GENERATION OF NEW SOLUTIONS FROM FUTURESOLS
//#############################################################

typedef struct {
    int thread_id;
    int n_threads;
    const problem *prob;
    int n_fsols;
    void *futuresols;
    size_t fsol_size;
    solution **out_sols;
} expand_thread_args;

void *expand_thread_execution(void *arg){
    expand_thread_args *args = (expand_thread_args *) arg;
    for(int r=args->thread_id;r<args->n_fsols;r+=args->n_threads){
        // Generate a new solution from the fsol, and then check if it passes filtering.
        futuresol *fsol = (futuresol *)(args->futuresols+args->fsol_size*r);
        solution *new_sol = solution_copy(args->prob,fsol->origin);
        solution_add(args->prob,new_sol,fsol->newf);
        if(args->prob->filter >= BETTER_THAN_SUBSETS){
            // TODO: implement.
        }else if(args->prob->filter > NO_FILTER){
            if(new_sol->value > fsol->origin->value){
                args->out_sols[r] = new_sol;
            }else{
                solution_free(new_sol);
                args->out_sols[r] = NULL;
            }
        }
    }
    return NULL;
}

//#############################################################
// EXPANSION
//#############################################################

solution **new_expand_solutions(const problem *prob,
        solution **sols, int n_sols, int *out_n_sols){
    // TODO: Remove this once the missing filter is added
    assert(prob->filter<BETTER_THAN_SUBSETS);
    // Get the corrent size of the solutions on this expansion:
    int current_size = n_sols>0? sols[0]->n_facs : 0;
    // Compute the size of each futuresol (flexible array member must be added):
    size_t fsol_size = sizeof(futuresol)+sizeof(int)*(current_size+1);
    // Allocate enough memory for the maximium amount of futuresols that can appear:   
    void *futuresols = safe_malloc(fsol_size*n_sols*prob->n_facs);
    // Get the candidates to future solutions
    int n_futuresols = 0;
    for(int i=0;i<n_sols;i++){
        assert(sols[i]->n_facs==current_size); // All solutions are expected to have the same size.
        for(short f=0;f<prob->n_facs;f++){
            futuresol *fsol = (futuresol *)(futuresols+fsol_size*n_futuresols);
            n_futuresols += futuresol_init_from(fsol,sols[i],f);
        }
    }
    
    { // Sort the futuresols in order to detect the similar ones faster:
        qsort(futuresols,n_futuresols,fsol_size,futuresol_cmp);
        int n_futuresols2 = 0;
        if(n_futuresols>0){
            futuresol *last_fsol = (futuresol *)(futuresols);
            for(int r=1;r<n_futuresols;r++){
                futuresol *fsol = (futuresol *)(futuresols+fsol_size*r);
                // Compare fsol with the last_fsol:
                int ftsol_cmp = futuresol_cmp(last_fsol,fsol);
                // Check if fsol creates a brave new solution.
                if(ftsol_cmp!=0){
                    futuresol *next_pos = (futuresol *)(futuresols+fsol_size*n_futuresols2);
                    memmove(next_pos,fsol,fsol_size);
                    last_fsol = next_pos;
                    n_futuresols2 += 1;
                }else{
                    /* If fsol doesn't create a new solution but creates it from a better (worst) one, in that case fsol replaces last_fsol. 
                    Because, depending on the filter, the new solution should be better that the better (worst) one that generates it. */
                    int is_better = fsol->origin->value>last_fsol->origin->value;
                    if((prob->filter>=BETTER_THAN_ALL_PARENTS) == is_better){
                        memmove(last_fsol,fsol,fsol_size);
                    }
                }
            }
        }
        // Update the futuresols, and realloc to reduce memory usage
        n_futuresols = n_futuresols2;
        futuresols = realloc(futuresols,fsol_size*n_futuresols);
    }

    solution **out_sols = safe_malloc(sizeof(solution*)*n_futuresols);
    { // Create new solutions [in parallel]
        #if THREADS>0
            const int n_threads = THREADS;
        #else
            const int n_threads = 1;
        #endif
        expand_thread_args *targs = safe_malloc(sizeof(expand_thread_args)*n_threads);
        for(int i=0;i<n_threads;i++){
            targs[i].thread_id = i;
            targs[i].n_threads = n_threads;
            targs[i].prob = prob;
            targs[i].n_fsols = n_futuresols;
            targs[i].futuresols = futuresols;
            targs[i].fsol_size = fsol_size;
            targs[i].out_sols = out_sols;
        }
        #if THREADS>0
            // Generate threads in order to expand the solutions
            pthread_t *threads = safe_malloc(sizeof(pthread_t)*n_threads);
            for(int i=0;i<n_threads;i++){
                int rc = pthread_create(&threads[i],NULL,expand_thread_execution,&targs[i]);
                if(rc){
                    printf("Error %d on thread creation\n",rc);
                    exit(1);
                }
            }
            // Join threads
            for(int i=0;i<THREADS;i++){
                pthread_join(threads[i],NULL);
            }
            //
            free(threads);
        #else
            expand_thread_execution(&targs[0]);
        #endif
        free(targs);
    }

    { // Eliminate NULLed out solutions
        int n_sols = 0;
        for(int r=0;r<n_futuresols;r++){
            if(out_sols[r]!=NULL){
                out_sols[n_sols] = out_sols[r];
                n_sols += 1;
            }
        }
        out_sols = realloc(out_sols,sizeof(solution*)*(n_sols));
        *out_n_sols = n_sols;
    }

    free(futuresols);
    return out_sols;
}