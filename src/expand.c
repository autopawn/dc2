#include "expand.h"

//#############################################################
// FUTURESOLS
//#############################################################

// Possible future solution that results from another one.
typedef struct {
    solution *origin;
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
int futuresol_init_from(futuresol *fsol, solution *sol, int newf){
    assert(sol!=NULL);
    fsol->origin = sol;
    fsol->newf = newf;
    fsol->n_facs = sol->n_facs;
    // Copy facilitites, check if f already exists, init hash.
    fsol -> hash = hash_int(newf);
    for(int k=0;k<sol->n_facs;k++){
        if(sol->facs[k]==newf) return 0; // If the solution candidate is not new // TODO: make faster when sols->nfacs is near n
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
    const rundata *run;
    int n_fsols;
    void *futuresols;
    size_t fsol_size;
    solution **out_sols;
} expand_thread_args;

// Check if a solution passes the value of which it would be filtered.
// Equal valued solutions are not filtered if the size restriction has not yet ben reached
// Or if the other value is -INFINITY.
int is_filtered(const problem *prob, const solution *sol, double other){
    int not_equality = sol->n_facs<=prob->size_restriction_minimum || other<=-INFINITY;
    if(not_equality){
        return sol->value < other;
    }else{
        return sol->value <= other;
    }
}

void *expand_thread_execution(void *arg){
    expand_thread_args *args = (expand_thread_args *) arg;
    const problem *prob = args->run->prob;

    // Auxiliary arrays that could be useful
    int *phi2 = NULL;
    double *v = NULL;

    for(int r=args->thread_id;r<args->n_fsols;r+=args->run->n_threads){
        // Generate a new solution from the fsol, and then check if it passes filtering.
        futuresol *fsol = (futuresol *)(args->futuresols+args->fsol_size*r);
        solution *new_sol = solution_copy(prob,fsol->origin);
        solution_add(prob,new_sol,fsol->newf,NULL);
        int filtered = 0;
        // Must be better than any other subset (minus 1 facility)
        if(args->run->filter >= BETTER_THAN_SUBSETS){
            // Initialize useful arrays if they aren't already
            if(v==NULL){
                v = safe_malloc(sizeof(double)*prob->n_facs);
                for(int i=0;i<prob->n_facs;i++) v[i] = -INFINITY;
            }
            if(phi2==NULL) phi2 = safe_malloc(sizeof(int)*prob->n_clis);
            // Find second nearest facility for each client
            for(int i=0;i<prob->n_clis;i++){
                phi2[i] = solution_client_2nd_nearest(prob,new_sol,i);
            }
            // Check if there's profit after picking the best facility for removal
            int f_rem;
            double delta_profit,delta_profit_worem;
            solution_findout(prob,new_sol,-1,v,phi2,NULL,&f_rem,&delta_profit,&delta_profit_worem);
            filtered = is_filtered(prob,new_sol,new_sol->value+delta_profit);
        }
        // Filters that use fsol->origin
        // Notice that if the filter is BETTER_THAN_ONE_PARENT, then fsol->origin is the worst parent
        // If it is BETTER_THAN_ALL_PARENTS, then fsol->origin is the best parent
        else if(args->run->filter >= BETTER_THAN_ONE_PARENT){
            filtered = is_filtered(prob,new_sol,fsol->origin->value);
        }
        // If it must be better than the empty solution
        else if(args->run->filter == BETTER_THAN_EMPTY){
            filtered = is_filtered(prob,new_sol,-INFINITY);
        }
        // Delete the solution
        if(filtered){
            solution_free(new_sol);
            args->out_sols[r] = NULL;
        }else{
            args->out_sols[r] = new_sol;
        }
    }
    // Free auxilary arrays if they were allocated
    if(v!=NULL) free(v);
    if(phi2!=NULL) free(phi2);
    //
    return NULL;
}

//#############################################################
// EXPANSION
//#############################################################

solution **new_expand_solutions(const rundata *run,
        solution **sols, int n_sols, int *out_n_sols, int pool_size){
    const problem *prob = run->prob;
    // Get the corrent size of the solutions on this expansion:
    int current_size = n_sols>0? sols[0]->n_facs : 0;
    // Compute the size of each futuresol (flexible array member must be added):
    size_t fsol_size = sizeof(futuresol)+sizeof(int)*(current_size+1);

    // ==== Generate futuresols depending on the branching factor
    assert(run->branching_factor>=-1);

    float branchingf;
    if(run->branching_factor==-1){
        // All children are generated
        branchingf = prob->n_facs-current_size;
    }else if(run->branching_factor==0){
        // Generate according to Resende & Werneck's formula
        branchingf = ceil(log2f((float)prob->n_facs/fmaxf(1.0,current_size)));
    }else{
        // Use the branching factor given
        branchingf = run->branching_factor;
    }
    // Branching factor correction ensures that enough solutions are created at the first generations
    if(run->branching_correction){
        float expected_sols = 1;
        for(int i=0;i<current_size;i++){
            expected_sols *= prob->n_facs;
            if(expected_sols>pool_size) break;
        }
        if(expected_sols>pool_size) expected_sols = pool_size;
        branchingf *= pool_size/expected_sols;
    }
    // Final branching factor
    int branching = ceilf(branchingf);
    if(branching>(prob->n_facs-current_size)) branching = prob->n_facs-current_size;

    // Allocate enough memory for the maximium amount of futuresols that can appear:
    void *futuresols = safe_malloc(fsol_size*(n_sols*branching+1));
    int n_futuresols = 0;

    // Get the candidates to future solutions
    if(branching >= prob->n_facs-current_size){
        // Full branching, all children
        for(int i=0;i<n_sols;i++){
            assert(sols[i]->n_facs==current_size); // All solutions are expected to have the same size.
            for(int f=0;f<prob->n_facs;f++){
                futuresol *fsol = (futuresol *)(futuresols+fsol_size*n_futuresols);
                n_futuresols += futuresol_init_from(fsol,sols[i],f);
            }
        }
    }else{
        // Pick children at random
        shuffler *shuf = shuffler_init(prob->n_facs);
        for(int i=0;i<n_sols;i++){
            shuffler_reshuffle(shuf);
            int n_childs = 0;
            while(n_childs<branching){
                int f = (int) shuffler_next(shuf);
                futuresol *fsol = (futuresol *)(futuresols+fsol_size*n_futuresols);
                int new_found = futuresol_init_from(fsol,sols[i],f);
                n_childs     += new_found;
                n_futuresols += new_found;
            }
        }
        shuffler_free(shuf);
    }

    { // Sort the futuresols in order to detect the similar ones faster:
        qsort(futuresols,n_futuresols,fsol_size,futuresol_cmp);
        int n_futuresols2 = 0;
        if(n_futuresols>0){
            futuresol *last_fsol = (futuresol *)(futuresols);
            n_futuresols2 = 1;
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
                    if((run->filter>=BETTER_THAN_ALL_PARENTS) == is_better){
                        memmove(last_fsol,fsol,fsol_size);
                    }
                }
            }
        }
        // Update the futuresols, and realloc to reduce memory usage
        n_futuresols = n_futuresols2;
        futuresols = safe_realloc(futuresols,fsol_size*n_futuresols);
    }

    solution **out_sols = safe_malloc(sizeof(solution*)*n_futuresols);
    { // Create new solutions [in parallel]
        expand_thread_args *targs = safe_malloc(sizeof(expand_thread_args)*run->n_threads);
        for(int i=0;i<run->n_threads;i++){
            targs[i].thread_id = i;
            targs[i].run = run;
            targs[i].n_fsols = n_futuresols;
            targs[i].futuresols = futuresols;
            targs[i].fsol_size = fsol_size;
            targs[i].out_sols = out_sols;
        }
        // Generate threads in order to expand the solutions
        pthread_t *threads = safe_malloc(sizeof(pthread_t)*run->n_threads);
        for(int i=0;i<run->n_threads;i++){
            int rc = pthread_create(&threads[i],NULL,expand_thread_execution,&targs[i]);
            if(rc){
                fprintf(stderr,"Error %d on thread creation\n",rc);
                exit(1);
            }
        }
        // Join threads
        for(int i=0;i<run->n_threads;i++){
            pthread_join(threads[i],NULL);
        }
        //
        free(threads);
        free(targs);
    }

    // Set the terminal flag for the original solutions (revert with futuresols origins)
    for(int i=0;i<n_sols;i++) sols[i]->terminal = 1;

    { // Eliminate NULLed out solutions
        int n_sols = 0;
        for(int r=0;r<n_futuresols;r++){
            if(out_sols[r]!=NULL){
                futuresol *fsol = (futuresol *)(futuresols+fsol_size*r);
                fsol->origin->terminal = 0; // origin is not terminal because it had a child.
                out_sols[n_sols] = out_sols[r];
                n_sols += 1;
            }
        }
        out_sols = safe_realloc(out_sols,sizeof(solution*)*(n_sols));
        *out_n_sols = n_sols;
    }


    free(futuresols);
    return out_sols;
}
