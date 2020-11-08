#include "localsearch.h"

void update_phi1_and_phi2(const problem *prob, const solution *sol, int f_ins, int f_rem,
        int *phi1, int *phi2, int *affected_mask){
    for(int i=0;i<prob->n_clis;i++){
        #ifdef DEBUG
            int before_phi1 = phi1[i];
            int before_phi2 = phi2[i];
        #endif
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
        #ifdef DEBUG
            assert(phi1[i]!=phi2[i] || phi1[i]==-1);
            assert(problem_assig_value(prob,phi1[i],i) >= problem_assig_value(prob,phi2[i],i));
            assert(affected_mask==NULL || affected_mask[i] || (before_phi1==phi1[i] && before_phi2==phi2[i]));
        #endif
    }
}

void solutions_sort_and_delete_repeated(solution **sols, int *n_sols){
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
            // Delete sols[i] as it is the same as sols[i-1]
            solution_free(sols[i]);
        }
    }
    *n_sols = n_final;
}


// ============================================================================
// ======== LOCAL SEARCH
// ============================================================================

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
                args->n_moves += solution_resendewerneck_hill_climbing(args->run,&args->sols[r],NULL,mat);
            }
            fastmat_free(mat);
        }
    }else{
        for(int r=args->thread_id;r<args->n_sols;r+=args->run->n_threads){
            // Perform local search on the given solution
            args->n_moves += solution_whitaker_hill_climbing(args->run,&args->sols[r],NULL,args->shuff);
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
    run->run_inf->n_local_searches += n_sols;
    run->run_inf->n_local_search_movements += n_moves;
    run->run_inf->local_search_seconds += seconds;
}


// ============================================================================
// ======== PATH RELINKING
// ============================================================================

typedef struct {
    int thread_id;
    const rundata *run;
    solution **pool;
    int n_pool;
    solution **result;
    shuffler *shuff;
} path_relinking_thread_args;

void *path_relinking_thread_execution(void *arg){
    path_relinking_thread_args *args = (path_relinking_thread_args *) arg;

    // Allocate a unique reusable fastmat if SWAP_RESENDE_WERNECK
    fastmat *mat = NULL;
    if(args->run->local_search_pr==SWAP_RESENDE_WERNECK){
        mat = fastmat_init(args->run->prob->n_facs,args->run->prob->n_facs);
    }

    int c_pair = 0;
    for(int i=0;i<args->n_pool;i++){
        for(int j=i+1;j<args->n_pool;j++){
            // Check if this pair should be linked by the current thread
            if(c_pair%args->run->n_threads==args->thread_id){

                // Pick initial and ending solution from the pair according to solution value
                const solution *sol_ini, *sol_end;
                if(args->pool[i]->value >= args->pool[j]->value){
                    sol_ini = args->pool[i];
                    sol_end = args->pool[j];
                }else{
                    sol_ini = args->pool[j];
                    sol_end = args->pool[i];
                }

                // Perform path relinking
                solution *sol = solution_copy(args->run->prob,sol_ini);
                if(args->run->local_search_pr==SWAP_RESENDE_WERNECK){
                    solution_resendewerneck_hill_climbing(args->run,&sol,sol_end,mat);
                }else{
                    solution_whitaker_hill_climbing(args->run,&sol,sol_end,args->shuff);
                }

                args->result[c_pair] = sol;

                // Chack that path relinking was performed correctly
                #ifdef DEBUG
                    // Check that all facilitites in sol came from one of the solutions
                    for(int k=0; k<sol->n_facs; k++){
                        int in_ini = elem_in_sorted(sol_ini->facs,sol_ini->n_facs,sol->facs[k]);
                        int in_end = elem_in_sorted(sol_end->facs,sol_end->n_facs,sol->facs[k]);
                        assert(in_ini || in_end);
                    }
                    // Check that all facilities in both solutions remain in sol
                    for(int k=0; k<sol_ini->n_facs; k++){
                        int f = sol_ini->facs[k];
                        if(elem_in_sorted(sol_end->facs, sol_end->n_facs, f)){
                            assert(elem_in_sorted(sol->facs,sol->n_facs,f));
                        }
                    }
                #endif
            }
            c_pair += 1;
        }
    }

    assert(c_pair== args->n_pool*(args->n_pool-1)/2);

    if(mat) fastmat_free(mat);


    return NULL;
}

// Perform path relinking searches from the current solutions that won't be modified (in parallel).
void solutions_path_relinking(rundata *run, solution ***sols, int *n_sols){
    // Start measuring time
    clock_t start = clock();
    // Allocate memory for resulting set of solutions
    int n_resulting = (*n_sols)*((*n_sols)-1)/2;
    solution **resulting = safe_malloc(sizeof(solution*)*n_resulting);

    // Allocate memory for threads and arguments
    pthread_t *threads = safe_malloc(sizeof(pthread_t)*run->n_threads);
    path_relinking_thread_args *targs = safe_malloc(sizeof(path_relinking_thread_args)*run->n_threads);

    // Call all threads to perform local search
    for(int i=0;i<run->n_threads;i++){
        // Set arguments for the thread
        targs[i].thread_id = i;
        targs[i].run = run;
        targs[i].pool = (*sols);
        targs[i].n_pool = (*n_sols);
        targs[i].result = resulting;
        if(run->local_search_pr==SWAP_FIRST_IMPROVEMENT){
            targs[i].shuff = shuffler_init(run->prob->n_facs);
        }else{
            targs[i].shuff = NULL;
        }

        // Create thread
        int rc = pthread_create(&threads[i],NULL,path_relinking_thread_execution,&targs[i]);
        if(rc){
            fprintf(stderr,"ERROR: Error %d on pthread_create\n",rc);
            exit(1);
        }
    }

    // Join threads
    for(int i=0;i<run->n_threads;i++){
        pthread_join(threads[i],NULL);
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
    run->run_inf->path_relinking_seconds += seconds;

    // Delete similar solutions
    solutions_sort_and_delete_repeated(resulting,&n_resulting);

    // Delete original solutions
    for(int i=0;i<(*n_sols);i++){
        solution_free((*sols)[i]);
    }
    free(*sols);

    // Replace the original solution array for the array of resulting solutions
    *sols = resulting;
    *n_sols = n_resulting;
}


// ============================================================================
// ======== AVAIL MOVES
// ============================================================================

availmoves *availmoves_init(const problem *prob, const solution *sol, const solution *tgt){
    // Allowed insertions
    int *inss = safe_malloc(sizeof(int)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++) inss[i] = 1;
    for(int k=0;k<sol->n_facs;k++)  inss[sol->facs[k]] = 0; // can't insert if already in solution

    // Allowed removals
    int *rems = safe_malloc(sizeof(int)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++) rems[i] = 0;
    for(int k=0;k<sol->n_facs;k++)  rems[sol->facs[k]] = 1; // can only remove if already present

    int *used = safe_malloc(sizeof(int)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++) used[i] = 0;
    for(int k=0;k<sol->n_facs;k++)  used[sol->facs[k]] = 1;

    // Restrict movements more if tgt
    if(tgt){
        // Create array to directly know if tgt has a facility
        int *tgt_used = safe_malloc(sizeof(int)*prob->n_facs);
        for(int i=0;i<prob->n_facs;i++) tgt_used[i] = 0;
        for(int k=0;k<tgt->n_facs;k++)  tgt_used[tgt->facs[k]] = 1;
        // Update inss and rems
        for(int i=0;i<prob->n_facs;i++){
            inss[i] = inss[i] && tgt_used[i]; // can only insert if tgt has it
            rems[i] = rems[i] && !tgt_used[i]; // can only remove if tgt doens't have it
        }
        //
        free(tgt_used);
    }

    // == Initialize availmoves
    availmoves *av = safe_malloc(sizeof(availmoves));

    av->avail_inss = inss;
    av->avail_rems = rems;

    av->n_insertions = 0;
    av->insertions = safe_malloc(sizeof(int)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++){
        if(inss[i]){
            av->insertions[av->n_insertions] = i;
            av->n_insertions++;
        }
    }

    av->n_removals = 0;
    av->removals = safe_malloc(sizeof(int)*prob->n_facs);
    for(int i=0;i<prob->n_facs;i++){
        if(rems[i]){
            av->removals[av->n_removals] = i;
            av->n_removals++;
        }
    }

    av->used = used;

    // Won't recover indexes unless tgt is NULL
    av->path_relinking = (tgt!=NULL);

    return av;
}

void availmoves_register_move(availmoves *av, int f_ins, int f_rem){
    // Update used array
    if(f_ins>=0){
        assert(!av->used[f_ins]);
        av->used[f_ins] = 1;
    }
    if(f_rem>=0){
        assert(av->used[f_rem]);
        av->used[f_rem] = 0;
    }

    // Remove insertion option
    if(f_ins>=0){

        assert(av->avail_inss[f_ins]);
        av->avail_inss[f_ins] = 0;

        int ins_indx = -1;
        for(int i=0;i<av->n_insertions;i++){
            if(av->insertions[i]==f_ins){
                ins_indx = i;
                break;
            }
        }
        assert(ins_indx>=0);
        av->insertions[ins_indx] = av->insertions[av->n_insertions-1];
        av->n_insertions -= 1;
    }

    // Remove removal option
    if(f_rem>=0){
        assert(av->avail_rems[f_rem]);
        av->avail_rems[f_rem] = 0;

        int rem_indx = -1;
        for(int i=0;i<av->n_removals;i++){
            if(av->removals[i]==f_rem){
                rem_indx = i;
                break;
            }
        }
        assert(rem_indx>=0);
        av->removals[rem_indx] = av->removals[av->n_removals-1];
        av->n_removals -= 1;

    }

    // Add new options if recovering indexes is allowed
    if(!av->path_relinking){
        // Add removal option for insertion
        if(f_ins>=0){
            assert(!av->avail_rems[f_ins]);
            av->avail_rems[f_ins] = 1;

            av->removals[av->n_removals] = f_ins;
            av->n_removals += 1;
        }
        // Add insertion option for removal
        if(f_rem>=0){
            assert(!av->avail_inss[f_rem]);
            av->avail_inss[f_rem] = 1;

            av->insertions[av->n_insertions] = f_rem;
            av->n_insertions += 1;
        }
    }
}

void availmoves_free(availmoves *av){
    free(av->removals);
    free(av->insertions);
    free(av->used);
    free(av->avail_rems);
    free(av->avail_inss);
    free(av);
}
