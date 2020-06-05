#include "construction.h"

typedef struct {
    // Best solutions found so far, to be retrieved.
    solution **final; // NOTE: This array is expected to have a size of 2*prob->target_sols
    int n_final;
    // Saved solution for PR.
    solution **prpool;
    int n_prpool;
} solmemory;

void solmemory_merge_with_final(rundata *run, solmemory *solmem, solution **sols, int n_sols, int restart){
    // Save current restart best solution found
    for(int i=0;i<n_sols;i++){
        if(sols[i]->value > run->restart_values[restart]){
            run->restart_values[restart] = sols[i]->value;
        }
    }

    // Pick the best prob->target_sols candidates
    reduction_bests(run,sols,&n_sols,run->target_sols);
    // Merge candidates with the final solutions
    for(int i=0;i<n_sols;i++){
        solmem->final[solmem->n_final] = sols[i];
        solmem->n_final += 1;
    }
    // Remove repeated solutions
    solutions_delete_repeated(solmem->final,&solmem->n_final);
    // Pick the best prob->target_sols for final solutions
    reduction_bests(run,solmem->final,&solmem->n_final,run->target_sols);
    // Update the lower bound
    if(solmem->final[0]->value > run->lower_bound) run->lower_bound = solmem->final[0]->value;


    // Solution cands is not longer needed.
    free(sols);
}

// final is expected to have a size >= 2*prob->target_sols
// cands is liberated.
void update_final_solutions(rundata *run, solmemory *solmem,
        solution **cands, int n_cands, int csize, int restart){
    // The candidates should have the same number of solutions that csize
    assert(n_cands==0 || cands[0]->n_facs==csize);
    // Perform local search on the candidate solutions
    if(run->local_search){
        if(run->local_search_only_terminal && csize<run->prob->size_restriction_maximum){
            // Run local search just on the terminal solutions
            solution **terminals = safe_malloc(sizeof(solution *)*n_cands);
            int n_terminal = 0;
            for(int i=0;i<n_cands;i++){
                if(cands[i]->terminal){
                    terminals[n_terminal] = cands[i];
                    n_terminal += 1;
                }
            }
            if(n_terminal>0){
                if(run->verbose) printf("Performing LS on \033[34;1m%d\033[0m terminal solutions.\n",n_terminal);
                solutions_hill_climbing(run,terminals,n_terminal);
            }
            free(terminals);
        }else{
            // Run local search on all the solutions
            if(n_cands>0){
                if(run->verbose) printf("Performing LS on \033[34;1m%d\033[0m solutions.\n",n_cands);
                solutions_hill_climbing(run,cands,n_cands);
            }
        }
        // Delete repeated solutions after local search
        int n_cands0 = n_cands;
        solutions_delete_repeated(cands,&n_cands);
        if(n_cands0>n_cands){
            if(run->verbose) printf(
                "Reduced \033[34;1m%d\033[0m solutions to \033[34;1m%d\033[0m local optima.\n",
                n_cands0,n_cands);
        }
    }

    // Copy terminal solutions on the PR pool
    if(run->path_relinking){
        // Count number of terminal solutions
        int n_terminal = 0;
        for(int i=0;i<n_cands;i++){
            solution *sol = cands[i];
            if(sol->terminal || csize==run->prob->size_restriction_maximum) n_terminal += 1;
        }
        // Add terminal solutions to the PR pool
        solmem->prpool = safe_realloc(solmem->prpool,sizeof(solution *)*(solmem->n_prpool+n_terminal));
        for(int i=0;i<n_cands;i++){
            solution *sol = cands[i];
            if(sol->terminal || csize==run->prob->size_restriction_maximum){
                solmem->prpool[solmem->n_prpool] = solution_copy(run->prob,sol);
                solmem->n_prpool += 1;
            }
        }
    }

    // Save number of local optima found
    if(restart==0 && run->local_search){
        run->firstr_per_size_n_local_optima[csize] = n_cands;
    }

    // == Update array of best solutions
    solmemory_merge_with_final(run,solmem,cands,n_cands,restart);
}

solution **new_find_best_solutions(rundata *run, redstrategy *rstrats, int n_rstrats,
        int *out_n_sols){

    const problem *prob = run->prob;

    // Set random seed
    srand(run->random_seed);

    // The final solutions:
    solmemory solmem;
    solmem.n_final   = 0;
    solmem.final     = safe_malloc(sizeof(solution *)*run->target_sols*2);

    for(int r=0;r<run->n_restarts;r++){

        // Initialize PR pool
        if(run->path_relinking){
            solmem.n_prpool  = 0;
            solmem.prpool    = safe_malloc(sizeof(solution *)*1);
        }

        // Measure restart time
        clock_t restart_start = clock();

        int first_restart = r==0;

        if(run->verbose) printf("\n== RESTART %d/%d ==\n",r+1,run->n_restarts);

        // The previous generation
        int prev_n_sols = 1;
        solution **prev_sols = safe_malloc(sizeof(solution *)*prev_n_sols);
        prev_sols[0] = solution_empty(prob);

        int csize = 0; // Last base computed
        while(prev_n_sols>0){

            if(run->verbose) printf("\nBase has \033[31;1m%d\033[0m solutions of size \033[32;1m%d\033[0m.\n",prev_n_sols,csize);
            // Apply branch and bound
            if(run->branch_and_bound){
                int n_sols0 = prev_n_sols;
                branch_and_bound(run,prev_sols,&prev_n_sols);
                if(prev_n_sols < n_sols0){
                    if(run->verbose) printf("Pruned \033[31;1m%d\033[0m -> \033[31;1m%d\033[0m solutions, by B&B.\n",n_sols0,prev_n_sols);
                }
            }

            // Save number of solutions after expansion
            if(first_restart && run->local_search){
                run->firstr_per_size_n_sols[csize] = prev_n_sols;
            }

            // Apply the reduction strategies
            for(int i=0;i<n_rstrats;i++){
                reduce_by_redstrategy(run,rstrats[i],prev_sols,&prev_n_sols);
            }


            // Save number of solutions after reduction
            if(first_restart && run->local_search){
                run->firstr_per_size_n_sols_after_red[csize] = prev_n_sols;
            }

            int next_n_sols = 0;
            solution **next_sols = NULL;

            // Expand solutions from the previous generation
            if(csize<prob->n_facs && prev_n_sols>0){
                if(prob->size_restriction_maximum==-1 || csize<prob->size_restriction_maximum){
                    if(run->verbose) printf("Expanding \033[31;1m%d\033[0m solutions.\n",prev_n_sols);
                    next_n_sols = prev_n_sols;
                    int pool_size = n_rstrats>0? rstrats[n_rstrats-1].n_target : prob->n_facs;
                    next_sols = new_expand_solutions(run,prev_sols,prev_n_sols,&next_n_sols,pool_size);
                }
            }


            // Put the prev generation after LS in the final solutions
            if(prob->size_restriction_minimum==-1 || csize>=prob->size_restriction_minimum){
                update_final_solutions(run,&solmem,prev_sols,prev_n_sols,csize,r);
            }else{
                // Free solution memory
                for(int i=0;i<prev_n_sols;i++) solution_free(prev_sols[i]);
                free(prev_sols);
            }

            // Now the current gen is the previous one
            prev_n_sols = next_n_sols;
            prev_sols = next_sols;

            // Increase csize
            csize += 1;
            // Update number of iterations
            if(first_restart) run->firstr_n_iterations = csize;
        }
        run->total_n_iterations += csize;

        // Turn the PR pool into
        if(run->path_relinking){

            // Perform path relinking on the terminal solutions
            if(run->verbose) printf("\nPerforming Path Relinking on \033[34;1m%d\033[0m terminal solutions.\n",solmem.n_prpool);
            solutions_path_relinking(run,&solmem.prpool,&solmem.n_prpool);

            if(run->verbose) printf("PR resulted in \033[34;1m%d\033[0m different solutions.\n",solmem.n_prpool);

            // Perform local search on the resulting solutions
            if(run->local_search){
                if(run->verbose) printf("Performing LS on \033[34;1m%d\033[0m resulting solutions.\n",solmem.n_prpool);
                solutions_hill_climbing(run,solmem.prpool,solmem.n_prpool);
                // Delete repeated solutions after local search
                int n_prpool0 = solmem.n_prpool;
                solutions_delete_repeated(solmem.prpool,&solmem.n_prpool);
                if(n_prpool0>solmem.n_prpool){
                    if(run->verbose) printf(
                        "Reduced \033[34;1m%d\033[0m temrinal solutions to \033[34;1m%d\033[0m local optima.\n",
                        n_prpool0,solmem.n_prpool);
                }
            }


            // == Update array of best solutions
            solmemory_merge_with_final(run,&solmem,solmem.prpool,solmem.n_prpool,r);

        }

        clock_t restart_end = clock();
        run->restart_times[r] = (double)(restart_end - restart_start) / (double)CLOCKS_PER_SEC;

        if(run->verbose) printf("\n");

    }


    // Retrieve the final solutions:
    *out_n_sols = solmem.n_final;
    return solmem.final;
}
