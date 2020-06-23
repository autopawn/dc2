#include "construction.h"

typedef struct {
    // Best solutions found so far, to be retrieved.
    solution **final; // NOTE: This array is expected to have a size of 2*prob->target_sols
    int n_final;
    // Saved solution for PR.
    solution **prpool;
    int n_prpool;
} solmemory;

// free solutions in the sols array, but some may be saved in the solmem final solutions
// if they are better than the current ones on it.
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
    solutions_sort_and_delete_repeated(solmem->final,&solmem->n_final);
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
        if(run->local_search_only_terminal && (csize<run->prob->size_restriction_maximum || run->prob->size_restriction_maximum==-1)){
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
        solutions_sort_and_delete_repeated(cands,&n_cands);
    }

    // Copy terminal solutions on the PR pool
    if(run->path_relinking!=NO_PATH_RELINKING){
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
        if(run->path_relinking!=NO_PATH_RELINKING){
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
                if(!rstrats[i].for_path_relinking){ // Only apply reductions not intended for path relinking
                    reduce_by_redstrategy(run,rstrats[i],prev_sols,&prev_n_sols);
                }
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
                    // Get the final pool size expected at the end of each iteration
                    int pool_size = prob->n_facs;
                    for(int s=0;s<n_rstrats;s++){
                        if(!rstrats[s].for_path_relinking) pool_size = rstrats[s].n_target;
                    }
                    // Expand solutions to get the next generation
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

        // Apply post-optimization to the PR pool
        if(run->path_relinking!=NO_PATH_RELINKING){

            if(run->verbose) printf("\nPerforming Post-optimization on \033[34;1m%d\033[0m terminal solutions:\n",solmem.n_prpool);

            while(1){

                if(run->verbose) printf("\n");

                // Find current best solution value
                double prev_best_sol_value = -INFINITY;
                for(int i=0;i<solmem.n_prpool;i++){
                    if(solmem.prpool[i]->value > prev_best_sol_value) prev_best_sol_value = solmem.prpool[i]->value;
                }

                // Apply the reduction strategies
                for(int i=0;i<n_rstrats;i++){
                    if(rstrats[i].for_path_relinking){ // Only apply reductions not intended for path relinking
                        reduce_by_redstrategy(run,rstrats[i],solmem.prpool,&solmem.n_prpool);
                    }
                }

                // Perform path relinking on the terminal solutions
                int n_prpool_before_pr = solmem.n_prpool;
                if(run->verbose) printf("Performing Path Relinking on \033[34;1m%d\033[0m solutions.\n",solmem.n_prpool);
                solutions_path_relinking(run,&solmem.prpool,&solmem.n_prpool);

                if(run->verbose) printf("PR resulted in \033[34;1m%d\033[0m different solutions.\n",solmem.n_prpool);

                // Perform local search on the resulting solutions
                if(run->local_search){
                    if(run->verbose) printf("Performing LS on \033[34;1m%d\033[0m resulting solutions.\n",solmem.n_prpool);
                    solutions_hill_climbing(run,solmem.prpool,solmem.n_prpool);
                    // Delete repeated solutions after local search
                    int n_prpool0 = solmem.n_prpool;
                    solutions_sort_and_delete_repeated(solmem.prpool,&solmem.n_prpool);
                    if(n_prpool0>solmem.n_prpool){
                        if(run->verbose) printf(
                            "Reduced \033[34;1m%d\033[0m terminal solutions to \033[34;1m%d\033[0m local optima.\n",
                            n_prpool0,solmem.n_prpool);
                    }
                }

                // // Delete a fraction of the worst solutions
                // float remov_factor = 0.33333;
                // int n_prpool_before = solmem.n_prpool;
                // reduction_remove_worst(solmem.prpool,&solmem.n_prpool,remov_factor);
                // if(run->verbose) printf("Deleted the worst %d%% of the solutions: \033[34;1m%d\033[0m -> \033[34;1m%d\033[0m.\n",
                //     (int)(100*remov_factor),n_prpool_before,solmem.n_prpool);

                // Find new best solution value
                solution *new_best_solution = NULL;
                double new_best_sol_value   = -INFINITY;
                for(int i=0;i<solmem.n_prpool;i++){
                    if(solmem.prpool[i]->value > new_best_sol_value){
                        new_best_sol_value = solmem.prpool[i]->value;
                        new_best_solution  = solmem.prpool[i];
                    }
                }
                // Add best solution to final solutions, just in case the reduction process deletes it
                if(new_best_solution){
                    solution **singleton_best = safe_malloc(sizeof(solution *)*1);
                    singleton_best[0] = solution_copy(run->prob,new_best_solution);
                    solmemory_merge_with_final(run,&solmem,singleton_best,1,r);
                }
                printf("current: %f\n",new_best_sol_value);

                // Check for terminating conditions and save best solution found on this iteration
                if(solmem.n_prpool<=1) break;
                if(solmem.n_prpool<=n_prpool_before_pr-1) break;
                if(run->path_relinking==PATH_RELINKING_1_STEP){
                    // Make PR only happen once with PATH_RELINKING_1_STEP
                    break;
                }else if(run->path_relinking==PATH_RELINKING_UNTIL_NO_BETTER){
                    // No better solution was found
                    if(new_best_sol_value<=prev_best_sol_value) break;
                }

            }

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
