#include "construction.h"

typedef struct {
    // Best solutions found so far, to be retrieved.
    solution **final; // NOTE: This array is expected to have a size of 2*prob->target_sols
    int n_final;
    // Solutions that have been selected in the current iteration
    solution **selectpool;
    int n_selectpool;
} solmemory;

// Apply local search on selected solutions and then delete repeated
void solmem_local_search_and_delete_repeated(rundata *run, solmemory *solmem){
    if(run->verbose) printf("Performing LS on \033[34;1m%d\033[0m selected solutions.\n",solmem->n_selectpool);
    solutions_hill_climbing(run,solmem->selectpool,solmem->n_selectpool);
    // Delete repeated solutions after local search
    int n_prpool0 = solmem->n_selectpool;
    solutions_sort_and_delete_repeated(solmem->selectpool,&solmem->n_selectpool);
    if(n_prpool0>solmem->n_selectpool){
        if(run->verbose) printf(
            "Reduced \033[34;1m%d\033[0m selected solutions to \033[34;1m%d\033[0m local optima.\n",
            n_prpool0,solmem->n_selectpool);
    }
}

void update_selected_solutions(rundata *run, solmemory *solmem,
        solution **cands, int n_cands, int csize, int restart){
    // The candidates should have the same number of solutions, csize.
    assert(n_cands==0 || cands[0]->n_facs==csize);
    // Check if the solution meets the conditions to be selected
    int min_size_ok = run->prob->size_restriction_minimum==-1 || csize>=run->prob->size_restriction_minimum;
    int max_size_ok = run->prob->size_restriction_maximum==-1 || csize<=run->prob->size_restriction_maximum;
    if(min_size_ok && max_size_ok){
        // If solution is selected
        int *sel = safe_malloc(sizeof(int)*n_cands);
        // Count number of new selected solutions
        int n_sel = 0;
        for(int i=0;i<n_cands;i++){
            solution *sol = cands[i];
            sel[i] = !run->select_only_terminal || sol->terminal || csize==run->prob->size_restriction_maximum;
            n_sel += sel[i];
        }
        // Add terminal solutions to the selected pool
        solmem->selectpool = safe_realloc(solmem->selectpool,sizeof(solution *)*(solmem->n_selectpool+n_sel));
        for(int i=0;i<n_cands;i++){
            solution *sol = cands[i];
            if(sel[i]){
                solmem->selectpool[solmem->n_selectpool] = sol;
                solmem->n_selectpool += 1;
            }else{
                solution_free(sol);
            }
        }
        free(sel);
    }else{
        // Free all solutions
        for(int i=0;i<n_cands;i++) solution_free(cands[i]);
    }
    // Free candidate list
    free(cands);
}

// Save selected solutions (copying) in the final array if they are better than the current ones on it.
// Retrieves value of the best solution
double solmemory_register_in_final(rundata *run, solmemory *solmem, int restart, int only_best){
    // Find new best solution value
    solution *new_best_solution = NULL;
    double new_best_sol_value   = -INFINITY;
    for(int i=0;i<solmem->n_selectpool;i++){
        if(solmem->selectpool[i]->value > new_best_sol_value){
            new_best_sol_value = solmem->selectpool[i]->value;
            new_best_solution  = solmem->selectpool[i];
        }
    }

    if(new_best_solution!=NULL){
        // Save value of current restart best solution found
        if(new_best_solution->value > run->run_inf->restart_values[restart]){
            run->run_inf->restart_values[restart] = new_best_solution->value;
        }

        if(only_best){
            // Add best solution to the final solutions
            solmem->final[solmem->n_final] = solution_copy(run->prob,new_best_solution);
            solmem->n_final += 1;
        }else{
            // Copy all selected solutions
            int n_sel_cpy = solmem->n_selectpool;
            solution **sel_cpy = safe_malloc(sizeof(solution *)*solmem->n_selectpool);
            for(int i=0;i<n_sel_cpy;i++){
                sel_cpy[i] = solution_copy(run->prob,solmem->selectpool[i]);
            }

            // Pick the best prob->target_sols candidates
            reduction_bests(run,sel_cpy,&n_sel_cpy,run->target_sols);

            // Merge candidates with the final solutions
            for(int i=0;i<n_sel_cpy;i++){
                solmem->final[solmem->n_final] = sel_cpy[i];
                solmem->n_final += 1;
            }

            // Copy of selecteds is not longer needed.
            free(sel_cpy);
        }

        // Remove repeated solutions
        solutions_sort_and_delete_repeated(solmem->final,&solmem->n_final);
        // Pick the best prob->target_sols for final solutions
        reduction_bests(run,solmem->final,&solmem->n_final,run->target_sols);
        // Update the lower bound
        if(solmem->final[0]->value > run->bnb_lower_bound) run->bnb_lower_bound = solmem->final[0]->value;
    }

    return new_best_sol_value;
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

        // Initialize terminal pool
        solmem.n_selectpool  = 0;
        solmem.selectpool    = safe_malloc(sizeof(solution *)*1);

        // Measure restart time
        clock_t restart_start = clock();

        int first_restart = r==0;

        if(run->verbose) printf("\n== RESTART %d/%d ==\n",r+1,run->n_restarts);

        // The previous generation
        int prev_n_sols = 1;
        solution **prev_sols = safe_malloc(sizeof(solution *)*prev_n_sols);
        prev_sols[0] = solution_empty(prob);

        int csize = 0; // Current solution size, last base computed
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
                run->run_inf->firstr_per_size_n_sols[csize] = prev_n_sols;
            }

            // Apply the reduction strategies
            for(int i=0;i<n_rstrats;i++){
                if(!rstrats[i].for_selected_sols){ // Only apply reductions not intended for selected solutions
                    reduce_by_redstrategy(run,rstrats[i],prev_sols,&prev_n_sols);
                }
            }

            // Save number of solutions after reduction
            if(first_restart && run->local_search){
                run->run_inf->firstr_per_size_n_sols_after_red[csize] = prev_n_sols;
            }

            // Expand solutions from the previous generation to create the next one
            int next_n_sols = 0;
            solution **next_sols = NULL;

            if(csize<prob->n_facs && prev_n_sols>0){
                if(prob->size_restriction_maximum==-1 || csize<prob->size_restriction_maximum){
                    if(run->verbose) printf("Expanding \033[31;1m%d\033[0m solutions.\n",prev_n_sols);
                    next_n_sols = prev_n_sols;
                    // Get the final pool size expected at the end of each iteration
                    int pool_size = prev_n_sols*prob->n_facs;
                    for(int s=0;s<n_rstrats;s++){
                        if(!rstrats[s].for_selected_sols) pool_size = rstrats[s].n_target;
                    }
                    // Expand solutions to get the next generation
                    next_sols = new_expand_solutions(run,prev_sols,prev_n_sols,&next_n_sols,pool_size);
                }
            }

            // Add the prev generation solutions in the selected solutions and free their memory
            update_selected_solutions(run,&solmem,prev_sols,prev_n_sols,csize,r);

            // Now the current gen is the previous one
            prev_n_sols = next_n_sols;
            prev_sols = next_sols;

            // Increase csize
            csize += 1;
            // Update number of iterations
            if(first_restart) run->run_inf->firstr_n_iterations = csize;

        }
        run->run_inf->total_n_iterations += csize;

        if(run->verbose) printf("\nStarting final reduction with \033[34;1m%d\033[0m selected solutions:\n",solmem.n_selectpool);

        { /* Operate on selected solutions */

            // Apply local search to selected solutions if we are before final selections
            if(run->local_search!=NO_LOCAL_SEARCH && run->local_search_before_select){
                solmem_local_search_and_delete_repeated(run,&solmem);
            }

            // Save best solution, just in case it is deleted in the reductions
            solmemory_register_in_final(run,&solmem,r,1);

            // Apply reduction to selected solutions
            for(int i=0;i<n_rstrats;i++){
                if(rstrats[i].for_selected_sols){
                    reduce_by_redstrategy(run,rstrats[i],solmem.selectpool,&solmem.n_selectpool);
                }
            }

            // Apply local search to selected solutions if we are before final selections
            if(run->local_search!=NO_LOCAL_SEARCH && !run->local_search_before_select){
                solmem_local_search_and_delete_repeated(run,&solmem);
            }
        }

        // Apply post-optimization to the selected solutions
        if(run->path_relinking!=NO_PATH_RELINKING){

            if(run->verbose) printf("\nStarting Post-optimization with \033[34;1m%d\033[0m selected solutions:\n",solmem.n_selectpool);

            while(1){

                if(run->verbose) printf("\n");

                // Find current best solution value
                double prev_best_sol_value = solmemory_register_in_final(run,&solmem,r,1);

                // Apply the reduction strategies
                for(int i=0;i<n_rstrats;i++){
                    if(rstrats[i].for_selected_sols){ // Only apply reductions not intended for path relinking
                        reduce_by_redstrategy(run,rstrats[i],solmem.selectpool,&solmem.n_selectpool);
                    }
                }

                // Perform path relinking on the terminal solutions
                if(run->verbose) printf("Performing Path Relinking on \033[34;1m%d\033[0m solutions.\n",solmem.n_selectpool);
                solutions_path_relinking(run,&solmem.selectpool,&solmem.n_selectpool);

                if(run->verbose) printf("PR resulted in \033[34;1m%d\033[0m different solutions.\n",solmem.n_selectpool);


                // Perform local search on the resulting solutions
                if(run->local_search != NO_LOCAL_SEARCH){
                    if(run->verbose) printf("Performing LS on \033[34;1m%d\033[0m resulting solutions.\n",solmem.n_selectpool);
                    solutions_hill_climbing(run,solmem.selectpool,solmem.n_selectpool);

                    // Delete repeated solutions after local search
                    int n_prpool0 = solmem.n_selectpool;
                    solutions_sort_and_delete_repeated(solmem.selectpool,&solmem.n_selectpool);
                    if(n_prpool0>solmem.n_selectpool){
                        if(run->verbose) printf(
                            "Reduced \033[34;1m%d\033[0m terminal solutions to \033[34;1m%d\033[0m local optima.\n",
                            n_prpool0,solmem.n_selectpool);
                    }
                }

                // Find new best solution value
                double new_best_sol_value = solmemory_register_in_final(run,&solmem,r,1);

                // Check for terminating conditions and save best solution found on this iteration
                if(solmem.n_selectpool<=1) break;
                // if(solmem.n_selectpool<=n_prpool_before_pr-1) break; // NOTE: Can this help?
                if(run->path_relinking==PATH_RELINKING_1_ITER){
                    // Make PR only happen once with PATH_RELINKING_1_STEP
                    break;
                }else if(run->path_relinking==PATH_RELINKING_UNTIL_NO_BETTER){
                    // No better solution was found
                    if(new_best_sol_value<=prev_best_sol_value) break;
                }

            }
        }

        // Save best solutions in the final ones
        solmemory_register_in_final(run,&solmem,r,0);

        // Free selected pool
        for(int i=0;i<solmem.n_selectpool;i++){
            solution_free(solmem.selectpool[i]);
        }
        free(solmem.selectpool);

        clock_t restart_end = clock();
        run->run_inf->restart_times[r] = (double)(restart_end - restart_start) / (double)CLOCKS_PER_SEC;

        if(run->verbose) printf("\n");


    }
    // Retrieve the final solutions:
    *out_n_sols = solmem.n_final;

    #ifdef DEBUG
        // Check final solutions integrity
        for(int i=0;i<solmem.n_final;i++){
            solution *sol = solmem.final[i];
            assert(solution_check_integrity(prob,sol));
        }
    #endif

    return solmem.final;
}
