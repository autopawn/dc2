#include "construction.h"

// final is expected to have a size >= 2*prob->target_sols 
// cands is liberated.
void update_final_solutions(problem *prob, solution **final, int *n_final,
        solution **cands, int n_cands){
    // Perform local search on the candidate solutions
    if(n_cands>0){
        printf("Performing LS on \033[31;1m%d\033[0m solutions of size \033[32;1m%d\033[0m, saving the bests for final.\n",
            n_cands,cands[0]->n_facs);
    }
    for(int i=0;i<n_cands;i++){
        solution_whitaker_hill_climbing(prob,cands[i]);
    }
    // Delete repeated solutions after local search
    int n_cands0 = n_cands;
    solutions_delete_repeated(prob,cands,&n_cands);
    if(n_cands0>n_cands){
        printf("Reduced \033[31;1m%d\033[0m solutions to \033[32;1m%d\033[0m local optima.\n",
            n_cands0,n_cands);
    }
    // Pick the best prob->target_sols candidates
    reduction_bests(prob,cands,&n_cands,prob->target_sols);
    // Merge candidates with the final solutions
    for(int i=0;i<n_cands;i++){
        final[*n_final] = cands[i];
        *n_final += 1;
    }
    // Pick the best prob->target_sols for final solutions
    reduction_bests(prob,final,n_final,prob->target_sols);
    // Update the lower bound
    if(final[0]->value > prob->lower_bound) prob->lower_bound = final[0]->value;
    // Solution cands is not longer needed.
    free(cands);
}

solution **new_find_best_solutions(problem *prob, redstrategy *rstrats, int n_rstrats,
        int *out_n_sols, int *out_n_iterations){

    // Perform problem precomputations
    printf("\n");
    printf("Performing precomputations.\n");
    problem_precompute(prob,rstrats,n_rstrats);
    
    // The final solutions:
    int final_n_sols = 0;
    solution **final_sols = safe_malloc(sizeof(solution *)*prob->target_sols*2);
    
    // The previous generation
    int prev_n_sols = 1;
    solution **prev_sols = safe_malloc(sizeof(solution *)*prev_n_sols);
    prev_sols[0] = solution_empty(prob);

    int csize = 0; // Last size computed
    while(csize<prob->n_facs && (prob->size_restriction==-1 || csize<prob->size_restriction)){
        printf("\n");
        printf("Expanding \033[31;1m%d\033[0m solutions of size \033[32;1m%d\033[0m\n",prev_n_sols,csize);
        
        // Expand solutions from the previous generation
        int next_n_sols = prev_n_sols;
        solution **next_sols = new_expand_solutions(prob,prev_sols,prev_n_sols,&next_n_sols);

        // Put the prev generation after LS in the final solutions
        update_final_solutions(prob,final_sols,&final_n_sols,prev_sols,prev_n_sols);

        // Increase csize
        csize += 1;

        // Apply branch and bound
        if(prob->branch_and_bound){
            int n_sols0 = next_n_sols;
            branch_and_bound(prob,next_sols,&next_n_sols);
            if(next_n_sols < n_sols0){
                printf("Pruned \033[31;1m%d\033[0m -> \033[31;1m%d\033[0m solutions, by B&B.\n",n_sols0,next_n_sols);
            }
        }

        // Apply the reduction strategies
        for(int i=0;i<n_rstrats;i++){
            reduce_by_redstrategy(prob,rstrats[i],next_sols,&next_n_sols);
        }

        // Now the current gen is the previous one
        prev_n_sols = next_n_sols;
        prev_sols = next_sols;
        
        if(prev_n_sols==0) break;
    }

    printf("\n");
    // Put the last generation after LS in the final solutions
    update_final_solutions(prob,final_sols,&final_n_sols,prev_sols,prev_n_sols);
    printf("\n");

    // Retrieve the final solutions:
    *out_n_sols = final_n_sols;
    *out_n_iterations = csize;
    return final_sols;
}