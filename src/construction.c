#include "construction.h"

solution **new_find_best_solutions(problem *prob, redstrategy *rstrats, int n_rstrats,
        int *out_n_sols, int *out_n_iterations){
    solution ***sols_by_size = safe_malloc(sizeof(solution **)*(prob->n_facs+2));
    int *n_sols = safe_malloc(sizeof(int)*(prob->n_facs+2));

    // Initialize first array with just the empty solution
    n_sols[0] = 1;
    sols_by_size[0] = safe_malloc(sizeof(solution *)*1);
    sols_by_size[0][0] = solution_empty(prob);


    int csize = 0; // Current size
    while((prob->size_restriction==-1 || csize<prob->size_restriction) && csize<prob->n_facs){
        printf("\n");
        // Expand solutions from the previous generation
        printf("Expanding %5d solutions of size %d\n",n_sols[csize],csize);
        sols_by_size[csize+1] = new_expand_solutions(
            prob,sols_by_size[csize],n_sols[csize],&n_sols[csize+1]);
        // Once the expansion is done, delete some solutions of the previous gen to save memory
        reduction_bests(prob,sols_by_size[csize],&n_sols[csize],prob->target_sols);

        // Sort new solutions by decreasing value
        csize += 1;
        qsort(sols_by_size[csize],n_sols[csize],sizeof(solution *),solutionp_value_cmp_inv);

        // Apply the reduction strategies
        for(int i=0;i<n_rstrats;i++){
            redstrategy_reduce(prob,rstrats[i],sols_by_size[csize],&n_sols[csize]);
        }
        if(n_sols[csize]==0) break;
    }
    // Delete some solutions of the previous gen for consistency
    reduction_bests(prob,sols_by_size[csize],&n_sols[csize],prob->target_sols);

    printf("\n");
    // Merge pools
    printf("Picking the best %5d solutions.\n",prob->target_sols);
    solution **merged_sols = safe_malloc(sizeof(solution *)*prob->target_sols*(csize+1));
    int final_n_sols = 0;
    for(int i=0;i<=csize;i++){
        for(int j=0;j<n_sols[i];j++){
            merged_sols[final_n_sols] = sols_by_size[i][j];
            final_n_sols += 1;
        }
        // Release array of references.
        free(sols_by_size[i]);
    }
    free(n_sols);
    free(sols_by_size);

    // Sort merged pool, pick best solutions.
    qsort(merged_sols,final_n_sols,sizeof(solution *),solutionp_value_cmp_inv);
    reduction_bests(prob,merged_sols,&final_n_sols,prob->target_sols);

    *out_n_sols = final_n_sols;
    *out_n_iterations = csize;
    return merged_sols;
}