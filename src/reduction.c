#include "reduction.h"

void reduce_by_redstrategy(problem *prob, const redstrategy rstrat, 
        solution **sols, int *n_sols){
    if(*n_sols<=rstrat.n_target) return;

    printf("Reducing \033[31;1m%d\033[0m -> \033[31;1m%d\033[0m solutions, ",
        *n_sols,rstrat.n_target);
    if(rstrat.method==REDUCTION_BESTS){
        printf("selecting bests.\n");
        reduction_bests(prob,sols,n_sols,rstrat.n_target);
    }
    else if(rstrat.method==REDUCTION_RANDOM_UNIFORM){
        printf("randomly (uniform).\n");
        reduction_random_uniform(prob,sols,n_sols,rstrat.n_target);
    }
    else if(rstrat.method==REDUCTION_RANDOM_RANK){
        printf("randomly (by rank).\n");
        reduction_random_rank(prob,sols,n_sols,rstrat.n_target);
    }
    else if(rstrat.method==REDUCTION_GLOVER_SDCE){
        printf("simple diversity-based clustering.\n");
        reduction_diversity_starting(prob,sols,n_sols,rstrat.n_target,
            rstrat.soldis,rstrat.facdis,0);
    }
    else if(rstrat.method==REDUCTION_GLOVER_SDCE_BESTS){
        printf("simple diversity-based clustering (bests).\n");
        reduction_diversity_starting(prob,sols,n_sols,rstrat.n_target,
            rstrat.soldis,rstrat.facdis,1);
    }
    else{
        printf("???.\n");
        fprintf(stderr,"ERROR: Reduction method not yet implemented.\n");
        exit(1);
    }
}

void reduction_bests(const problem *prob, solution **sols, int *n_sols, int n_target){
    // Sort solutions in decreasing order
    qsort(sols,*n_sols,sizeof(solution *),solutionp_value_cmp_inv);
    // If there are already less solutions, do nothing.
    if(*n_sols<=n_target) return;
    // Free other solutions:
    for(int i=n_target; i<*n_sols; i++){
        solution_free(sols[i]);
    }
    // Set the amount of solutions right.
    *n_sols = n_target;
}

void reduction_random_uniform(const problem *prob, solution **sols, int *n_sols, int n_target){
    if(*n_sols<=n_target) return;
    // Put target_n randomly selected solutions first on the array, but keep the best so far.
    for(int i=1;i<n_target;i++){ // Fisher-Yates shuffle
        int choice = i+rand()%(*n_sols-i);
        solution *aux = sols[i];
        sols[i] = sols[choice];
        sols[choice] = aux;
    }
    // Free other solutions:
    for(int i=n_target;i<*n_sols;i++){
        solution_free(sols[i]);
    }
    // Set the amount of solutions right.
    *n_sols = n_target;
}

