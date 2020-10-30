#include "reduction.h"

void reduce_by_redstrategy(const rundata *run, const redstrategy rstrat,
        solution **sols, int *n_sols){
    if(*n_sols<=rstrat.n_target) return;

    if(run->verbose) printf(
            "Reducing \033[31;1m%d\033[0m -> \033[31;1m%d\033[0m solutions, ",*n_sols,rstrat.n_target);
    if(rstrat.method==REDUCTION_BESTS){
        if(run->verbose) printf("selecting bests.\n");
        reduction_bests(run,sols,n_sols,rstrat.n_target);
    }
    else if(rstrat.method==REDUCTION_RANDOM_UNIFORM){
        if(run->verbose){
            if(rstrat.elitist) printf("randomly (uniform, elitist).\n");
            else printf("randomly (uniform).\n");
        }
        reduction_random_uniform(run,sols,n_sols,rstrat.n_target,rstrat.elitist);
    }
    else if(rstrat.method==REDUCTION_RANDOM_RANK){
        if(run->verbose){
            if(rstrat.elitist) printf("randomly (by rank, elitist).\n");
            else printf("randomly (by rank).\n");
        }
        reduction_random_rank(run,sols,n_sols,rstrat.n_target,rstrat.elitist);
    }
    else if(rstrat.method==REDUCTION_GLOVER_SDBS){
        if(run->verbose) printf("simple diversity-based starting method.\n");
        reduction_diversity_starting(run,sols,n_sols,rstrat.n_target,
            rstrat.soldis,rstrat.facdis,0);
    }
    else if(rstrat.method==REDUCTION_GLOVER_SDBS_BESTS){
        if(run->verbose) printf("simple diversity-based starting method (bests-of-cluster).\n");
        reduction_diversity_starting(run,sols,n_sols,rstrat.n_target,
            rstrat.soldis,rstrat.facdis,1);
    }
    else if(rstrat.method==REDUCTION_VRHEURISTIC){
        if(run->verbose) printf("VR heuristic (vision range: %d).\n",rstrat.arg);
        reduction_vr_heuristic(run,sols,n_sols,rstrat.n_target,
            rstrat.soldis,rstrat.facdis,rstrat.arg);
    }
    else{
        if(run->verbose) printf("???.\n");
        fprintf(stderr,"ERROR: Reduction method not yet implemented.\n");
        exit(1);
    }
}

void reduction_bests(const rundata *run, solution **sols, int *n_sols, int n_target){
    // Sort solutions in decreasing order
    qsort(sols,*n_sols,sizeof(solution *),solutionp_value_cmp_inv);
    assert(*n_sols<2 || sols[0]->value>=sols[1]->value);
    // If there are already less solutions, do nothing.
    if(*n_sols<=n_target) return;
    // Free other solutions:
    for(int i=n_target; i<*n_sols; i++){
        solution_free(sols[i]);
    }
    // Set the amount of solutions right.
    *n_sols = n_target;
}

void reduction_random_uniform(const rundata *run, solution **sols, int *n_sols, int n_target, int elitist){
    if(*n_sols<=n_target) return;
    // Put target_n randomly selected solutions first on the array, but keep the best so far if elitist.
    for(int i=elitist;i<n_target;i++){ // Fisher-Yates shuffle
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

void reduction_remove_worst(solution **sols, int *n_sols, double factor){
    // Sort solutions in decreasing value
    qsort(sols,*n_sols,sizeof(solution *),solutionp_value_cmp_inv);
    // How many solutions to delete?
    int n_delete = (int) floorf((*n_sols)*factor);
    assert(n_delete <= *n_sols);
    // Delete last solutions in the array
    for(int i=0;i<n_delete;i++){
        int idx = *n_sols -1 - i;
        assert(idx>=0);
        solution_free(sols[idx]);
    }
    *n_sols = *n_sols - n_delete;
}
