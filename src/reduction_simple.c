#include "reduction.h"

void reduction_bests(const problem *prob, solution **sols, int *n_sols, int n_target){
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
    // If there are already less solutions, do nothing.
    if(*n_sols<=n_target) return;
    // Put target_n randomly selected solutions first on the array, but keep the best so far.
    for(int i=1;i<n_target;i++){
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

void reduction_random_rank(const problem *prob, solution **sols, int *n_sols, int n_target){
}