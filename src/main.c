#include "load.h"
#include "problem.h"
#include "expand.h"
#include "reduction.h"
#include "construction.h"

int main(int argc, const char **argv){
    // Print information if arguments are invalid
    if(argc<3){
        fprintf(stderr,"usage: %s {strategy:n} <input_file> <output_file>\n",argv[0]);
        exit(1);
    }
    
    int n_strategies = argc-3;
    redstrategy *strategies = redstrategy_init_from_nomenclatures(&argv[1],&n_strategies);

    // Read problem
    problem *prob = new_problem_load(argv[argc-2]);
    
    // Print current run info:
    problem_print(prob,stdout);
    printf("# REDUCTION:");
    for(int i=0;i<n_strategies;i++) printf(" %s",strategies[i].nomenclature);
    printf("\n");

    // Final solutions
    int final_n_sols;
    int n_iterations;
    solution **final_sols = new_find_best_solutions(prob,strategies,n_strategies,
        &final_n_sols, &n_iterations);

    // Print solutions
    for(int i=0;i<final_n_sols;i++){
        solution_print(prob,final_sols[i],stdout);
    }
    // Free memory
    for(int i=0;i<final_n_sols;i++){
        solution_free(final_sols[i]);
    }
    free(final_sols);
    problem_free(prob);
    free(strategies);
}