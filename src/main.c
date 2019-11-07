#include "load.h"
#include "problem.h"
#include "expand.h"
#include "reduction.h"

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

    // Free memory
    problem_free(prob);
    free(strategies);
}