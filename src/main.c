#include "load.h"
#include "problem.h"
#include "expand.h"
#include "reduction.h"
#include "construction.h"

int main(int argc, const char **argv){
    // Print information if arguments are invalid
    if(argc<3){
        fprintf(stderr,"usage: %s [-n<target_n>] [-f<filter>] {strategy:n} <input> <output>\n",argv[0]);
        exit(1);
    }

    // Set a fixed seed
    srand(0xc0ffee);

    // Strategy nomenclature arguments
    int n_strategies = 0;
    const char **strategy_args = safe_malloc(sizeof(const char *)*argc);

    // Target number of solutions
    int target_n = -1;
    int filter_n = -1;

    // Parse arguments
    for(int i=1;i<argc-2;i++){
        if(argv[i][0]=='-'){
            if(argv[i][1]=='n'){
                int n_read = sscanf(argv[i],"-n%d",&target_n);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected number of target sols. on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
            }else if(argv[i][1]=='f'){
                int n_read = sscanf(argv[i],"-f%d",&filter_n);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected filter level on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
                assert(0<=filter_n && filter_n<=MAX_FILTER);
            }else{
                fprintf(stderr,"ERROR: argument \"%s\" not recognized.\n",argv[i]);
                exit(1);
            }
        }else{
            strategy_args[n_strategies] = argv[i];
            n_strategies++;
        }
    }

    redstrategy *strategies = redstrategy_init_from_nomenclatures(strategy_args,&n_strategies);

    // Read problem and use console arguments
    problem *prob = new_problem_load(argv[argc-2]);
    if(target_n>=0) prob->target_sols = target_n;
    if(filter_n>=0) prob->filter = filter_n;
    
    printf("\n");
    printf("Performing precomputations.\n");
    problem_precompute(prob,strategies,n_strategies);
    
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
    free(strategy_args);
}