#include "load.h"
#include "problem.h"
#include "expand.h"
#include "reduction.h"
#include "construction.h"
#include "output.h"

#include <time.h>
#include <sys/time.h>

int main(int argc, const char **argv){
    // Print information if arguments are invalid
    if(argc<3){
        fprintf(stderr,"usage: %s [-r<n>] [-n<n>] [-t<n>] [-f<n>] [-s<n>] [-S<n>] [-b] [-l] {strategy:n} <input> <output>\n",argv[0]);
        exit(1);
    }

    // Strategy nomenclature arguments
    int n_strategies = 0;
    const char **strategy_args = safe_malloc(sizeof(const char *)*argc);

    // Filenames
    const char *input_fname = argv[argc-2];
    const char *output_fname = argv[argc-1];

    // Problem arguments to be changed by command line arguments
    int random_seed = -1;
    int target_n = -1;
    int filter_n = -1;
    int max_size = -1;
    int min_size = -1;
    int n_threads = -1;
    int local_search = -1;
    int bnb = -1;

    // Parse arguments
    for(int i=1;i<argc-2;i++){
        if(argv[i][0]=='-'){
            if(argv[i][1]=='r'){
                // Number of target solutions
                int n_read = sscanf(argv[i],"-r%d",&random_seed);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected seed on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
            }else if(argv[i][1]=='n'){
                // Number of target solutions
                int n_read = sscanf(argv[i],"-n%d",&target_n);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected number of target sols. on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
            }else if(argv[i][1]=='t'){
                // Number of target solutions
                int n_read = sscanf(argv[i],"-t%d",&n_threads);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected number of threads on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
            }else if(argv[i][1]=='s'){
                // Minimum solution size
                int n_read = sscanf(argv[i],"-s%d",&min_size);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected minimum size on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
            }else if(argv[i][1]=='S'){
                // Number of target solutions
                int n_read = sscanf(argv[i],"-S%d",&max_size);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected maximum size on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
            }else if(argv[i][1]=='f'){
                // Select kind of filter
                int n_read = sscanf(argv[i],"-f%d",&filter_n);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected filter level on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
                assert(0<=filter_n && filter_n<=MAX_FILTER);
            }else if(argv[i][1]=='b' && strcmp(argv[i],"-b")==0){
                // Disable branch and bound
                bnb = 0;
            }else if(argv[i][1]=='l' && strcmp(argv[i],"-l")==0){
                // Disable local search
                local_search = 0;
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
    problem *prob = new_problem_load(input_fname);
    if(target_n>=0) prob->target_sols = target_n;
    if(filter_n>=0) prob->filter = filter_n;
    if(bnb>=0) prob->branch_and_bound = bnb;
    if(min_size>=0) prob->size_restriction_minimum = min_size;
    if(max_size>=0) prob->size_restriction_maximum = max_size;
    if(n_threads>0) prob->n_threads = n_threads;
    if(local_search>0) prob->local_search = local_search; 
    
    // Set random seed
    if(random_seed==-1) random_seed = (int) time(NULL);
    prob->random_seed = random_seed; 

    // Start counting time (cpu and elapsed)
    clock_t start = clock();
    struct timeval elapsed_start;
    gettimeofday(&elapsed_start,NULL);
    // ---@>
    
    // Perform precomputations
    // (derivated attributes of the problem, like facility-facility distanceS)
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

    // End counting time
    clock_t end = clock();
    struct timeval elapsed_end;
    gettimeofday(&elapsed_end,NULL);
    float seconds = (float)(end - start) / (float)CLOCKS_PER_SEC;
    float elapsed_seconds = get_delta_seconds(elapsed_start,elapsed_end);
    // ---@>

    // Save output
    save_solutions(output_fname,
        prob,final_sols,final_n_sols,input_fname,
        seconds,elapsed_seconds,n_iterations,strategies,n_strategies);

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