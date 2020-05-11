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
        fprintf(stderr,"usage: %s [-r<n>] [-n<n>] [-R<n>] [-B<n>] [-t<n>] [-f<n>] [-s<n>] [-S<n>] [-b] [-l] [-A] [-x] {strategy:n} <input> <output>\n",argv[0]);
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
    int restarts = -1;
    int target_n = -1;
    int filter_n = -1;
    int max_size = -1;
    int min_size = -1;
    int n_threads = -1;
    int local_search = -1;
    int local_search_only_on_terminal = -1;
    int local_search_rem_movement = -1;
    int local_search_add_movement = -1;
    int bnb = -1;
    int verbose = -1;
    int branching = -2;

    // Parse arguments
    for(int i=1;i<argc-2;i++){
        if(argv[i][0]=='-'){
            if(argv[i][1]=='r'){
                // Random seed
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
            }else if(argv[i][1]=='R'){
                // Number of restarts
                int n_read = sscanf(argv[i],"-R%d",&restarts);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected number of restarts on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
            }else if(argv[i][1]=='B'){
                // Set branching factor
                int n_read = sscanf(argv[i],"-B%d",&branching);
                if(n_read<1){
                   fprintf(stderr,"ERROR: expected branching factor on argument \"%s\".\n",argv[i]);
                   exit(1);
                }
                assert(branching>=-2);
            }else if(argv[i][1]=='t'){
                // Number of threads
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
                local_search = NO_LOCAL_SEARCH;
            }else if(argv[i][1]=='L' && strcmp(argv[i],"-L")==0){
                // First improvement local search
                local_search = SWAP_FIRST_IMPROVEMENT;
            }else if(argv[i][1]=='A' && strcmp(argv[i],"-A")==0){
                // Perform local search on all nodes
                local_search_only_on_terminal = 0;
            }else if(argv[i][1]=='x' && strcmp(argv[i],"-x")==0){
                // Don't allow local search to perform movements that change the size
                local_search_rem_movement = 0;
                local_search_add_movement = 0;
            }else if(argv[i][1]=='V' && strcmp(argv[i],"-V")==0){
                // Non verbose mode
                verbose = 0;
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

    if(restarts<0) restarts = 1;

    // Read problem and set size restrictions
    problem *prob = new_problem_load(input_fname);
    if(min_size>=0) prob->size_restriction_minimum = min_size;
    if(max_size>=0) prob->size_restriction_maximum = max_size;

    // Create rundata (perform precomputations)
    if(verbose!=0) printf("\nPerforming precomputations.\n");
    rundata *run = rundata_init(prob, strategies,n_strategies,restarts);

    // Free problem (rundata kepps a copy)
    problem_free(prob);
    prob = NULL;

    // Set console arguments
    if(target_n>=0) run->target_sols = target_n;
    if(filter_n>=0) run->filter = filter_n;
    if(bnb>=0) run->branch_and_bound = bnb;
    if(n_threads>0) run->n_threads = n_threads;
    if(local_search>=0) run->local_search = local_search;
    if(local_search_only_on_terminal>=0) run->local_search_only_terminal = local_search_only_on_terminal;
    if(local_search_rem_movement>=0) run->local_search_rem_movement = local_search_rem_movement;
    if(local_search_add_movement>=0) run->local_search_add_movement = local_search_add_movement;
    if(verbose>=0) run->verbose = verbose;
    if(branching>-2) run->branching_factor = branching;

    // Set random seed
    if(random_seed==-1) random_seed = (int) time(NULL);
    run->random_seed = random_seed;

    // Start counting time (cpu and elapsed)
    clock_t start = clock();
    struct timeval elapsed_start;
    gettimeofday(&elapsed_start,NULL);
    // ---@>

    // Print current run info:
    if(run->verbose){
        rundata_print(run,stdout);
        printf("# REDUCTION:");
        for(int i=0;i<n_strategies;i++) printf(" %s",strategies[i].nomenclature);
        printf("\n");
        printf("\n");
    }

    // Final solutions
    int final_n_sols;
    solution **final_sols = new_find_best_solutions(run,strategies,n_strategies,
        &final_n_sols);

    // End counting time
    clock_t end = clock();
    struct timeval elapsed_end;
    gettimeofday(&elapsed_end,NULL);
    float seconds = (float)(end - start) / (float)CLOCKS_PER_SEC;
    float elapsed_seconds = get_delta_seconds(elapsed_start,elapsed_end);
    // ---@>

    if(run->verbose) printf("\n");

    int virt_mem_usage_peak;
    get_memory_usage(NULL,NULL,NULL,&virt_mem_usage_peak);

    // Save output
    save_solutions(output_fname,
        run,final_sols,final_n_sols,input_fname,
        seconds,elapsed_seconds,virt_mem_usage_peak,
        strategies,n_strategies);

    // Print output
    save_solutions(NULL,
        run,final_sols,final_n_sols,input_fname,
        seconds,elapsed_seconds,virt_mem_usage_peak,
        strategies,n_strategies);

    // Free memory
    for(int i=0;i<final_n_sols;i++){
        solution_free(final_sols[i]);
    }
    rundata_free(run);
    free(final_sols);
    free(strategies);
    free(strategy_args);
}
