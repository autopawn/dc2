#include "load.h"
#include "problem.h"
#include "expand.h"
#include "reduction.h"
#include "construction.h"
#include "output.h"

#include <time.h>
#include <sys/time.h>

#define VALUE_NOT_SET

int main(int argc, const char **argv){
    // Print information if arguments are invalid
    if(argc<3){
        fprintf(stderr,"usage: %s [OPTIONS] {strategy:n} <input> <output>\n",argv[0]);
        exit(1);
    }

    // Strategy nomenclature arguments
    int n_strategies = 0;
    const char **strategy_args = safe_malloc(sizeof(const char *)*argc);

    // Filenames
    const char *input_fname = argv[argc-2];
    const char *output_fname = argv[argc-1];

    // Problem arguments to be changed by command line arguments
    const int UNSET = -2;
    int random_seed = UNSET;
    int restarts = UNSET;
    int target_n = UNSET;
    int filter_n = UNSET;
    int max_size = UNSET;
    int min_size = UNSET;
    int n_threads = UNSET;
    int local_search = UNSET;
    int local_search_only_on_terminal = UNSET;
    int local_search_rem_movement = UNSET;
    int local_search_add_movement = UNSET;
    int bnb = UNSET;
    int verbose = UNSET;
    int branching = UNSET;
    int path_relinking = UNSET;

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
                assert(local_search==UNSET);
                local_search = NO_LOCAL_SEARCH;
            }else if(argv[i][1]=='L' && strcmp(argv[i],"-L")==0){
                // First improvement local search
                assert(local_search==UNSET);
                local_search = SWAP_FIRST_IMPROVEMENT;
            }else if(argv[i][1]=='W' && strcmp(argv[i],"-W")==0){
                // Resende and Werneck's local search
                assert(local_search==UNSET);
                local_search = SWAP_RESENDE_WERNECK;
            }else if(argv[i][1]=='A' && strcmp(argv[i],"-A")==0){
                // Perform local search on all nodes
                local_search_only_on_terminal = 0;
            }else if(argv[i][1]=='x' && strcmp(argv[i],"-x")==0){
                // Don't allow local search to perform movements that change the size
                local_search_rem_movement = 0;
                local_search_add_movement = 0;
            }else if(argv[i][1]=='P' && strcmp(argv[i],"-P")==0){
                // Enable 1 path relinking
                assert(path_relinking==UNSET);
                path_relinking = PATH_RELINKING_1_ITER;
            }else if(argv[i][1]=='M' && strcmp(argv[i],"-M")==0){
                // Enable path relinking until no better solution is found made
                assert(path_relinking==UNSET);
                path_relinking = PATH_RELINKING_UNTIL_NO_BETTER;
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
    int precomp_nearly_indexes = (local_search==SWAP_RESENDE_WERNECK);
    rundata *run = rundata_init(prob, strategies,n_strategies,restarts,precomp_nearly_indexes,n_threads,verbose);

    // Free problem (rundata kepps a copy)
    problem_free(prob);
    prob = NULL;

    // Set console arguments
    if(target_n!=UNSET) run->target_sols = target_n;
    if(filter_n!=UNSET) run->filter = filter_n;
    if(bnb!=UNSET) run->branch_and_bound = bnb;
    if(n_threads>0) run->n_threads = n_threads;
    if(local_search!=UNSET) run->local_search = local_search;
    if(local_search_only_on_terminal!=UNSET) run->local_search_only_terminal = local_search_only_on_terminal;
    if(local_search_rem_movement!=UNSET) run->local_search_rem_movement = local_search_rem_movement;
    if(local_search_add_movement!=UNSET) run->local_search_add_movement = local_search_add_movement;
    if(verbose!=UNSET) run->verbose = verbose;
    if(branching!=UNSET) run->branching_factor = branching;
    if(path_relinking!=UNSET) run->path_relinking = path_relinking;

    // Check that there is no reduction for path relinking if we are not using it
    if(run->path_relinking==NO_PATH_RELINKING){
        for(int i=0;i<n_strategies;i++){
            if(strategies[i].for_path_relinking){
                fprintf(stderr,"ERROR: \"%s\" reduction but no path relinking.\n",strategies[i].nomenclature);
                exit(1);
            }
        }
    }

    // Simple coeficients version of dc2 check:
    #ifdef DC2_110
       for(int c=0;c<run->prob->n_clis;c++){
           if(run->prob->client_weight[c]!=1){
               fprintf(stderr,"ERROR: %s expects problems with client_weight[%d] = 1\n",argv[0],c);
	       exit(1);
	   }
       }
       if(run->prob->transport_cost!=1){
	       fprintf(stderr,"ERROR: %s expects problems with transport_cost = 1\n",argv[0]);
	       exit(1);
       }
       if(run->prob->client_gain!=0){
	       fprintf(stderr,"ERROR: %s expects problems with client_gain = 0\n",argv[0]);
	       exit(1);
       }
    #endif

    // Set random seed
    if(random_seed==UNSET) random_seed = (int) time(NULL);
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
