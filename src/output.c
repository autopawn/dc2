#include "output.h"

float get_delta_seconds(struct timeval tv1, struct timeval tv2){
    struct timeval tvdiff = {tv2.tv_sec-tv1.tv_sec,tv2.tv_usec-tv1.tv_usec};
    if(tvdiff.tv_usec<0){
        tvdiff.tv_usec += 1000000;
        tvdiff.tv_sec -= 1;
    }
    float secs = (float)(tvdiff.tv_sec)+1e-6*tvdiff.tv_usec;
    return secs;
}

void save_solutions(const char *file,
        const rundata *run, solution **sols, int n_sols,
        const char *input_file, float seconds, float elapsed, int mem_usage,
        const redstrategy *strategies, int n_strategies, int only_1_output_sol){
    FILE *fp;
    // Open output file
    if(file){
        printf("Opening file \"%s\" for saving results...\n",file);
        fp = fopen(file,"w");
        if(fp==NULL){
            fprintf(stderr,"ERROR: couldn't open file \"%s\"!\n",file);
            exit(1);
        }
        printf("\n");
    }else{
        fp = stdout;
    }
    // Print problem data
    rundata_print(run,fp);
    // Print reduction scheme
    fprintf(fp,"\n");
    fprintf(fp,"# REDUCTION:");
    for(int i=0;i<n_strategies;i++) fprintf(fp," %s",strategies[i].nomenclature);
    fprintf(fp,"\n");
    // Print other useful information
    fprintf(fp,"# INPUT_FILE: \"%s\"\n",input_file);
    fprintf(fp,"# CPU_TIME: %f\n",seconds);
    fprintf(fp,"# ELAPSED: %f\n",elapsed);
    fprintf(fp,"# VIRT_MEM_PEAK_KB: %d\n",mem_usage);
    fprintf(fp,"# TOTAL_ITERATIONS: %d\n",run->run_inf->total_n_iterations);
    fprintf(fp,"\n");

    /* LOCAL SEARCH INFO */
    fprintf(fp,"== LOCAL SEARCH INFO ==\n");
    fprintf(fp,"# LOCAL_SEARCH_CPU_TIME: %f\n",run->run_inf->local_search_seconds);
    fprintf(fp,"# N_LOCAL_SEARCHES: %lld\n",run->run_inf->n_local_searches);
    fprintf(fp,"# AVG_LOCAL_SEARCH_MOVES: %f\n",(double)run->run_inf->n_local_search_movements/(double)run->run_inf->n_local_searches);
    fprintf(fp,"\n");

    /* PATH RELINKING INFO */
    fprintf(fp,"== PATH RELINKING INFO ==\n");
    fprintf(fp,"# PATH_RELINKING_CPU_TIME: %f\n",run->run_inf->path_relinking_seconds);

    /* FIRST RESTART DATA */
    fprintf(fp,"== FIRST RESTART INFO ==\n");
    fprintf(fp,"# FIRST_ITERATIONS: %d\n",run->run_inf->firstr_n_iterations);
    //
    fprintf(fp,"# FIRST_PER_SIZE_SOLS:       ");
    for(int i=0;i<run->run_inf->firstr_n_iterations;i++){
        fprintf(fp," %6d",run->run_inf->firstr_per_size_n_sols[i]);
    }
    fprintf(fp,"\n");
    //
    fprintf(fp,"# FIRST_PER_SIZE_SOLS_RED:   ");
    for(int i=0;i<run->run_inf->firstr_n_iterations;i++){
        fprintf(fp," %6d",run->run_inf->firstr_per_size_n_sols_after_red[i]);
    }
    fprintf(fp,"\n");
    //
    fprintf(fp,"# FIRST_PER_SIZE_LOCAL_OPT:  ");
    for(int i=0;i<run->run_inf->firstr_n_iterations;i++){
        fprintf(fp," %6d",run->run_inf->firstr_per_size_n_local_optima[i]);
    }
    fprintf(fp,"\n");
    //
    fprintf(fp,"\n");

    /* EACH RESTART DATA */
    fprintf(fp,"== PER RESTART INFO ==\n");
    for(int i=0;i<run->n_restarts;i++){
        fprintf(fp,"# RST %8d %lf %lf\n",i,run->run_inf->restart_values[i],run->run_inf->restart_times[i]);
    }
    fprintf(fp,"\n");


    /* SOLUTIONS DATA */
    fprintf(fp,"== SOLUTIONS DATA ==\n");
    fprintf(fp,"# OUTPUT_SOLUTIONS: %d\n",n_sols);
    for(int i=0;i<n_sols;i++){
        if(i==0 || !only_1_output_sol){
            fprintf(fp,"\n");
            solution_print(run->prob,sols[i],fp);
        }
    }
    // Close file descriptor
    if(fp!=stdout){
        fclose(fp);
    }
}
