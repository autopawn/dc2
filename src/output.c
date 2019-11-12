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
        const problem *prob, solution **sols, int n_sols,
        const char *input_file, float seconds, float elapsed, int n_iters,
        const redstrategy *strategies, int n_strategies){
    FILE *fp;
    // Open output file
    printf("Opening file \"%s\" for saving results...\n",file);
    fp = fopen(file,"w");
    if(fp==NULL){
        fprintf(stderr,"ERROR: couldn't open file \"%s\"!\n",file);
        exit(1);
    }
    // Print problem data
    problem_print(prob,fp);
    // Print reduction scheme
    fprintf(fp,"\n");
    fprintf(fp,"# REDUCTION:");
    for(int i=0;i<n_strategies;i++) fprintf(fp," %s",strategies[i].nomenclature);
    fprintf(fp,"\n");
    // Print other useful information
    fprintf(fp,"# CPU_TIME: %f\n",seconds);
    fprintf(fp,"# ELAPSED: %f\n",elapsed);
    fprintf(fp,"# INPUT_FILE: \"%s\"\n",input_file);
    fprintf(fp,"# ITERATIONS: %d\n",n_iters);
    fprintf(fp,"\n");
    fprintf(fp,"# SOLUTIONS: %d\n",n_sols);
    for(int i=0;i<n_sols;i++){
        fprintf(fp,"\n");
        solution_print(prob,sols[i],fp);
    }
    fclose(fp);
}