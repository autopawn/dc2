#include "load.h"
#include "redstrategy.h"
#include "problem.h"
#include "solution.h"

/*
The opt_checker is a simple program that reads the facilitites
present in the solution on a .opt or .bub file and uses those
indexes to build its own solution and retrieve it in the same
format.
*/

#include <time.h>
#include <sys/time.h>

int main(int argc, const char **argv){
    // Print information if arguments are invalid
    if(argc!=3){
        fprintf(stderr,"usage: %s <input> <input_opt>\n",argv[0]);
        exit(1);
    }

    const char *input_fname = argv[1];
    const char *opt_fname = argv[2];

    // Read problem
    problem *prob = new_problem_load(input_fname);

    // Perform precomputations
    problem_precompute(prob,NULL,0);

    // Create empty solution
    solution *solution = solution_empty(prob);

    // Read opt file
    FILE *fp = fopen(opt_fname,"r");
    assert(fp!=NULL);
    for(int i=0;i<prob->n_clis;i++){
        int fac;
        int n_read = fscanf(fp,"%d",&fac);
        assert(n_read==1);
        solution_add(prob,solution,fac);
        // Get the assign cost by the assignment made by the algorithm
        double old_assign_cost = problem_assig_value(prob,solution->assigns[i],i);
        // Get the assign cost by the assignment on the .opt solution
        double opt_assign_cost = problem_assig_value(prob,fac,i);
        // Assert that both are the same
        assert(old_assign_cost==opt_assign_cost);
        // Replace algorithm assignment with .opt assignment
        solution->assigns[i] = fac;
    }
    fclose(fp);

    // Print result
    for(int i=0;i<prob->n_clis;i++){
        printf("%d ",solution->assigns[i]);
    }
    printf("%.9lf\n",-solution->value);

    // Free memory
    solution_free(solution);
    problem_free(prob);
}